#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <vector>

// Macros de utilidad
#define min(A, B) ((A) < (B) ? (A) : (B))
#define max(A, B) ((A) > (B) ? (A) : (B))

// Constantes del problema (idénticas en ambas versiones)
const int imax = 80;
const int kmax = 80;
const int itmax = 20000;
const double eps = 1.0e-08;

// Tags MPI
#define TAG_UP   100
#define TAG_DOWN 101

// Estructura para reportar tiempos MPI
struct TiemposMPI {
    double total;
    double computo;
    double comunicacion;
};

// ======================================================
// FUNCIÓN 1: SOLUCIONADOR MPI (Tu version2.cpp adaptada)
// ======================================================
TiemposMPI ejecutar_mpi(int rank, int size) {
    double dx, dy, dx2, dy2, dx2i, dy2i, dt;
    int k, it;

    // Parámetros físicos
    dx = 1.0 / kmax;
    dy = 1.0 / imax;
    dx2 = dx * dx;
    dy2 = dy * dy;
    dx2i = 1.0 / dx2;
    dy2i = 1.0 / dy2;
    dt = min(dx2, dy2) / 4.0;

    // Descomposición de dominio
    int filas_totales = imax - 1;
    int filas_base = filas_totales / size;
    int filas_extra = filas_totales % size;
    int filas_locales;
    int fila_inicio_global;

    if (rank < filas_extra) {
        filas_locales = filas_base + 1;
        fila_inicio_global = rank * filas_locales + 1;
    } else {
        filas_locales = filas_base;
        fila_inicio_global = rank * filas_base + filas_extra + 1;
    }

    // Allocación
    double** phi_local = new double*[filas_locales + 2];
    double** phin_local = new double*[filas_locales + 2];
    for (int i = 0; i < filas_locales + 2; i++) {
        phi_local[i] = new double[kmax + 1];
        phin_local[i] = new double[kmax + 1];
    }

    // Inicialización
    for (int i = 0; i < filas_locales + 2; i++) {
        for (k = 0; k <= kmax; k++) {
            phi_local[i][k] = 0.0;
            phin_local[i][k] = 0.0;
        }
        phi_local[i][kmax] = 1.0; // Frontera superior
    }
    
    // Fronteras laterales
    if (rank == 0) phi_local[0][0] = 0.0;
    if (rank == size - 1) phi_local[filas_locales + 1][0] = 0.0;

    if (rank == 0) {
        for (k = 1; k < kmax; k++) phi_local[0][k] = phi_local[0][k - 1] + dx;
    }
    if (rank == size - 1) {
        for (k = 1; k < kmax; k++) phi_local[filas_locales + 1][k] = phi_local[filas_locales + 1][k - 1] + dx;
    }

    // --- INICIO MEDICIÓN MPI ---
    MPI_Barrier(MPI_COMM_WORLD); // Sincronizar antes de medir
    double t_inicio = MPI_Wtime();
    double t_comunicacion_acum = 0.0;
    double t_aux_start;

    for (it = 1; it <= itmax; it++) {
        double dphimax_local = 0.0;
        double dphimax_global = 0.0;

        MPI_Request requests[4];
        int req_count = 0;

        // >> MEDIR COMUNICACIÓN (ISEND/IRECV)
        t_aux_start = MPI_Wtime();
        if (rank > 0 && filas_locales > 0)
            MPI_Isend(phi_local[1], kmax + 1, MPI_DOUBLE, rank - 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
        if (rank < size - 1 && filas_locales > 0)
            MPI_Isend(phi_local[filas_locales], kmax + 1, MPI_DOUBLE, rank + 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
        if (rank > 0 && filas_locales > 0)
            MPI_Irecv(phi_local[0], kmax + 1, MPI_DOUBLE, rank - 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
        if (rank < size - 1 && filas_locales > 0)
            MPI_Irecv(phi_local[filas_locales + 1], kmax + 1, MPI_DOUBLE, rank + 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
        t_comunicacion_acum += (MPI_Wtime() - t_aux_start);

        // >> CÁLCULO (INTERIOR)
        int i_inicio = 2, i_fin = filas_locales - 1;
        if (i_fin >= i_inicio) {
            for (k = 1; k < kmax; k++) {
                for (int i_local = i_inicio; i_local <= i_fin; i_local++) {
                    double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k] - 2.0 * phi_local[i_local][k]) * dy2i +
                                  (phi_local[i_local][k + 1] + phi_local[i_local][k - 1] - 2.0 * phi_local[i_local][k]) * dx2i;
                    dphi *= dt;
                    dphimax_local = max(dphimax_local, std::fabs(dphi));
                    phin_local[i_local][k] = phi_local[i_local][k] + dphi;
                }
            }
        }

        // >> MEDIR COMUNICACIÓN (WAITALL)
        t_aux_start = MPI_Wtime();
        if (req_count > 0) MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
        t_comunicacion_acum += (MPI_Wtime() - t_aux_start);

        // >> CÁLCULO (BORDES HALO)
        if (filas_locales >= 1) {
            int i_local = 1;
            for (k = 1; k < kmax; k++) {
                double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k] - 2.0 * phi_local[i_local][k]) * dy2i +
                              (phi_local[i_local][k + 1] + phi_local[i_local][k - 1] - 2.0 * phi_local[i_local][k]) * dx2i;
                dphi *= dt;
                dphimax_local = max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }
        if (filas_locales >= 2) {
            int i_local = filas_locales;
            for (k = 1; k < kmax; k++) {
                double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k] - 2.0 * phi_local[i_local][k]) * dy2i +
                              (phi_local[i_local][k + 1] + phi_local[i_local][k - 1] - 2.0 * phi_local[i_local][k]) * dx2i;
                dphi *= dt;
                dphimax_local = max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }

        // Actualización
        for (int i_local = 1; i_local <= filas_locales; i_local++) {
            for (k = 1; k < kmax; k++) phi_local[i_local][k] = phin_local[i_local][k];
        }

        // >> MEDIR COMUNICACIÓN (ALLREDUCE)
        t_aux_start = MPI_Wtime();
        MPI_Allreduce(&dphimax_local, &dphimax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        t_comunicacion_acum += (MPI_Wtime() - t_aux_start);

        if (dphimax_global < eps) break;
    }
    double t_final = MPI_Wtime();

    // Limpieza
    for (int i = 0; i < filas_locales + 2; i++) {
        delete[] phi_local[i];
        delete[] phin_local[i];
    }
    delete[] phi_local;
    delete[] phin_local;

    // Retornar resultados
    TiemposMPI res;
    res.total = t_final - t_inicio;
    res.comunicacion = t_comunicacion_acum;
    res.computo = res.total - res.comunicacion;
    return res;
}

// ======================================================
// FUNCIÓN 2: SOLUCIONADOR OPENMP (Tu version1.cpp adaptada)
// ======================================================
double ejecutar_omp() {
    double phi[imax + 1][kmax + 1], phin[imax][kmax];
    double dx, dy, dx2, dy2, dx2i, dy2i, dt, dphi, dphimax;
    int i, k, it;

    dx = 1.0 / kmax;
    dy = 1.0 / imax;
    dx2 = dx * dx;
    dy2 = dy * dy;
    dx2i = 1.0 / dx2;
    dy2i = 1.0 / dy2;
    dt = min(dx2, dy2) / 4.0;
    
    // Inicialización
    for (k = 0; k < kmax; k++) for (i = 1; i < imax; i++) phi[i][k] = 0.0;
    for (i = 0; i <= imax; i++) phi[i][kmax] = 1.0;
    phi[0][0] = 0.0; phi[imax][0] = 0.0;
    for (k = 1; k < kmax; k++) {
        phi[0][k] = phi[0][k - 1] + dx;
        phi[imax][k] = phi[imax][k - 1] + dx;
    }

    double t_inicio = omp_get_wtime(); // Tiempo OMP

    for (it = 1; it <= itmax; it++) {
        dphimax = 0.;
        
        #pragma omp parallel for private(i, dphi) reduction(max:dphimax) schedule(static)
        for (k = 1; k < kmax; k++) {
            for (i = 1; i < imax; i++) {
                dphi = (phi[i + 1][k] + phi[i - 1][k] - 2. * phi[i][k]) * dy2i + 
                       (phi[i][k + 1] + phi[i][k - 1] - 2. * phi[i][k]) * dx2i;
                dphi = dphi * dt;
                dphimax = max(dphimax, dphi);
                phin[i][k] = phi[i][k] + dphi;
            }
        }

        #pragma omp parallel for private(i) schedule(static)
        for (k = 1; k < kmax; k++) {
            for (i = 1; i < imax; i++) {
                phi[i][k] = phin[i][k];
            }
        }
        
        if (dphimax < eps) break;
    }

    double t_fin = omp_get_wtime();
    return t_fin - t_inicio;
}

// ======================================================
// MAIN: CONTROLADOR PRINCIPAL
// ======================================================
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv); // Inicializar MPI

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("\n========================================\n");
        printf("  BENCHMARK: MPI vs OpenMP\n");
        printf("  Grid: %dx%d | Iteraciones Max: %d\n", imax, kmax, itmax);
        printf("========================================\n\n");
    }

    // 1. Ejecutar MPI en todos los nodos
    if (rank == 0) printf("[MPI] Ejecutando con %d procesos...\n", size);
    TiemposMPI tiempos_mpi = ejecutar_mpi(rank, size);

    // Recolectar el tiempo máximo de todos los procesos MPI (para ser justos, tomamos el más lento)
    double mpi_total_max, mpi_comm_max, mpi_comp_max;
    MPI_Reduce(&tiempos_mpi.total, &mpi_total_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tiempos_mpi.comunicacion, &mpi_comm_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tiempos_mpi.computo, &mpi_comp_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); // Sincronizar antes de pasar a OMP

    // 2. Ejecutar OpenMP SOLAMENTE en el nodo maestro (Rank 0)
    double tiempo_omp = 0.0;
    int num_threads = 0;
    if (rank == 0) {
        num_threads = omp_get_max_threads();
        printf("[OMP] Ejecutando con %d hilos (Rank 0 unicamente)...\n", num_threads);
        tiempo_omp = ejecutar_omp();
    }

    // 3. Reporte Final (Solo Rank 0 imprime)
    if (rank == 0) {
        printf("\n========================================\n");
        printf("  RESULTADOS COMPARATIVOS\n");
        printf("========================================\n");
        
        printf("1. MPI (%d Procesos):\n", size);
        printf("   - Tiempo Total        : %12.6f s\n", mpi_total_max);
        printf("   - Tiempo Computo      : %12.6f s\n", mpi_comp_max);
        printf("   - Tiempo Comunicacion : %12.6f s (%.2f%%)\n", mpi_comm_max, (mpi_comm_max/mpi_total_max)*100.0);
        
        printf("\n2. OpenMP (%d Hilos):\n", num_threads);
        printf("   - Tiempo Total        : %12.6f s\n", tiempo_omp);

        printf("\n----------------------------------------\n");
        if (tiempo_omp < mpi_total_max) {
            printf("GANADOR: OpenMP (%.2fx mas rapido)\n", mpi_total_max / tiempo_omp);
            printf("Nota: Para grids pequeños (%dx%d), la latencia de red en MPI suele dominar.\n", imax, kmax);
        } else {
            printf("GANADOR: MPI (%.2fx mas rapido)\n", tiempo_omp / mpi_total_max);
        }
        printf("========================================\n");
    }

    MPI_Finalize();
    return 0;
}