#include <iostream>
#include <mpi.h>
#include <cmath>

#define min(A, B) ((A) < (B) ? (A) : (B))
#define max(A, B) ((A) > (B) ? (A) : (B))

const int imax = 80;
const int kmax = 80;
const int itmax = 20000;

// Tags para comunicación MPI
#define TAG_UP   100
#define TAG_DOWN 101

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double eps = 1.0e-08;
    double dx, dy, dx2, dy2, dx2i, dy2i, dt;
    int k, it;

    // Parámetros físicos
    dx   = 1.0 / kmax;
    dy   = 1.0 / imax;
    dx2  = dx * dx;
    dy2  = dy * dy;
    dx2i = 1.0 / dx2;
    dy2i = 1.0 / dy2;
    dt   = min(dx2, dy2) / 4.0;

    // ============================================
    // DESCOMPOSICIÓN DEL DOMINIO EN 1D (POR FILAS)
    // ============================================
    int filas_totales = imax - 1;           // filas interiores globales: i = 1..imax-1
    int filas_base    = filas_totales / size;
    int filas_extra   = filas_totales % size;

    int filas_locales;
    int fila_inicio_global; // i_global de la primera fila interior de este proceso

    if (rank < filas_extra) {
        filas_locales      = filas_base + 1;
        fila_inicio_global = rank * filas_locales + 1;
    } else {
        filas_locales      = filas_base;
        fila_inicio_global = rank * filas_base + filas_extra + 1;
    }

    // ============================================
    // ALLOCACIÓN DE MEMORIA CON HALOS
    // ============================================
    // phi_local[0]              = halo superior (i_global = fila_inicio_global-1 o borde i=0)
    // phi_local[1..filas_local] = filas interiores propias
    // phi_local[filas_local+1]  = halo inferior (i_global = fila_inicio_global+filas_locales o borde i=imax)
    double** phi_local  = new double*[filas_locales + 2];
    double** phin_local = new double*[filas_locales + 2];

    for (int i = 0; i < filas_locales + 2; i++) {
        phi_local[i]  = new double[kmax + 1];
        phin_local[i] = new double[kmax + 1];
    }

    // ============================================
    // INICIALIZACIÓN DE CONDICIONES DE FRONTERA
    // Igual que el código secuencial, pero distribuido
    // ============================================

    // Primero, todo a 0
    for (int i = 0; i < filas_locales + 2; i++) {
        for (k = 0; k <= kmax; k++) {
            phi_local[i][k]  = 0.0;
            phin_local[i][k] = 0.0;
        }
    }

    // Borde superior (k = kmax): phi[i][kmax] = 1.0 para todo i (incluye halos)
    for (int i = 0; i < filas_locales + 2; i++) {
        phi_local[i][kmax] = 1.0;
    }

    // Esquinas i=0, i=imax en k=0
    if (rank == 0) {
        // i_global = 0 está en phi_local[0]
        phi_local[0][0] = 0.0;
    }
    if (rank == size - 1) {
        // i_global = imax está en phi_local[filas_locales+1]
        phi_local[filas_locales + 1][0] = 0.0;
    }

    // Bordes laterales (i=0 y i=imax) crecen con dx en k (como en secuencial)
    if (rank == 0) {
        // Borde i=0: phi[0][k] = phi[0][k-1] + dx
        for (k = 1; k < kmax; k++) {
            phi_local[0][k] = phi_local[0][k - 1] + dx;
        }
    }

    if (rank == size - 1) {
        // Borde i=imax: phi[imax][k] = phi[imax][k-1] + dx
        for (k = 1; k < kmax; k++) {
            phi_local[filas_locales + 1][k] = phi_local[filas_locales + 1][k - 1] + dx;
        }
    }

    // ============================================
    // INFORMACIÓN INICIAL
    // ============================================
    if (rank == 0) {
        printf("\n========================================\n");
        printf("  Transmision de calor 2D - Version 2\n");
        printf("  Paralelizacion: MPI No-Bloqueante, 1D por filas\n");
        printf("========================================\n");
        printf("\ndx = %12.4g, dy = %12.4g, dt = %12.4g, eps = %12.4g\n",
               dx, dy, dt, eps);
        printf("Numero de procesos: %d\n", size);
        printf("Filas interiores totales: %d (1..%d)\n", filas_totales, imax - 1);
        printf("\nDistribucion del dominio:\n");
    }

    // Cada proceso reporta su rango de filas globales
    for (int r = 0; r < size; r++) {
        if (rank == r) {
            int fila_fin_global = fila_inicio_global + filas_locales - 1;
            printf("  Proceso %d: filas globales %d-%d (%d filas)\n",
                   rank, fila_inicio_global, fila_fin_global, filas_locales);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    double t_inicio = MPI_Wtime();

    // ============================================
    // BUCLE TEMPORAL PRINCIPAL
    // ============================================
    for (it = 1; it <= itmax; it++) {
        double dphimax_local  = 0.0;
        double dphimax_global = 0.0;

        // ============================================
        // COMUNICACIÓN NO-BLOQUEANTE DE HALOS
        // ============================================
        MPI_Request requests[4];
        int req_count = 0;

        // Enviar primera fila interior (i_local = 1) al vecino superior
        if (rank > 0 && filas_locales > 0) {
            MPI_Isend(phi_local[1], kmax + 1, MPI_DOUBLE,
                      rank - 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
        }

        // Enviar última fila interior (i_local = filas_locales) al vecino inferior
        if (rank < size - 1 && filas_locales > 0) {
            MPI_Isend(phi_local[filas_locales], kmax + 1, MPI_DOUBLE,
                      rank + 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
        }

        // Recibir halo superior desde vecino superior
        if (rank > 0 && filas_locales > 0) {
            MPI_Irecv(phi_local[0], kmax + 1, MPI_DOUBLE,
                      rank - 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
        }

        // Recibir halo inferior desde vecino inferior
        if (rank < size - 1 && filas_locales > 0) {
            MPI_Irecv(phi_local[filas_locales + 1], kmax + 1, MPI_DOUBLE,
                      rank + 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
        }

        // ============================================
        // CÁLCULO DE PUNTOS INTERIORES (sin usar halos)
        // i_local = 2..filas_locales-1
        // ============================================
        int i_inicio = 2;
        int i_fin    = filas_locales - 1;

        if (i_fin >= i_inicio) {
            for (k = 1; k < kmax; k++) {
                for (int i_local = i_inicio; i_local <= i_fin; i_local++) {
                    double dphi =
                        (phi_local[i_local + 1][k] + phi_local[i_local - 1][k]
                         - 2.0 * phi_local[i_local][k]) * dy2i +
                        (phi_local[i_local][k + 1] + phi_local[i_local][k - 1]
                         - 2.0 * phi_local[i_local][k]) * dx2i;

                    dphi *= dt;
                    dphimax_local = max(dphimax_local, std::fabs(dphi));
                    phin_local[i_local][k] = phi_local[i_local][k] + dphi;
                }
            }
        }

        // ============================================
        // ESPERAR A QUE LLEGUEN LOS HALOS
        // ============================================
        if (req_count > 0) {
            MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
        }

        // ============================================
        // CÁLCULO DE PUNTOS QUE USAN HALOS (primer y último fila local)
        // ============================================
        // Primera fila interior local (i_local = 1), si existe
        if (filas_locales >= 1) {
            int i_local = 1;
            for (k = 1; k < kmax; k++) {
                double dphi =
                    (phi_local[i_local + 1][k] + phi_local[i_local - 1][k]
                     - 2.0 * phi_local[i_local][k]) * dy2i +
                    (phi_local[i_local][k + 1] + phi_local[i_local][k - 1]
                     - 2.0 * phi_local[i_local][k]) * dx2i;

                dphi *= dt;
                dphimax_local = max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }

        // Última fila interior local (i_local = filas_locales), si es distinta
        if (filas_locales >= 2) {
            int i_local = filas_locales;
            for (k = 1; k < kmax; k++) {
                double dphi =
                    (phi_local[i_local + 1][k] + phi_local[i_local - 1][k]
                     - 2.0 * phi_local[i_local][k]) * dy2i +
                    (phi_local[i_local][k + 1] + phi_local[i_local][k - 1]
                     - 2.0 * phi_local[i_local][k]) * dx2i;

                dphi *= dt;
                dphimax_local = max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }

        // ============================================
        // ACTUALIZACIÓN DE VALORES (copiar phin -> phi)
        // ============================================
        for (int i_local = 1; i_local <= filas_locales; i_local++) {
            for (k = 1; k < kmax; k++) {
                phi_local[i_local][k] = phin_local[i_local][k];
            }
        }

        // ============================================
        // REDUCCIÓN GLOBAL PARA CONVERGENCIA
        // ============================================
        MPI_Allreduce(&dphimax_local, &dphimax_global,
                      1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (dphimax_global < eps) {
            break;
        }
    }

    double t_fin = MPI_Wtime();

    if (rank == 0) {
        printf("\n%d iteraciones completadas\n", it);
        printf("\nTiempo de ejecucion = %12.6f sec\n", t_fin - t_inicio);
        printf("========================================\n");
    }

    // Liberar memoria
    for (int i = 0; i < filas_locales + 2; i++) {
        delete[] phi_local[i];
        delete[] phin_local[i];
    }
    delete[] phi_local;
    delete[] phin_local;

    MPI_Finalize();
    return 0;
}
