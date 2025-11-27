#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <fstream>
#include <iomanip>

// Macros de utilidad
#define min(A, B) ((A) < (B) ? (A) : (B))
#define max(A, B) ((A) > (B) ? (A) : (B))

// Constantes del problema
const int imax = 80;
const int kmax = 80;
const int itmax = 20000;
const double eps = 1.0e-08;

// Tags MPI
#define TAG_UP   100
#define TAG_DOWN 101

// Estructura para métricas MPI
struct MetricasMPI {
    double tiempo_total;
    double tiempo_computo;
    double tiempo_comunicacion;
    int iteraciones;
    double gflops;
    double porcentaje_comunicacion;
};

// Estructura para métricas OpenMP
struct MetricasOMP {
    double tiempo_total;
    int iteraciones;
    double gflops;
};

// Estructura para métricas Secuencial
struct MetricasSecuencial {
    double tiempo_total;
    int iteraciones;
    double gflops;
};

// ======================================================
// FUNCIÓN: Calcular FLOPs
// ======================================================
long long calcular_flops(int filas, int columnas, int iteraciones) {
    // Por cada punto interno, por iteración:
    // Stencil en i: 2 sumas + 4 restas + 2 multiplicaciones = 8 ops
    // Stencil en k: 2 sumas + 4 restas + 2 multiplicaciones = 8 ops
    // dphi total: 1 suma + 1 multiplicación (dt) = 2 ops
    // Actualización: 1 suma = 1 op
    // Total: 19 FLOPs por punto por iteración
    
    long long puntos_internos = (long long)(filas - 1) * (columnas - 1);
    return puntos_internos * iteraciones * 19LL;
}

// ======================================================
// FUNCIÓN 1: SOLUCIONADOR MPI CON MEDICIÓN DETALLADA
// ======================================================
MetricasMPI ejecutar_mpi(int rank, int size) {
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
        phi_local[i][kmax] = 1.0;
    }
    
    if (rank == 0) phi_local[0][0] = 0.0;
    if (rank == size - 1) phi_local[filas_locales + 1][0] = 0.0;

    if (rank == 0) {
        for (k = 1; k < kmax; k++) 
            phi_local[0][k] = phi_local[0][k - 1] + dx;
    }
    if (rank == size - 1) {
        for (k = 1; k < kmax; k++) 
            phi_local[filas_locales + 1][k] = phi_local[filas_locales + 1][k - 1] + dx;
    }

    // Sincronizar antes de medir
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Variables de medición
    double t_inicio = MPI_Wtime();
    double t_comunicacion_total = 0.0;
    double t_isend_irecv = 0.0;
    double t_waitall = 0.0;
    double t_allreduce = 0.0;
    double t_aux;

    // Bucle temporal
    for (it = 1; it <= itmax; it++) {
        double dphimax_local = 0.0;
        double dphimax_global = 0.0;

        MPI_Request requests[4];
        int req_count = 0;

        // ==========================================
        // COMUNICACIÓN: Isend/Irecv
        // ==========================================
        t_aux = MPI_Wtime();
        
        if (rank > 0 && filas_locales > 0) {
            MPI_Isend(phi_local[1], kmax + 1, MPI_DOUBLE, 
                     rank - 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
        }
        if (rank < size - 1 && filas_locales > 0) {
            MPI_Isend(phi_local[filas_locales], kmax + 1, MPI_DOUBLE, 
                     rank + 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
        }
        if (rank > 0 && filas_locales > 0) {
            MPI_Irecv(phi_local[0], kmax + 1, MPI_DOUBLE, 
                     rank - 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
        }
        if (rank < size - 1 && filas_locales > 0) {
            MPI_Irecv(phi_local[filas_locales + 1], kmax + 1, MPI_DOUBLE, 
                     rank + 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
        }
        
        t_isend_irecv += (MPI_Wtime() - t_aux);

        // ==========================================
        // CÓMPUTO: Puntos interiores
        // ==========================================
        int i_inicio = 2;
        int i_fin = filas_locales - 1;

        if (i_fin >= i_inicio) {
            for (k = 1; k < kmax; k++) {
                for (int i_local = i_inicio; i_local <= i_fin; i_local++) {
                    double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k] 
                                  - 2.0 * phi_local[i_local][k]) * dy2i +
                                  (phi_local[i_local][k + 1] + phi_local[i_local][k - 1] 
                                  - 2.0 * phi_local[i_local][k]) * dx2i;
                    dphi *= dt;
                    dphimax_local = max(dphimax_local, std::fabs(dphi));
                    phin_local[i_local][k] = phi_local[i_local][k] + dphi;
                }
            }
        }

        // ==========================================
        // COMUNICACIÓN: Waitall
        // ==========================================
        t_aux = MPI_Wtime();
        if (req_count > 0) {
            MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
        }
        t_waitall += (MPI_Wtime() - t_aux);

        // ==========================================
        // CÓMPUTO: Puntos frontera
        // ==========================================
        if (filas_locales >= 1) {
            int i_local = 1;
            for (k = 1; k < kmax; k++) {
                double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k] 
                              - 2.0 * phi_local[i_local][k]) * dy2i +
                              (phi_local[i_local][k + 1] + phi_local[i_local][k - 1] 
                              - 2.0 * phi_local[i_local][k]) * dx2i;
                dphi *= dt;
                dphimax_local = max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }

        if (filas_locales >= 2) {
            int i_local = filas_locales;
            for (k = 1; k < kmax; k++) {
                double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k] 
                              - 2.0 * phi_local[i_local][k]) * dy2i +
                              (phi_local[i_local][k + 1] + phi_local[i_local][k - 1] 
                              - 2.0 * phi_local[i_local][k]) * dx2i;
                dphi *= dt;
                dphimax_local = max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }

        // Actualización
        for (int i_local = 1; i_local <= filas_locales; i_local++) {
            for (k = 1; k < kmax; k++) {
                phi_local[i_local][k] = phin_local[i_local][k];
            }
        }

        // ==========================================
        // COMUNICACIÓN: Allreduce
        // ==========================================
        t_aux = MPI_Wtime();
        MPI_Allreduce(&dphimax_local, &dphimax_global, 
                      1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        t_allreduce += (MPI_Wtime() - t_aux);

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

    // Calcular métricas
    t_comunicacion_total = t_isend_irecv + t_waitall + t_allreduce;
    
    MetricasMPI metricas;
    metricas.tiempo_total = t_final - t_inicio;
    metricas.tiempo_comunicacion = t_comunicacion_total;
    metricas.tiempo_computo = metricas.tiempo_total - t_comunicacion_total;
    metricas.iteraciones = it;
    
    long long flops = calcular_flops(imax, kmax, it);
    metricas.gflops = (double)flops / (metricas.tiempo_total * 1e9);
    metricas.porcentaje_comunicacion = (t_comunicacion_total / metricas.tiempo_total) * 100.0;

    return metricas;
}

// ======================================================
// FUNCIÓN 0: SOLUCIONADOR SECUENCIAL (BASELINE)
// ======================================================
MetricasSecuencial ejecutar_secuencial() {
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
    for (k = 0; k < kmax; k++) 
        for (i = 1; i < imax; i++) 
            phi[i][k] = 0.0;
    
    for (i = 0; i <= imax; i++) 
        phi[i][kmax] = 1.0;
    
    phi[0][0] = 0.0; 
    phi[imax][0] = 0.0;
    
    for (k = 1; k < kmax; k++) {
        phi[0][k] = phi[0][k - 1] + dx;
        phi[imax][k] = phi[imax][k - 1] + dx;
    }

    double t_inicio = omp_get_wtime();

    // Bucle temporal sin paralelización
    for (it = 1; it <= itmax; it++) {
        dphimax = 0.;
        
        for (k = 1; k < kmax; k++) {
            for (i = 1; i < imax; i++) {
                dphi = (phi[i + 1][k] + phi[i - 1][k] - 2. * phi[i][k]) * dy2i + 
                       (phi[i][k + 1] + phi[i][k - 1] - 2. * phi[i][k]) * dx2i;
                dphi = dphi * dt;
                dphimax = max(dphimax, dphi);
                phin[i][k] = phi[i][k] + dphi;
            }
        }

        for (k = 1; k < kmax; k++) {
            for (i = 1; i < imax; i++) {
                phi[i][k] = phin[i][k];
            }
        }
        
        if (dphimax < eps) break;
    }

    double t_fin = omp_get_wtime();

    MetricasSecuencial metricas;
    metricas.tiempo_total = t_fin - t_inicio;
    metricas.iteraciones = it;
    
    long long flops = calcular_flops(imax, kmax, it);
    metricas.gflops = (double)flops / (metricas.tiempo_total * 1e9);

    return metricas;
}

// ======================================================
// FUNCIÓN 2: SOLUCIONADOR OPENMP CON MEDICIÓN
// ======================================================
MetricasOMP ejecutar_omp() {
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
    for (k = 0; k < kmax; k++) 
        for (i = 1; i < imax; i++) 
            phi[i][k] = 0.0;
    
    for (i = 0; i <= imax; i++) 
        phi[i][kmax] = 1.0;
    
    phi[0][0] = 0.0; 
    phi[imax][0] = 0.0;
    
    for (k = 1; k < kmax; k++) {
        phi[0][k] = phi[0][k - 1] + dx;
        phi[imax][k] = phi[imax][k - 1] + dx;
    }

    double t_inicio = omp_get_wtime();

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

    MetricasOMP metricas;
    metricas.tiempo_total = t_fin - t_inicio;
    metricas.iteraciones = it;
    
    long long flops = calcular_flops(imax, kmax, it);
    metricas.gflops = (double)flops / (metricas.tiempo_total * 1e9);

    return metricas;
}

// ======================================================
// FUNCIÓN: Exportar CSV para análisis
// ======================================================
void exportar_csv(int num_procesos, int num_threads, 
                  const MetricasSecuencial& seq,
                  const MetricasMPI& mpi, const MetricasOMP& omp,
                  double speedup_mpi, double eficiencia_mpi,
                  double speedup_omp, double eficiencia_omp) {
    std::ofstream archivo("resultados_benchmark.csv", std::ios::app);
    
    // Si el archivo está vacío, escribir encabezado
    archivo.seekp(0, std::ios::end);
    if (archivo.tellp() == 0) {
        archivo << "Procesos,Threads,T_Secuencial,T_MPI,T_OMP,"
                << "Speedup_MPI,Eficiencia_MPI,Speedup_OMP,Eficiencia_OMP,"
                << "GFlops_Seq,GFlops_MPI,GFlops_OMP,Comunicacion_%,Iteraciones\n";
    }
    
    archivo << std::fixed << std::setprecision(6);
    archivo << num_procesos << ","
            << num_threads << ","
            << seq.tiempo_total << ","
            << mpi.tiempo_total << ","
            << omp.tiempo_total << ","
            << speedup_mpi << ","
            << eficiencia_mpi << ","
            << speedup_omp << ","
            << eficiencia_omp << ","
            << seq.gflops << ","
            << mpi.gflops << ","
            << omp.gflops << ","
            << mpi.porcentaje_comunicacion << ","
            << mpi.iteraciones << "\n";
    
    archivo.close();
}

// ======================================================
// MAIN: BENCHMARK Y ANÁLISIS
// ======================================================
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("\n========================================\n");
        printf("  BENCHMARK DETALLADO: Seq vs MPI vs OMP\n");
        printf("  Grid: %dx%d | Max Iteraciones: %d\n", imax, kmax, itmax);
        printf("========================================\n\n");
    }

    // ==========================================
    // FASE 0: EJECUTAR SECUENCIAL (Baseline)
    // ==========================================
    MetricasSecuencial metricas_seq;
    
    if (rank == 0) {
        printf("[SEQ] Ejecutando versión secuencial (baseline)...\n");
        metricas_seq = ejecutar_secuencial();
        printf("      Tiempo: %.6f s | Iteraciones: %d\n\n", 
               metricas_seq.tiempo_total, metricas_seq.iteraciones);
    }

    // Broadcast del tiempo secuencial a todos los procesos
    MPI_Bcast(&metricas_seq.tiempo_total, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&metricas_seq.iteraciones, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&metricas_seq.gflops, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // ==========================================
    // FASE 1: EJECUTAR MPI
    // ==========================================
    if (rank == 0) printf("[MPI] Ejecutando con %d procesos...\n", size);
    MetricasMPI metricas_mpi_local = ejecutar_mpi(rank, size);

    // Recolectar tiempos máximos de todos los procesos
    MetricasMPI metricas_mpi;
    MPI_Reduce(&metricas_mpi_local.tiempo_total, &metricas_mpi.tiempo_total, 
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&metricas_mpi_local.tiempo_comunicacion, &metricas_mpi.tiempo_comunicacion, 
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&metricas_mpi_local.tiempo_computo, &metricas_mpi.tiempo_computo, 
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        metricas_mpi.iteraciones = metricas_mpi_local.iteraciones;
        long long flops = calcular_flops(imax, kmax, metricas_mpi.iteraciones);
        metricas_mpi.gflops = (double)flops / (metricas_mpi.tiempo_total * 1e9);
        metricas_mpi.porcentaje_comunicacion = 
            (metricas_mpi.tiempo_comunicacion / metricas_mpi.tiempo_total) * 100.0;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ==========================================
    // FASE 2: EJECUTAR OPENMP (solo Rank 0)
    // ==========================================
    MetricasOMP metricas_omp;
    int num_threads = 0;
    
    if (rank == 0) {
        num_threads = omp_get_max_threads();
        printf("[OMP] Ejecutando con %d threads...\n", num_threads);
        metricas_omp = ejecutar_omp();
        printf("      Tiempo: %.6f s | Iteraciones: %d\n\n", 
               metricas_omp.tiempo_total, metricas_omp.iteraciones);
    }

    // ==========================================
    // FASE 3: ANÁLISIS Y REPORTE
    // ==========================================
    if (rank == 0) {
        printf("========================================\n");
        printf("  RESULTADOS DETALLADOS\n");
        printf("========================================\n\n");
        
        // --- SECUENCIAL ---
        printf("0. SECUENCIAL (Baseline):\n");
        printf("   Tiempo Total        : %12.6f s\n", metricas_seq.tiempo_total);
        printf("   Iteraciones         : %12d\n", metricas_seq.iteraciones);
        printf("   Rendimiento         : %12.6f GFlops\n\n", metricas_seq.gflops);
        
        // --- MPI ---
        printf("1. MPI (%d Procesos):\n", size);
        printf("   Tiempo Total        : %12.6f s\n", metricas_mpi.tiempo_total);
        printf("   Tiempo Cómputo      : %12.6f s (%5.2f%%)\n", 
               metricas_mpi.tiempo_computo,
               (metricas_mpi.tiempo_computo / metricas_mpi.tiempo_total) * 100.0);
        printf("   Tiempo Comunicación : %12.6f s (%5.2f%%)\n", 
               metricas_mpi.tiempo_comunicacion,
               metricas_mpi.porcentaje_comunicacion);
        printf("   Iteraciones         : %12d\n", metricas_mpi.iteraciones);
        printf("   Rendimiento         : %12.6f GFlops\n\n", metricas_mpi.gflops);
        
        // --- OpenMP ---
        printf("2. OpenMP (%d Threads):\n", num_threads);
        printf("   Tiempo Total        : %12.6f s\n", metricas_omp.tiempo_total);
        printf("   Iteraciones         : %12d\n", metricas_omp.iteraciones);
        printf("   Rendimiento         : %12.6f GFlops\n\n", metricas_omp.gflops);
        
        // --- Métricas Comparativas (CORRECTO: vs Secuencial) ---
        double speedup_mpi = metricas_seq.tiempo_total / metricas_mpi.tiempo_total;
        double eficiencia_mpi = (speedup_mpi / size) * 100.0;
        
        double speedup_omp = metricas_seq.tiempo_total / metricas_omp.tiempo_total;
        double eficiencia_omp = (speedup_omp / num_threads) * 100.0;
        
        printf("========================================\n");
        printf("  MÉTRICAS DE RENDIMIENTO\n");
        printf("========================================\n");
        printf("MPI:\n");
        printf("  Speedup    (Tseq/Tmpi) : %.4fx\n", speedup_mpi);
        printf("  Eficiencia (Sp/p×100)  : %.2f%%\n", eficiencia_mpi);
        printf("\nOpenMP:\n");
        printf("  Speedup    (Tseq/Tomp) : %.4fx\n", speedup_omp);
        printf("  Eficiencia (Sp/p×100)  : %.2f%%\n", eficiencia_omp);
        
        printf("\n----------------------------------------\n");
        printf("Comparación MPI vs OpenMP:\n");
        if (metricas_mpi.tiempo_total < metricas_omp.tiempo_total) {
            double ventaja = metricas_omp.tiempo_total / metricas_mpi.tiempo_total;
            printf("✓ MPI es %.2fx más rápido que OpenMP\n", ventaja);
        } else {
            double ventaja = metricas_mpi.tiempo_total / metricas_omp.tiempo_total;
            printf("✓ OpenMP es %.2fx más rápido que MPI\n", ventaja);
        }
        
        if (metricas_mpi.porcentaje_comunicacion > 30.0) {
            printf("\n⚠ ADVERTENCIA: Overhead de comunicación alto (%.1f%%)\n", 
                   metricas_mpi.porcentaje_comunicacion);
            printf("  Recomendaciones:\n");
            printf("  - Aumentar tamaño de malla (imax/kmax)\n");
            printf("  - Reducir número de procesos MPI\n");
        }
        printf("========================================\n\n");
        
        // Exportar CSV
        exportar_csv(size, num_threads, metricas_seq, metricas_mpi, metricas_omp, 
                     speedup_mpi, eficiencia_mpi, speedup_omp, eficiencia_omp);
        printf("✓ Resultados exportados a: resultados_benchmark.csv\n\n");
    }

    MPI_Finalize();
    return 0;
}
