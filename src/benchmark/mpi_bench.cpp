// mpi_bench.cpp
#include "common.hpp"
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <iostream>

// Tags MPI
#define TAG_UP   100
#define TAG_DOWN 101

// ======================================================
// FUNCIÓN 1: SOLUCIONADOR MPI (idéntico al version3)
// ======================================================
MetricasMPI ejecutar_mpi(int rank, int size) {
    double dx, dy, dx2, dy2, dx2i, dy2i, dt;
    int k, it;

    // Parámetros físicos (igual que version3)
    dx = 1.0 / kmax;
    dy = 1.0 / imax;
    dx2 = dx * dx;
    dy2 = dy * dy;
    dx2i = 1.0 / dx2;
    dy2i = 1.0 / dy2;
    dt = _min(dx2, dy2) / 4.0;

    // Descomposición de dominio
    int filas_totales = imax - 1;
    int filas_base = filas_totales / size;
    int filas_extra = filas_totales % size;
    int filas_locales;

    if (rank < filas_extra) {
        filas_locales = filas_base + 1;
    } else {
        filas_locales = filas_base;
    }

    // Allocación (filas_locales + 2) x (kmax + 1)
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
        // COMUNICACIÓN: Isend / Irecv
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
                for (int i_local = i_inicio; i_local <= i_fin; ++i_local) {
                    double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k]
                                  - 2.0 * phi_local[i_local][k]) * dy2i +
                                  (phi_local[i_local][k + 1] + phi_local[i_local][k - 1]
                                  - 2.0 * phi_local[i_local][k]) * dx2i;
                    dphi *= dt;
                    dphimax_local = _max(dphimax_local, std::fabs(dphi));
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
            for (k = 1; k < kmax; ++k) {
                double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k]
                              - 2.0 * phi_local[i_local][k]) * dy2i +
                              (phi_local[i_local][k + 1] + phi_local[i_local][k - 1]
                              - 2.0 * phi_local[i_local][k]) * dx2i;
                dphi *= dt;
                dphimax_local = _max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }

        if (filas_locales >= 2) {
            int i_local = filas_locales;
            for (k = 1; k < kmax; ++k) {
                double dphi = (phi_local[i_local + 1][k] + phi_local[i_local - 1][k]
                              - 2.0 * phi_local[i_local][k]) * dy2i +
                              (phi_local[i_local][k + 1] + phi_local[i_local][k - 1]
                              - 2.0 * phi_local[i_local][k]) * dx2i;
                dphi *= dt;
                dphimax_local = _max(dphimax_local, std::fabs(dphi));
                phin_local[i_local][k] = phi_local[i_local][k] + dphi;
            }
        }

        // Actualización
        for (int i_local = 1; i_local <= filas_locales; ++i_local) {
            for (k = 1; k < kmax; ++k) {
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

    // Calcular métricas (igual que version3)
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

// MAIN
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) std::cout << "MPI benchmark starting with " << size << " ranks...\n";

    MetricasMPI local = ejecutar_mpi(rank, size);

    // Recolectar máximos a rank 0 (idéntico a version3)
    MetricasMPI global;
    MPI_Reduce(&local.tiempo_total, &global.tiempo_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local.tiempo_computo, &global.tiempo_computo, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local.tiempo_comunicacion, &global.tiempo_comunicacion, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        global.iteraciones = local.iteraciones;
        long long flops = calcular_flops(imax, kmax, global.iteraciones);
        global.gflops = (double)flops / (global.tiempo_total * 1e9);
        global.porcentaje_comunicacion = (global.tiempo_comunicacion / global.tiempo_total) * 100.0;

        // same CSV output name as your split workflow expects
        export_mpi_csv("plots/mpi_results.csv", size, global);
        std::cout << "MPI done. Wrote plots/mpi_results.csv\n";
    }

    MPI_Finalize();
    return 0;
}
