#include "common.hpp"
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <iostream>

// Tags para comunicación MPI
#define TAG_UP   100
#define TAG_DOWN 101

MetricasMPI ejecutar_mpi(int rank, int size) {
    double dx, dy, dx2, dy2, dx2i, dy2i, dt;
    int k, it;

    // Parámetros físicos
    dx   = 1.0 / kmax;
    dy   = 1.0 / imax;
    dx2  = dx * dx;
    dy2  = dy * dy;
    dx2i = 1.0 / dx2;
    dy2i = 1.0 / dy2;
    dt   = std::min(dx2, dy2) / 4.0;

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
    // CREACIÓN DE TIPOS DERIVADOS MPI
    // ============================================
    // Tipo derivado para una fila completa (kmax+1 elementos double)
    // Esto permite comunicar filas completas de manera más eficiente
    MPI_Datatype MPI_ROW_TYPE;
    MPI_Type_contiguous(kmax + 1, MPI_DOUBLE, &MPI_ROW_TYPE);
    MPI_Type_commit(&MPI_ROW_TYPE);

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
        printf("  Tipos derivados MPI: MPI_ROW_TYPE\n");
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

    // Sincronizar antes de medir
    MPI_Barrier(MPI_COMM_WORLD);

    // Variables de medición
    double t_inicio = MPI_Wtime();
    double t_isend_irecv = 0.0;
    double t_waitall = 0.0;
    double t_allreduce = 0.0;

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

        double t_aux = MPI_Wtime();

        // If this rank has at least one interior row, we may send/recv halos.
        if (filas_locales >= 1) {

            // Enviar primera fila interior (i_local = 1) al vecino superior
            if (rank > 0) {
                MPI_Isend(phi_local[1], 1, MPI_ROW_TYPE,
                          rank - 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
            }

            // Enviar última fila interior (i_local = filas_locales) al vecino inferior
            if (rank < size - 1) {
                MPI_Isend(phi_local[filas_locales], 1, MPI_ROW_TYPE,
                          rank + 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
            }

            // Recibir halo superior desde vecino superior
            if (rank > 0) {
                MPI_Irecv(phi_local[0], 1, MPI_ROW_TYPE,
                          rank - 1, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
            }

            // Recibir halo inferior desde vecino inferior
            if (rank < size - 1) {
                MPI_Irecv(phi_local[filas_locales + 1], 1, MPI_ROW_TYPE,
                          rank + 1, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
            }
        }

        t_isend_irecv += (MPI_Wtime() - t_aux);

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
                    dphimax_local = std::max(dphimax_local, std::fabs(dphi));
                    phin_local[i_local][k] = phi_local[i_local][k] + dphi;
                }
            }
        }

        // ============================================
        // ESPERAR A QUE LLEGUEN LOS HALOS
        // ============================================
        t_aux = MPI_Wtime();
        if (req_count > 0) {
            MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
        }
        t_waitall += (MPI_Wtime() - t_aux);

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
                dphimax_local = std::max(dphimax_local, std::fabs(dphi));
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
                dphimax_local = std::max(dphimax_local, std::fabs(dphi));
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
        t_aux = MPI_Wtime();
        MPI_Allreduce(&dphimax_local, &dphimax_global,
                      1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        t_allreduce += (MPI_Wtime() - t_aux);

        if (dphimax_global < eps) {
            break;
        }
    } // for it

    double t_fin = MPI_Wtime();

    // ============================================
    // LIMPIEZA DE MEMORIA Y TIPOS
    // ============================================
    for (int i = 0; i < filas_locales + 2; i++) {
        delete[] phi_local[i];
        delete[] phin_local[i];
    }
    delete[] phi_local;
    delete[] phin_local;

    MPI_Type_free(&MPI_ROW_TYPE);

    // ============================================
    // CÁLCULO DE MÉTRICAS
    // ============================================
    double t_comunicacion_total = t_isend_irecv + t_waitall + t_allreduce;

    MetricasMPI metricas;
    metricas.tiempo_total        = t_fin - t_inicio;
    metricas.tiempo_comunicacion = t_comunicacion_total;
    metricas.tiempo_computo      = metricas.tiempo_total - t_comunicacion_total;
    metricas.iteraciones         = it;

    long long flops = calcular_flops(imax, kmax, it);
    metricas.gflops = (double)flops / (metricas.tiempo_total * 1e9);
    metricas.porcentaje_comunicacion = (metricas.tiempo_comunicacion / metricas.tiempo_total) * 100.0;

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

        // Export CSV using the same helper as the other bench pieces
        export_mpi_csv("plots/mpi_results.csv", size, global);
        std::cout << "MPI done. Wrote plots/mpi_results.csv\n";
    }

    MPI_Finalize();
    return 0;
}
