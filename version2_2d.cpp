#include <iostream>
#include <mpi.h>
#include <cmath>

#define min(A, B) ((A) < (B) ? (A) : (B))
#define max(A, B) ((A) > (B) ? (A) : (B))

const int imax  = 80;
const int kmax  = 80;
const int itmax = 20000;

// Tags para comunicación MPI
#define TAG_UP     100
#define TAG_DOWN   101
#define TAG_LEFT   102
#define TAG_RIGHT  103

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double eps = 1.0e-08;
    double dx, dy, dx2, dy2, dx2i, dy2i, dt;

    // Parámetros físicos
    dx   = 1.0 / kmax;
    dy   = 1.0 / imax;
    dx2  = dx * dx;
    dy2  = dy * dy;
    dx2i = 1.0 / dx2;
    dy2i = 1.0 / dy2;
    dt   = min(dx2, dy2) / 4.0;

    // Descomposición 2D: calcular grid px × py ≈ √size
    int px = (int)std::sqrt((double)size);
    while (px > 1 && (size % px != 0)) px--;
    int py = size / px;

    int rank_x = rank % px;
    int rank_y = rank / px;

    // Distribución de dominio: i = 1..imax-1, k = 1..kmax-1
    int filas_totales = imax - 1;
    int cols_totales  = kmax - 1;

    // Distribución de filas (dimensión i)
    int filas_base  = filas_totales / py;
    int filas_extra = filas_totales % py;
    int filas_locales = (rank_y < filas_extra) ? filas_base + 1 : filas_base;
    int fila_inicio_global = (rank_y < filas_extra)
        ? rank_y * (filas_base + 1) + 1
        : filas_extra * (filas_base + 1) + (rank_y - filas_extra) * filas_base + 1;

    // Distribución de columnas (dimensión k)
    int cols_base  = cols_totales / px;
    int cols_extra = cols_totales % px;
    int cols_locales = (rank_x < cols_extra) ? cols_base + 1 : cols_base;
    int col_inicio_global = (rank_x < cols_extra)
        ? rank_x * (cols_base + 1) + 1
        : cols_extra * (cols_base + 1) + (rank_x - cols_extra) * cols_base + 1;

    // Vecinos en el grid 2D
    int vecino_arriba = (rank_y > 0)        ? rank - px : -1;
    int vecino_abajo  = (rank_y < py - 1)   ? rank + px : -1;
    int vecino_izq    = (rank_x > 0)        ? rank - 1  : -1;
    int vecino_der    = (rank_x < px - 1)   ? rank + 1  : -1;

    // Allocación con halos: phi_local[0..filas_locales+1][0..cols_locales+1]
    // Interiores: i_local=1..filas_locales, k_local=1..cols_locales
    double** phi_local  = new double*[filas_locales + 2];
    double** phin_local = new double*[filas_locales + 2];
    for (int i = 0; i < filas_locales + 2; i++) {
        phi_local[i]  = new double[cols_locales + 2];
        phin_local[i] = new double[cols_locales + 2];
    }

    // Tipos derivados MPI: fila contigua y columna con stride
    MPI_Datatype MPI_ROW_TYPE, MPI_COL_TYPE;
    MPI_Type_contiguous(cols_locales + 2, MPI_DOUBLE, &MPI_ROW_TYPE);
    MPI_Type_commit(&MPI_ROW_TYPE);
    MPI_Type_vector(filas_locales + 2, 1, cols_locales + 2, MPI_DOUBLE, &MPI_COL_TYPE);
    MPI_Type_commit(&MPI_COL_TYPE);

    // Inicialización: todo a 0
    for (int k_local = 0; k_local <= cols_locales + 1; k_local++) {
        for (int i_local = 0; i_local <= filas_locales + 1; i_local++) {
            phi_local[i_local][k_local]  = 0.0;
            phin_local[i_local][k_local] = 0.0;
        }
    }

    // Borde superior (k = kmax): phi[i][kmax] = 1.0
    int k_global_max = col_inicio_global + cols_locales;
    if (k_global_max == kmax) {
        for (int i_local = 0; i_local <= filas_locales + 1; i_local++) {
            phi_local[i_local][cols_locales + 1] = 1.0;
        }
    }

    // Esquinas: phi[0][0] = 0.0, phi[imax][0] = 0.0
    int i_global_min = fila_inicio_global - 1;
    int i_global_max = fila_inicio_global + filas_locales;
    int k_global_min = col_inicio_global - 1;
    
    if (i_global_min == 0 && k_global_min == 0) {
        phi_local[0][0] = 0.0;
    }
    if (i_global_max == imax && k_global_min == 0) {
        phi_local[filas_locales + 1][0] = 0.0;
    }

    // Bordes laterales: phi[0][k] = phi[0][k-1] + dx, phi[imax][k] = phi[imax][k-1] + dx
    if (i_global_min == 0) {
        for (int k_local = 1; k_local <= cols_locales; k_local++) {
            int k_global = col_inicio_global + k_local - 1;
            if (k_global > 0 && k_global < kmax) {
                phi_local[0][k_local] = phi_local[0][k_local - 1] + dx;
            }
        }
    }
    if (i_global_max == imax) {
        for (int k_local = 1; k_local <= cols_locales; k_local++) {
            int k_global = col_inicio_global + k_local - 1;
            if (k_global > 0 && k_global < kmax) {
                phi_local[filas_locales + 1][k_local] = phi_local[filas_locales + 1][k_local - 1] + dx;
            }
        }
    }

    // Inicializar halos que son fronteras físicas (cuando no hay vecino)
    if (vecino_izq < 0 && col_inicio_global == 1) {
        for (int i_local = 0; i_local <= filas_locales + 1; i_local++) {
            phi_local[i_local][0] = 0.0;
        }
    }
    if (vecino_der < 0 && col_inicio_global + cols_locales == kmax) {
        for (int i_local = 0; i_local <= filas_locales + 1; i_local++) {
            phi_local[i_local][cols_locales + 1] = 1.0;
        }
    }
    if (vecino_arriba < 0 && fila_inicio_global == 1) {
        phi_local[0][0] = 0.0;
        for (int k_local = 1; k_local <= cols_locales; k_local++) {
            int k_global = col_inicio_global + k_local - 1;
            if (k_global > 0 && k_global < kmax) {
                phi_local[0][k_local] = phi_local[0][k_local - 1] + dx;
            }
        }
        if (col_inicio_global + cols_locales == kmax) {
            phi_local[0][cols_locales + 1] = 1.0;
        }
    }
    if (vecino_abajo < 0 && fila_inicio_global + filas_locales == imax) {
        phi_local[filas_locales + 1][0] = 0.0;
        for (int k_local = 1; k_local <= cols_locales; k_local++) {
            int k_global = col_inicio_global + k_local - 1;
            if (k_global > 0 && k_global < kmax) {
                phi_local[filas_locales + 1][k_local] = phi_local[filas_locales + 1][k_local - 1] + dx;
            }
        }
        if (col_inicio_global + cols_locales == kmax) {
            phi_local[filas_locales + 1][cols_locales + 1] = 1.0;
        }
    }

    // Información inicial
    if (rank == 0) {
        printf("\n========================================\n");
        printf("  Transmision de calor 2D - Version 2D\n");
        printf("  Paralelizacion: MPI No-Bloqueante, 2D por bloques\n");
        printf("  Tipos derivados MPI: MPI_ROW_TYPE, MPI_COL_TYPE\n");
        printf("========================================\n");
        printf("\ndx = %12.4g, dy = %12.4g, dt = %12.4g, eps = %12.4g\n",
               dx, dy, dt, eps);
        printf("Numero de procesos: %d\n", size);
        printf("Grid de procesos: %d x %d (px x py)\n", px, py);
        printf("Filas interiores totales: %d (i=1..%d)\n", filas_totales, imax - 1);
        printf("Columnas interiores totales: %d (k=1..%d)\n", cols_totales, kmax - 1);
        printf("\nDistribucion del dominio:\n");
    }

    for (int r = 0; r < size; r++) {
        if (rank == r) {
            int fila_fin_global = fila_inicio_global + filas_locales - 1;
            int col_fin_global  = col_inicio_global + cols_locales - 1;
            printf("  Proceso %d (grid[%d,%d]): bloque i=[%d-%d], k=[%d-%d] (%dx%d)\n",
                   rank, rank_x, rank_y,
                   fila_inicio_global, fila_fin_global,
                   col_inicio_global, col_fin_global,
                   filas_locales, cols_locales);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    double t_inicio = MPI_Wtime();

    // Bucle temporal principal
    int it;
    for (it = 1; it <= itmax; it++) {
        double dphimax_local  = 0.0;
        double dphimax_global = 0.0;

        // Comunicación no-bloqueante de halos
        MPI_Request requests[8];
        int req_count = 0;

        // Enviar/recibir filas (arriba/abajo)
        if (vecino_arriba >= 0 && filas_locales > 0) {
            MPI_Isend(&phi_local[1][0], 1, MPI_ROW_TYPE, vecino_arriba, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&phi_local[0][0], 1, MPI_ROW_TYPE, vecino_arriba, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
        }
        if (vecino_abajo >= 0 && filas_locales > 0) {
            MPI_Isend(&phi_local[filas_locales][0], 1, MPI_ROW_TYPE, vecino_abajo, TAG_UP, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&phi_local[filas_locales + 1][0], 1, MPI_ROW_TYPE, vecino_abajo, TAG_DOWN, MPI_COMM_WORLD, &requests[req_count++]);
        }

        // Enviar/recibir columnas (izquierda/derecha)
        if (vecino_izq >= 0 && cols_locales > 0) {
            MPI_Isend(&phi_local[0][1], 1, MPI_COL_TYPE, vecino_izq, TAG_RIGHT, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&phi_local[0][0], 1, MPI_COL_TYPE, vecino_izq, TAG_LEFT, MPI_COMM_WORLD, &requests[req_count++]);
        }
        if (vecino_der >= 0 && cols_locales > 0) {
            MPI_Isend(&phi_local[0][cols_locales], 1, MPI_COL_TYPE, vecino_der, TAG_LEFT, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&phi_local[0][cols_locales + 1], 1, MPI_COL_TYPE, vecino_der, TAG_RIGHT, MPI_COMM_WORLD, &requests[req_count++]);
        }

        // Esperar halos antes de calcular
        if (req_count > 0) {
            MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
        }

        // Cálculo: solo puntos interiores (no fronteras físicas)
        for (int k_local = 1; k_local <= cols_locales; k_local++) {
            for (int i_local = 1; i_local <= filas_locales; i_local++) {
                int i_global = fila_inicio_global + (i_local - 1);
                int k_global = col_inicio_global + (k_local - 1);
                
                if (i_global >= 1 && i_global <= imax - 1 &&
                    k_global >= 1 && k_global <= kmax - 1) {
                    double dphi =
                        (phi_local[i_local + 1][k_local] + phi_local[i_local - 1][k_local] -
                         2.0 * phi_local[i_local][k_local]) * dy2i +
                        (phi_local[i_local][k_local + 1] + phi_local[i_local][k_local - 1] -
                         2.0 * phi_local[i_local][k_local]) * dx2i;

                    dphi = dphi * dt;
                    dphimax_local = max(dphimax_local, std::fabs(dphi));
                    phin_local[i_local][k_local] = phi_local[i_local][k_local] + dphi;
                }
            }
        }

        // Actualización: solo puntos interiores
        for (int k_local = 1; k_local <= cols_locales; k_local++) {
            for (int i_local = 1; i_local <= filas_locales; i_local++) {
                int i_global = fila_inicio_global + (i_local - 1);
                int k_global = col_inicio_global + (k_local - 1);
                
                if (i_global >= 1 && i_global <= imax - 1 &&
                    k_global >= 1 && k_global <= kmax - 1) {
                    phi_local[i_local][k_local] = phin_local[i_local][k_local];
                }
            }
        }

        // Reducción global para convergencia
        MPI_Allreduce(&dphimax_local, &dphimax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

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

    MPI_Type_free(&MPI_ROW_TYPE);
    MPI_Type_free(&MPI_COL_TYPE);

    MPI_Finalize();
    return 0;
}
