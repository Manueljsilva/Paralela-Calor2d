#include "common.hpp"
#include <omp.h>
#include <iostream>
#include <cmath>
#include <algorithm>

MetricasOMP ejecutar_omp() {
    const int NI = imax + 1;
    const int NK = kmax + 1;

    // Allocate arrays on heap to avoid stack overflow
    double *phi  = new double[(size_t)NI * NK];
    double *phin = new double[(size_t)NI * NK];
    auto idx = [&](int i, int k) -> size_t { return (size_t)i * NK + k; };

    // Grid spacing
    double dx = 1.0 / kmax;
    double dy = 1.0 / imax;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double dx2i = 1.0 / dx2;
    double dy2i = 1.0 / dy2;
    double dt = std::min(dx2, dy2) / 4.0;

    double eps = 1e-8;

    // ---- Initialization ----
    for (int i = 0; i < NI; ++i)
        for (int k = 0; k < NK; ++k)
            phi[idx(i,k)] = 0.0;

    for (int i = 0; i <= imax; ++i)
        phi[idx(i,kmax)] = 1.0;

    phi[idx(0,0)] = 0.0;
    phi[idx(imax,0)] = 0.0;

    for (int k = 1; k < kmax; ++k) {
        phi[idx(0,k)] = phi[idx(0,k-1)] + dx;
        phi[idx(imax,k)] = phi[idx(imax,k-1)] + dx;
    }

    // ---- Print header info like original ----
    std::cout << "  Transmision de calor 2D - Version 1\n";
    std::cout << "  Paralelizacion: OpenMP\n";
    std::cout << "\ndx = " << dx << ", dy = " << dy
              << ", dt = " << dt << ", eps = " << eps << "\n";

    int num_threads;
    #pragma omp parallel
    {
        #pragma omp single
        num_threads = omp_get_num_threads();
    }
    std::cout << "Numero de threads: " << num_threads << "\n";

    // ---- Iteration ----
    int it;
    double t_inicio = omp_get_wtime();
    for (it = 1; it <= itmax; ++it) {
        double dphimax = 0.0;

        #pragma omp parallel for reduction(max:dphimax) schedule(static)
        for (int k = 1; k < kmax; ++k) {
            for (int i = 1; i < imax; ++i) {
                double dphi = (phi[idx(i+1,k)] + phi[idx(i-1,k)] - 2.0 * phi[idx(i,k)]) * dy2i
                            + (phi[idx(i,k+1)] + phi[idx(i,k-1)] - 2.0 * phi[idx(i,k)]) * dx2i;
                dphi *= dt;
                dphimax = std::max(dphimax, std::fabs(dphi));
                phin[idx(i,k)] = phi[idx(i,k)] + dphi;
            }
        }

        // Copy phin -> phi
        #pragma omp parallel for schedule(static)
        for (int k = 1; k < kmax; ++k)
            for (int i = 1; i < imax; ++i)
                phi[idx(i,k)] = phin[idx(i,k)];

        if (dphimax < eps)
            break;
    }
    double t_fin = omp_get_wtime();

    // ---- Print results like original ----
    std::cout << "\n" << it << " iteraciones completadas\n";
    std::cout << "Tiempo de ejecucion = " << t_fin - t_inicio << " sec\n";

    // ---- Metrics struct for CSV ----
    MetricasOMP metricas;
    metricas.tiempo_total = t_fin - t_inicio;
    metricas.iteraciones = it;
    long long flops = calcular_flops(imax, kmax, it);
    metricas.gflops = (double)flops / (metricas.tiempo_total * 1e9);

    delete [] phi;
    delete [] phin;

    return metricas;
}

int main() {
    int threads = omp_get_max_threads();
    std::cout << "OMP benchmark starting with " << threads << " threads...\n";

    MetricasOMP m = ejecutar_omp();
    std::cout << "OMP done: time=" << m.tiempo_total
              << "s, it=" << m.iteraciones
              << ", gflops=" << m.gflops << "\n";

    export_omp_csv("plots/omp_results.csv", threads, m);
    std::cout << "Wrote plots/omp_results.csv\n";
    return 0;
}
