// omp_bench.cpp
#include "common.hpp"
#include <omp.h>
#include <iostream>
#include <cmath>
#include <algorithm>

// OpenMP solver (keeps behaviour of version3's OpenMP phase, but uses heap arrays)
MetricasOMP ejecutar_omp() {
    const int NI = imax + 1;
    const int NK = kmax + 1;
    // allocate on heap to avoid stack overflow for large imax/kmax
    double *phi = new double[(size_t)NI * NK];
    double *phin = new double[(size_t)NI * NK];
    auto idx = [&](int i, int k) -> size_t { return (size_t)i * NK + k; };

    double dx = 1.0 / kmax;
    double dy = 1.0 / imax;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double dx2i = 1.0 / dx2;
    double dy2i = 1.0 / dy2;
    // use std::min to compute dt (same formula as version3)
    double dt = std::min(dx2, dy2) / 4.0;

    // Initialization (same semantics as version3)
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

    double t_inicio = omp_get_wtime();
    int it;
    for (it = 1; it <= itmax; ++it) {
        double dphimax = 0.0;

        // compute interior points (same loop order, same reduction semantics)
        #pragma omp parallel for reduction(max:dphimax) schedule(static)
        for (int k = 1; k < kmax; ++k) {
            for (int i = 1; i < imax; ++i) {
                double dphi = (phi[idx(i+1,k)] + phi[idx(i-1,k)] - 2.0 * phi[idx(i,k)]) * dy2i
                            + (phi[idx(i,k+1)] + phi[idx(i,k-1)] - 2.0 * phi[idx(i,k)]) * dx2i;
                dphi *= dt;
                // mirror version3 behaviour: update dphimax with the new dphi
                if (dphi > dphimax) dphimax = dphi;
                phin[idx(i,k)] = phi[idx(i,k)] + dphi;
            }
        }

        // copy phin -> phi (same ordering, parallel)
        #pragma omp parallel for schedule(static)
        for (int k = 1; k < kmax; ++k) {
            for (int i = 1; i < imax; ++i) {
                phi[idx(i,k)] = phin[idx(i,k)];
            }
        }

        // convergence test identical to version3
        if (dphimax < eps) break;
    }
    double t_fin = omp_get_wtime();

    MetricasOMP metricas;
    metricas.tiempo_total = t_fin - t_inicio;
    metricas.iteraciones = it;
    long long flops = calcular_flops(imax, kmax, it);
    metricas.gflops = (double)flops / (metricas.tiempo_total * 1e9);

    delete [] phi;
    delete [] phin;
    return metricas;
}

int main(int argc, char** argv) {
    int threads = omp_get_max_threads();
    std::cout << "OMP benchmark starting with " << threads << " threads...\n";

    MetricasOMP m = ejecutar_omp();
    std::cout << "OMP done: time=" << m.tiempo_total << "s, it=" << m.iteraciones << ", gflops=" << m.gflops << "\n";

    // same CSV name & export helper as your other split files
    export_omp_csv("plots/omp_results.csv", threads, m);
    std::cout << "Wrote plots/omp_results.csv\n";
    return 0;
}
