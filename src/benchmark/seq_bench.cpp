#include "common.hpp"
#include <cmath>
#include <iostream>
#include <ctime>
#include <omp.h>

// Sequential implementation (adapted from your code, using heap arrays to avoid stack overflow)
MetricasSecuencial ejecutar_secuencial() {
    // allocate on heap
    const int NI = imax + 1;
    const int NK = kmax + 1;
    double *phi = new double[(size_t)NI * NK];
    double *phin = new double[(size_t)NI * NK];

    auto idx = [&](int i, int k) -> size_t { return (size_t)i * NK + k; };

    double dx = 1.0 / kmax;
    double dy = 1.0 / imax;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double dx2i = 1.0 / dx2;
    double dy2i = 1.0 / dy2;
    double dt = _min(dx2, dy2) / 4.0;

    // init
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

    int it;
    double t_start = omp_get_wtime();
    for (it = 1; it <= itmax; ++it) {
        double dphimax = 0.0;
        for (int k = 1; k < kmax; ++k) {
            for (int i = 1; i < imax; ++i) {
                double dphi =
                    (phi[idx(i+1,k)] + phi[idx(i-1,k)] - 2.0 * phi[idx(i,k)]) * dy2i +
                    (phi[idx(i,k+1)] + phi[idx(i,k-1)] - 2.0 * phi[idx(i,k)]) * dx2i;
                dphi *= dt;
                if (dphi > dphimax) dphimax = dphi;
                phin[idx(i,k)] = phi[idx(i,k)] + dphi;
            }
        }
        for (int k = 1; k < kmax; ++k)
            for (int i = 1; i < imax; ++i)
                phi[idx(i,k)] = phin[idx(i,k)];

        if (dphimax < eps) break;
    }
    double t_end = omp_get_wtime();

    MetricasSecuencial m;
    m.tiempo_total = t_end - t_start;
    m.iteraciones = it;
    long long flops = calcular_flops(imax, kmax, it);
    m.gflops = (double)flops / (m.tiempo_total * 1e9);

    delete [] phi;
    delete [] phin;
    return m;
}

int main(int argc, char** argv) {
    std::cout << "SEQ benchmark starting (this may take a while)...\n";
    MetricasSecuencial seq = ejecutar_secuencial();
    std::cout << "SEQ done: time=" << seq.tiempo_total << "s, it=" << seq.iteraciones << ", gflops=" << seq.gflops << "\n";

    export_seq_csv("plots/seq_results.csv", seq);
    std::cout << "Wrote plots/seq_results.csv\n";
    return 0;
}
