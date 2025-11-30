#ifndef COMMON_HPP
#define COMMON_HPP

#include <fstream>
#include <iomanip>
#include <string>
#include <filesystem>
#include <iostream>

// Macros de utilidad
#define _min(A, B) ((A) < (B) ? (A) : (B))
#define _max(A, B) ((A) > (B) ? (A) : (B))

// Constantes del problema
//#define CLUSTER_BENCHMARK

#ifdef CLUSTER_BENCHMARK
    constexpr int imax = 3840;
    constexpr int kmax = 3840;
    constexpr int itmax = 20000;
    constexpr double eps = 1e-10;
#else
    constexpr int imax = 80;
    constexpr int kmax = 80;
    constexpr int itmax = 20000;
    constexpr double eps = 1.0e-08;
#endif

// Metric structs
struct MetricasMPI {
    double tiempo_total = 0.0;
    double tiempo_computo = 0.0;
    double tiempo_comunicacion = 0.0;
    int iteraciones = 0;
    double gflops = 0.0;
    double porcentaje_comunicacion = 0.0;
};

struct MetricasOMP {
    double tiempo_total = 0.0;
    int iteraciones = 0;
    double gflops = 0.0;
};

struct MetricasSecuencial {
    double tiempo_total = 0.0;
    int iteraciones = 0;
    double gflops = 0.0;
};

// FLOP model
inline long long calcular_flops(int filas, int columnas, int iteraciones) {
    long long puntos_internos = (long long)(filas - 1) * (columnas - 1);
    return puntos_internos * (long long)iteraciones * 19LL;
}

// Ensure directory exists
inline void ensure_dir_for(const std::string &filepath) {
    namespace fs = std::filesystem;
    fs::path p(filepath);
    fs::path dir = p.parent_path();
    if (!dir.empty() && !fs::exists(dir)) {
        try { fs::create_directories(dir); }
        catch (const std::exception &e) {
            std::cerr << "WARNING: failed to create directory " << dir << " (" << e.what() << ")\n";
        }
    }
}

// Exporters (each program uses only the relevant one)

// seq: Tiempo, Iteraciones, GFlops
inline void export_seq_csv(const std::string &filepath, const MetricasSecuencial &seq) {
    ensure_dir_for(filepath);
    std::ofstream f(filepath, std::ios::app);
    if (!f.is_open()) { std::cerr << "ERROR: can't open " << filepath << "\n"; return; }
    f.seekp(0, std::ios::end);
    if (f.tellp() == 0) {
        f << "Tiempo,Iteraciones,GFlops\n";
    }
    f << std::fixed << std::setprecision(6)
      << seq.tiempo_total << "," << seq.iteraciones << "," << seq.gflops << "\n";
    f.close();
}

// mpi: Ranks, TiempoTotal, TiempoComputo, TiempoComunicacion, Iteraciones, GFlops, Comunicacion_%
inline void export_mpi_csv(const std::string &filepath, int num_ranks, const MetricasMPI &mpi) {
    ensure_dir_for(filepath);
    std::ofstream f(filepath, std::ios::app);
    if (!f.is_open()) { std::cerr << "ERROR: can't open " << filepath << "\n"; return; }
    f.seekp(0, std::ios::end);
    if (f.tellp() == 0) {
        f << "Ranks,TiempoTotal,TiempoComputo,TiempoComunicacion,Iteraciones,GFlops,Comunicacion_%\n";
    }
    f << std::fixed << std::setprecision(6)
      << num_ranks << ","
      << mpi.tiempo_total << ","
      << mpi.tiempo_computo << ","
      << mpi.tiempo_comunicacion << ","
      << mpi.iteraciones << ","
      << mpi.gflops << ","
      << mpi.porcentaje_comunicacion << "\n";
    f.close();
}

// omp: Threads, TiempoTotal, Iteraciones, GFlops
inline void export_omp_csv(const std::string &filepath, int threads, const MetricasOMP &omp) {
    ensure_dir_for(filepath);
    std::ofstream f(filepath, std::ios::app);
    if (!f.is_open()) { std::cerr << "ERROR: can't open " << filepath << "\n"; return; }
    f.seekp(0, std::ios::end);
    if (f.tellp() == 0) {
        f << "Threads,TiempoTotal,Iteraciones,GFlops\n";
    }
    f << std::fixed << std::setprecision(6)
      << threads << ","
      << omp.tiempo_total << ","
      << omp.iteraciones << ","
      << omp.gflops << "\n";
    f.close();
}

#endif // COMMON_HPP
