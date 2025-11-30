import pandas as pd

def merge_benchmarks(sec_file: str, omp_file: str, mpi_file: str) -> pd.DataFrame:
    df_sec = pd.read_csv(sec_file)       # Tiempo,Iteraciones,GFlops
    df_omp = pd.read_csv(omp_file)       # Threads,TiempoTotal,Iteraciones,GFlops
    df_mpi = pd.read_csv(mpi_file)       # Ranks,TiempoTotal,TiempoComputo,TiempoComunicacion,Iteraciones,GFlops,Comunicacion_%

    # Only keep rows where MPI Ranks == OMP Threads
    df_merged = pd.merge(df_mpi, df_omp, left_on="Ranks", right_on="Threads", suffixes=("_mpi","_omp"))

    # Build final DataFrame in old format
    rows = []
    for _, row in df_merged.iterrows():
        procesos = int(row["Ranks"])
        threads = int(row["Threads"])
        t_seq = float(df_sec["Tiempo"].iloc[0])
        t_mpi = float(row["TiempoTotal_mpi"])
        t_omp = float(row["TiempoTotal_omp"])
        iteraciones = int(row["Iteraciones_mpi"])

        speedup_mpi = t_seq / t_mpi
        eficiencia_mpi = (speedup_mpi / procesos) * 100
        speedup_omp = t_seq / t_omp
        eficiencia_omp = (speedup_omp / threads) * 100

        gflops_seq = float(df_sec["GFlops"].iloc[0])
        gflops_mpi = float(row["GFlops_mpi"])
        gflops_omp = float(row["GFlops_omp"])
        comunicacion = float(row.get("Comunicacion_%", 0.0))

        rows.append({
            "Procesos": procesos,
            "Threads": threads,
            "T_Secuencial": t_seq,
            "T_MPI": t_mpi,
            "T_OMP": t_omp,
            "Speedup_MPI": speedup_mpi,
            "Eficiencia_MPI": eficiencia_mpi,
            "Speedup_OMP": speedup_omp,
            "Eficiencia_OMP": eficiencia_omp,
            "GFlops_Seq": gflops_seq,
            "GFlops_MPI": gflops_mpi,
            "GFlops_OMP": gflops_omp,
            "Comunicacion_%": comunicacion,
            "Iteraciones": iteraciones
        })

    df_final = pd.DataFrame(rows)
    df_final.sort_values(["Procesos", "Threads"], inplace=True)
    df_final.reset_index(drop=True, inplace=True)
    return df_final
