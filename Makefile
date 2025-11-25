# Makefile para proyecto de Ecuación de Calor 2D
# Compilación de versiones secuencial y paralelas

# Compiladores
CXX = g++
MPICXX = mpicxx

# Flags de compilación
CXXFLAGS = -std=c++20 -Wall
OMPFLAGS = -fopenmp
MPIFLAGS =

# Directorios
SRCDIR = src
OUTDIR = out

# Ejecutables (targets kept as before; binaries placed in OUTDIR)
TARGETS = secuencial version1 version2 version3

# Archivos fuente (en SRCDIR)
SECUENCIAL_SRC = $(SRCDIR)/secuencial.cpp
VERSION1_SRC    = $(SRCDIR)/version1.cpp
VERSION2_SRC    = $(SRCDIR)/version2.cpp
VERSION3_SRC    = $(SRCDIR)/version3.cpp

# Rutas de salida (binarios en OUTDIR)
SECUENCIAL_BIN = $(OUTDIR)/secuencial
VERSION1_BIN    = $(OUTDIR)/version1
VERSION2_BIN    = $(OUTDIR)/version2
VERSION3_BIN    = $(OUTDIR)/version3

# Colores para output
VERDE = \033[0;32m
AZUL = \033[0;34m
RESET = \033[0m

# Regla por defecto: keep original target names, make them depend on binaries in out/
all: $(TARGETS)
	@printf "%b\n" "$(VERDE)✓ Compilación completada$(RESET)"
	@printf "%b\n" "$(AZUL)Ejecutables disponibles:$(RESET)"
	@printf "%b\n" "  - secuencial    (versión secuencial)"
	@printf "%b\n" "  - version1      (OpenMP)"
	@printf "%b\n" "  - version2      (MPI no-bloqueante)"
	@printf "%b\n" "  - version3      (MPI + OpenMP híbrido)"

# Make targets remain the same names but depend on out/<name>
secuencial: $(SECUENCIAL_BIN)

version1: $(VERSION1_BIN)

version2: $(VERSION2_BIN)

version3: $(VERSION3_BIN)

# Build rules (produce out/<name>; do NOT create OUTDIR — assume it exists)
$(SECUENCIAL_BIN): $(SECUENCIAL_SRC)
	@printf "%b\n" "$(AZUL)Compilando versión secuencial...$(RESET)"
	$(CXX) $(CXXFLAGS) -o $@ $<
	@printf "%b\n" "$(VERDE)✓ secuencial compilado$(RESET)"

$(VERSION1_BIN): $(VERSION1_SRC)
	@printf "%b\n" "$(AZUL)Compilando versión 1 (OpenMP)...$(RESET)"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<
	@printf "%b\n" "$(VERDE)✓ version1 compilado$(RESET)"

$(VERSION2_BIN): $(VERSION2_SRC)
	@printf "%b\n" "$(AZUL)Compilando versión 2 (MPI)...$(RESET)"
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $<
	@printf "%b\n" "$(VERDE)✓ version2 compilado$(RESET)"

$(VERSION3_BIN): $(VERSION3_SRC)
	@printf "%b\n" "$(AZUL)Compilando versión 3 (MPI + OpenMP híbrido)...$(RESET)"
	$(MPICXX) $(CXXFLAGS) $(OMPFLAGS) $(MPIFLAGS) -o $@ $<
	@printf "%b\n" "$(VERDE)✓ version3 compilado$(RESET)"

# Ejecutar versión secuencial
run-secuencial: secuencial
	@printf "%b\n" "$(AZUL)Ejecutando versión secuencial...$(RESET)"
	@./$(SECUENCIAL_BIN)

# Ejecutar versión 1 con diferentes números de threads
run-v1: version1
	@printf "%b\n" "$(AZUL)Ejecutando versión 1 (OpenMP) con 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 ./$(VERSION1_BIN)

run-v1-1: version1
	@printf "%b\n" "$(AZUL)Ejecutando versión 1 con 1 thread...$(RESET)"
	@OMP_NUM_THREADS=1 ./$(VERSION1_BIN)

run-v1-2: version1
	@printf "%b\n" "$(AZUL)Ejecutando versión 1 con 2 threads...$(RESET)"
	@OMP_NUM_THREADS=2 ./$(VERSION1_BIN)

run-v1-4: version1
	@printf "%b\n" "$(AZUL)Ejecutando versión 1 con 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 ./$(VERSION1_BIN)

run-v1-8: version1
	@printf "%b\n" "$(AZUL)Ejecutando versión 1 con 8 threads...$(RESET)"
	@OMP_NUM_THREADS=8 ./$(VERSION1_BIN)

# Ejecutar versión 2 con diferentes números de procesos
run-v2: version2
	@printf "%b\n" "$(AZUL)Ejecutando versión 2 (MPI) con 4 procesos...$(RESET)"
	@mpirun -np 4 ./$(VERSION2_BIN)

run-v2-1: version2
	@printf "%b\n" "$(AZUL)Ejecutando versión 2 con 1 proceso...$(RESET)"
	@mpirun -np 1 ./$(VERSION2_BIN)

run-v2-2: version2
	@printf "%b\n" "$(AZUL)Ejecutando versión 2 con 2 procesos...$(RESET)"
	@mpirun -np 2 ./$(VERSION2_BIN)

run-v2-4: version2
	@printf "%b\n" "$(AZUL)Ejecutando versión 2 con 4 procesos...$(RESET)"
	@mpirun -np 4 ./$(VERSION2_BIN)

run-v2-8: version2
	@printf "%b\n" "$(AZUL)Ejecutando versión 2 con 8 procesos...$(RESET)"
	@mpirun -np 8 ./$(VERSION2_BIN)

# Ejecutar versión 3 (hybrid MPI + OpenMP)
# Default: 4 MPI processes, 4 OpenMP threads per process
run-v3: version3
	@printf "%b\n" "$(AZUL)Ejecutando versión 3 (MPI+OpenMP) con 4 procesos y 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 4 ./$(VERSION3_BIN)

# Change only OpenMP threads (process count kept at 4)
run-v3-1: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 1 thread (per rank)...$(RESET)"
	@OMP_NUM_THREADS=1 mpirun -np 4 ./$(VERSION3_BIN)

run-v3-2: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 2 threads (per rank)...$(RESET)"
	@OMP_NUM_THREADS=2 mpirun -np 4 ./$(VERSION3_BIN)

run-v3-4: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 4 threads (per rank)...$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 4 ./$(VERSION3_BIN)

run-v3-8: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 8 threads (per rank)...$(RESET)"
	@OMP_NUM_THREADS=8 mpirun -np 4 ./$(VERSION3_BIN)

# Change only MPI process count (threads kept at 4)
run-v3-p1: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 1 proceso y 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 1 ./$(VERSION3_BIN)

run-v3-p2: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 2 procesos y 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 2 ./$(VERSION3_BIN)

run-v3-p4: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 4 procesos y 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 4 ./$(VERSION3_BIN)

run-v3-p8: version3
	@printf "%b\n" "$(AZUL)Ejecutando version3 con 8 procesos y 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 8 ./$(VERSION3_BIN)

# Benchmark - comparar todas las versiones
benchmark: secuencial version1 version2 version3
	@printf "%b\n" "$(VERDE)========================================$(RESET)"
	@printf "%b\n" "$(VERDE)  BENCHMARK - Ecuación de Calor 2D$(RESET)"
	@printf "%b\n" "$(VERDE)========================================$(RESET)"
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)1. Versión Secuencial:$(RESET)"
	@./$(SECUENCIAL_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)2. Versión OpenMP (1 thread):$(RESET)"
	@OMP_NUM_THREADS=1 ./$(VERSION1_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)3. Versión OpenMP (2 threads):$(RESET)"
	@OMP_NUM_THREADS=2 ./$(VERSION1_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)4. Versión OpenMP (4 threads):$(RESET)"
	@OMP_NUM_THREADS=4 ./$(VERSION1_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)5. Versión OpenMP (8 threads):$(RESET)"
	@OMP_NUM_THREADS=8 ./$(VERSION1_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)6. Versión MPI (2 procesos):$(RESET)"
	@mpirun -np 2 ./$(VERSION2_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)7. Versión MPI (4 procesos):$(RESET)"
	@mpirun -np 4 ./$(VERSION2_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)8. Versión MPI (8 procesos):$(RESET)"
	@mpirun -np 8 ./$(VERSION2_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)9. Versión Híbrida (MPI+OpenMP) - 2 processes x 4 threads:$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 2 ./$(VERSION3_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)10. Versión Híbrida (MPI+OpenMP) - 4 processes x 4 threads:$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 4 ./$(VERSION3_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)11. Versión Híbrida (MPI+OpenMP) - 8 processes x 4 threads:$(RESET)"
	@OMP_NUM_THREADS=4 mpirun -np 8 ./$(VERSION3_BIN)
	@printf "%b\n" ""
	@printf "%b\n" "$(VERDE)========================================$(RESET)"

# Limpiar ejecutables (removes binaries from out/)
clean:
	@printf "%b\n" "$(AZUL)Limpiando archivos compilados...$(RESET)"
	rm -f $(SECUENCIAL_BIN) $(VERSION1_BIN) $(VERSION2_BIN) $(VERSION3_BIN)
	@printf "%b\n" "$(VERDE)✓ Limpieza completada$(RESET)"

# Limpiar y recompilar
rebuild: clean all

# ---------------- Slurm submission helpers ----------------
# Paths to job scripts (edit or override on the make command-line)
JOB_DIR = jobs
JOB_OMP = $(JOB_DIR)/job_omp.sbatch
JOB_MPI = $(JOB_DIR)/job_mpi.sbatch
JOB_HYBRID = $(JOB_DIR)/job_hybrid.sbatch

# Ensure slurm-out exists for sbatch outputs
.PHONY: _ensure_slurm_out
_ensure_slurm_out:
	@mkdir -p slurm-out

# Submit OpenMP job
.PHONY: submit-omp
submit-omp: _ensure_slurm_out
	@printf "%b\n" "$(AZUL)Submitting OpenMP job ($(JOB_OMP))...$(RESET)"
	@sbatch $(JOB_OMP)

# Submit MPI job
.PHONY: submit-mpi
submit-mpi: _ensure_slurm_out
	@printf "%b\n" "$(AZUL)Submitting MPI job ($(JOB_MPI))...$(RESET)"
	@sbatch $(JOB_MPI)

# Submit Hybrid (MPI+OpenMP) job
.PHONY: submit-hybrid
submit-hybrid: _ensure_slurm_out
	@printf "%b\n" "$(AZUL)Submitting Hybrid job ($(JOB_HYBRID))...$(RESET)"
	@sbatch $(JOB_HYBRID)

# Submit all three jobs in sequence
.PHONY: submit-all
submit-all: submit-omp submit-mpi submit-hybrid
	@printf "%b\n" "$(VERDE)All jobs submitted.$(RESET)"

# Ayuda
help:
	@printf "%b\n" "$(VERDE)========================================$(RESET)"
	@printf "%b\n" "$(VERDE)  Makefile - Ecuación de Calor 2D$(RESET)"
	@printf "%b\n" "$(VERDE)========================================$(RESET)"
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)Objetivos disponibles:$(RESET)"
	@printf "%b\n" "  make                  - Compilar todas las versiones"
	@printf "%b\n" "  make all              - Compilar todas las versiones"
	@printf "%b\n" "  make secuencial       - Compilar solo versión secuencial"
	@printf "%b\n" "  make version1         - Compilar solo versión OpenMP"
	@printf "%b\n" "  make version2         - Compilar solo versión MPI"
	@printf "%b\n" "  make version3         - Compilar solo versión MPI+OpenMP híbrido"
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)Ejecución local:$(RESET)"
	@printf "%b\n" "  make run-secuencial   - Ejecutar versión secuencial"
	@printf "%b\n" "  make run-v1           - Ejecutar OpenMP con 4 threads"
	@printf "%b\n" "  make run-v1-1         - Ejecutar OpenMP con 1 thread"
	@printf "%b\n" "  make run-v1-2         - Ejecutar OpenMP with 2 threads"
	@printf "%b\n" "  make run-v1-4         - Ejecutar OpenMP with 4 threads"
	@printf "%b\n" "  make run-v1-8         - Ejecutar OpenMP with 8 threads"
	@printf "%b\n" "  make run-v2           - Ejecutar MPI con 4 procesos"
	@printf "%b\n" "  make run-v2-1         - Ejecutar MPI con 1 proceso"
	@printf "%b\n" "  make run-v2-2         - Ejecutar MPI con 2 procesos"
	@printf "%b\n" "  make run-v2-4         - Ejecutar MPI con 4 procesos"
	@printf "%b\n" "  make run-v2-8         - Ejecutar MPI con 8 procesos"
	@printf "%b\n" "  make run-v3           - Ejecutar Híbrido (4 procesos x 4 threads)"
	@printf "%b\n" "  make run-v3-1         - Ejecutar Híbrido (4 procesos x 1 thread)"
	@printf "%b\n" "  make run-v3-2         - Ejecutar Híbrido (4 procesos x 2 threads)"
	@printf "%b\n" "  make run-v3-4         - Ejecutar Híbrido (4 procesos x 4 threads)"
	@printf "%b\n" "  make run-v3-8         - Ejecutar Híbrido (4 processes x 8 threads)"
	@printf "%b\n" "  make run-v3-p1        - Ejecutar Híbrido (1 proceso x 4 threads)"
	@printf "%b\n" "  make run-v3-p2        - Ejecutar Híbrido (2 processes x 4 threads)"
	@printf "%b\n" "  make run-v3-p4        - Ejecutar Híbrido (4 processes x 4 threads)"
	@printf "%b\n" "  make run-v3-p8        - Ejecutar Híbrido (8 processes x 4 threads)"
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)Envio de slurms:$(RESET)"
	@printf "%b\n" "  make submit-omp       - Envio de jobs/job_omp.sbatch"
	@printf "%b\n" "  make submit-mpi       - Envio de jobs/job_mpi.sbatch"
	@printf "%b\n" "  make submit-hybrid    - Envio de jobs/job_hybrid.sbatch"
	@printf "%b\n" "  make submit-all       - Envio de los tres jobs"
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)Utilidades:$(RESET)"
	@printf "%b\n" "  make benchmark        - Ejecutar todas las versiones y comparar"
	@printf "%b\n" "  make clean            - Eliminar ejecutables"
	@printf "%b\n" "  make rebuild          - Limpiar y recompilar"
	@printf "%b\n" "  make help             - Mostrar esta ayuda"
	@printf "%b\n" ""

.PHONY: all clean rebuild run-secuencial run-v1 run-v1-1 run-v1-2 run-v1-4 run-v1-8 run-v2 run-v2-1 run-v2-2 run-v2-4 run-v2-8 run-v3 run-v3-1 run-v3-2 run-v3-4 run-v3-8 run-v3-p1 run-v3-p2 run-v3-p4 run-v3-p8 benchmark help submit-omp submit-mpi submit-hybrid submit-all _ensure_slurm_out secuencial version1 version2 version3
