# Makefile para Ecuación de Calor 2D

# Compiladores
CXX      := g++
MPICXX   := mpic++

# Flags
CXXFLAGS := -Wall
OMPFLAGS := -fopenmp
MPIFLAGS :=

# Salida y ficheros
FILES   := secuencial version1 version2

SRCDIR  := src
OUTDIR  := out
TARGETS := $(addprefix $(OUTDIR)/,$(FILES))
SRCFILES := $(addprefix $(SRCDIR)/,$(addsuffix .cpp,$(FILES)))

# Parámetro por defecto para MPI
NP ?= 4
# Parámetro por defecto para OpenMP
OMP_THREADS ?= 4

# Colores ANSI
VERDE := \033[0;32m
AZUL  := \033[0;34m
RESET := \033[0m

all: $(TARGETS)
	@printf "%b\n" "$(VERDE)✓ Compilación completada$(RESET)"
	@printf "%b\n" "$(AZUL)Ejecutables disponibles en ./out:$(RESET)"
	@printf "%b\n" "$(AZUL)  - secuencial    (versión secuencial)$(RESET)"
	@printf "%b\n" "$(AZUL)  - version1      (OpenMP)$(RESET)"
	@printf "%b\n" "$(AZUL)  - version2      (MPI 1D por filas)$(RESET)"

$(OUTDIR):
	@printf "%b\n" "$(AZUL)Creando directorio $(OUTDIR)...$(RESET)"
	@mkdir -p $(OUTDIR)

### Reglas de compilación ###
$(OUTDIR)/secuencial: $(SRCDIR)/secuencial.cpp | $(OUTDIR)
	@printf "%b\n" "$(AZUL)Compilando secuencial...$(RESET)"
	@$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

$(OUTDIR)/version1: $(SRCDIR)/version1.cpp | $(OUTDIR)
	@printf "%b\n" "$(AZUL)Compilando version1 (OpenMP)...$(RESET)"
	@$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

$(OUTDIR)/version2: $(SRCDIR)/version2.cpp | $(OUTDIR)
	@printf "%b\n" "$(AZUL)Compilando version2 (MPI 1D)...$(RESET)"
	@$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $<

### Reglas de ejecución ###
run-secuencial: $(OUTDIR)/secuencial
	@printf "%b\n" "$(AZUL)Ejecutando secuencial...$(RESET)"
	@./$(OUTDIR)/secuencial

run-v1: $(OUTDIR)/version1
	@printf "%b\n" "$(AZUL)Ejecutando version1 (OpenMP) con OMP_NUM_THREADS=$(OMP_THREADS)...$(RESET)"
	@OMP_NUM_THREADS=$(OMP_THREADS) ./$(OUTDIR)/version1

run-v1-threads: $(OUTDIR)/version1
	@printf "%b\n" "$(AZUL)Ejecutando version1 (OpenMP) con OMP_NUM_THREADS=$(OMP_NUM_THREADS)...$(RESET)"
	@OMP_NUM_THREADS=$(OMP_NUM_THREADS) ./$(OUTDIR)/version1

run-v2: $(OUTDIR)/version2
	@printf "%b\n" "$(AZUL)Ejecutando version2 (MPI) con $(NP) procesos...$(RESET)"
	@mpirun -np $(NP) ./$(OUTDIR)/version2

run-v2-quick: $(OUTDIR)/version2
	@printf "%b\n" "$(AZUL)Ejecutando version2 (MPI) con 4 procesos (quick)...$(RESET)"
	@mpirun -np 4 ./$(OUTDIR)/version2

### Utilidades de benchmark ###
benchmark-sec: $(OUTDIR)/secuencial
	@printf "%b\n" "$(VERDE)Benchmark secuencial:$(RESET)"
	@rm -f plots/sec_results.csv
	@mkdir -p plots
	@./$(OUTDIR)/secuencial
	@printf "%b\n" "$(VERDE)✓ CSV generado: plots/sec_results.csv$(RESET)"

benchmark-omp: $(OUTDIR)/version1
	@printf "%b\n" "$(VERDE)Benchmark OMP:$(RESET)"
	@rm -f plots/omp_results.csv
	@mkdir -p plots
	@for nt in 1 2 4 8 16; do \
		printf "%b%d%b\n" "$(AZUL)--- Threads = " $$nt " ---$(RESET)"; \
		OMP_NUM_THREADS=$$nt ./$(OUTDIR)/version1; \
	done
	@printf "%b\n" "$(VERDE)✓ CSV generado: plots/omp_results.csv$(RESET)"

benchmark-mpi: $(OUTDIR)/version2
	@printf "%b\n" "$(VERDE)Benchmark MPI:$(RESET)"
	@rm -f plots/mpi_results.csv
	@mkdir -p plots
	@for np in 1 2 4 8 16; do \
		printf "%b%d%b\n" "$(AZUL)--- Processes = " $$np " ---$(RESET)"; \
		mpirun -np $$np ./$(OUTDIR)/version2; \
	done
	@printf "%b\n" "$(VERDE)✓ CSV generado: plots/mpi_results.csv$(RESET)"

benchmark: benchmark-sec benchmark-omp benchmark-mpi
	@printf "%b\n" "$(VERDE)✓ Todos los benchmarks ejecutados$(RESET)"

benchmark-show: benchmark
	@printf "%b\n" "$(VERDE)Generando y mostrando grafica:$(RESET)"
	@python3 plots/main.py

### Limpiar ###
clean:
	@printf "%b\n" "$(AZUL)Limpiando ejecutables...$(RESET)"
	@rm -rf $(OUTDIR)
	@rm -f plots/sec_results.csv
	@rm -f plots/omp_results.csv
	@rm -f plots/mpi_results.csv
	@rm -f plots/graficas_rendimiento.png

rebuild: clean all

help:
	@printf "%b\n" "$(AZUL)Makefile para Ecuación de Calor 2D$(RESET)"
	@printf "%b\n" "$(AZUL)Compilación y ejecución de las diferentes versiones y benchmarks$(RESET)"
	@printf "%b\n" "$(AZUL)=========================================$(RESET)"
	@printf "%b\n" "$(AZUL)  make                  - Compila todos los ejecutables en ./out/$(RESET)"
	@printf "%b\n" "$(AZUL)  make clean            - Elimina ejecutables y CSVs/graficas generadas$(RESET)"
	@printf "%b\n" "$(AZUL)  make rebuild          - Limpia y recompila todo$(RESET)"
	@printf "%b\n" "$(AZUL)=========================================$(RESET)"
	@printf "%b\n" "$(AZUL)Ejecución de versiones individuales$(RESET)"
	@printf "%b\n" "$(AZUL)  make run-secuencial          - Ejecuta la versión secuencial$(RESET)"
	@printf "%b\n" "$(AZUL)  make run-v1 NP=4             - Ejecuta version1 (OpenMP) con NP threads (usa OMP_NUM_THREADS)\$(RESET)"
	@printf "%b\n" "$(AZUL)  make run-v1-threads          - Ejecuta version1 (OpenMP) con OMP_NUM_THREADS definido$(RESET)"
	@printf "%b\n" "$(AZUL)  make run-v2 NP=4             - Ejecuta version2 (MPI) con NP procesos$(RESET)"
	@printf "%b\n" "$(AZUL)  make run-v2-quick            - Ejecuta version2 (MPI) con 4 procesos (rápido para pruebas)$(RESET)"
	@printf "%b\n" "$(AZUL)=========================================$(RESET)"
	@printf "%b\n" "$(AZUL)Reglas de benchmark y generación de CSV$(RESET)"
	@printf "%b\n" "$(AZUL)  make benchmark-sec           - Ejecuta benchmark secuencial y genera plots/sec_results.csv$(RESET)"
	@printf "%b\n" "$(AZUL)  make benchmark-omp           - Ejecuta benchmark OpenMP (varios NP) y genera plots/omp_results.csv$(RESET)"
	@printf "%b\n" "$(AZUL)  make benchmark-mpi           - Ejecuta benchmark MPI (varios NP) y genera plots/mpi_results.csv$(RESET)"
	@printf "%b\n" "$(AZUL)  make benchmark               - Ejecuta benchmark-sec, benchmark-omp y benchmark-mpi$(RESET)"
	@printf "%b\n" "$(AZUL)  make benchmark-show          - Ejecuta todos los benchmarks y luego genera y muestra las gráficas con Python$(RESET)"
	@printf "%b\n" "$(AZUL)=========================================$(RESET)"
	@printf "%b\n" "$(AZUL)Opciones y variables comunes$(RESET)"
	@printf "%b\n" "$(AZUL)  NP=4                         - Número de procesos MPI para ejecución MPI$(RESET)"
	@printf "%b\n" "$(AZUL)  OMP_NUM_THREADS=4            - Número de threads OpenMP para version1$(RESET)"

.PHONY: all clean rebuild run-secuencial run-v1 run-v1-threads run-v2 run-v2-quick\
        benchmark benchmark-show help
