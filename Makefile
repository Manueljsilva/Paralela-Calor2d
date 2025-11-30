# Makefile minimal y claro para Ecuación de Calor 2D
# Compila los binarios en ./out y proporciona reglas para ejecutar.

# Compiladores
CXX := g++
MPICXX := mpic++

# Flags
CXXFLAGS := -Wall
OMPFLAGS := -fopenmp
MPIFLAGS :=

# Salida y ficheros
OUTDIR := out
SRCS := secuencial.cpp version1.cpp version2.cpp version3.cpp
TARGETS := $(OUTDIR)/secuencial $(OUTDIR)/version1 $(OUTDIR)/version2 $(OUTDIR)/version3

# Parámetro por defecto para mpirun (se puede sobreescribir: make run-v2 NP=6)
NP ?= 4

# Colores (opcionales)
VERDE := \033[0;32m
AZUL  := \033[0;34m
RESET := \033[0m

all: $(TARGETS)
	@echo "$(VERDE)✓ Compilación completada$(RESET)"
	@echo "$(AZUL)Ejecutables disponibles en ./out:$(RESET)"
	@echo "  - secuencial    (versión secuencial)"
	@echo "  - version1      (OpenMP)"
	@echo "  - version2      (MPI 1D por filas)"

	@echo "  - version3      (Benchmark completo: Seq + MPI + OMP)"

$(OUTDIR):
	mkdir -p $(OUTDIR)

### Reglas de compilación ###
$(OUTDIR)/secuencial: secuencial.cpp | $(OUTDIR)
	@echo "$(AZUL)Compilando secuencial...$(RESET)"
	$(CXX) $(CXXFLAGS) -o $@ $<

$(OUTDIR)/version1: version1.cpp | $(OUTDIR)
	@echo "$(AZUL)Compilando version1 (OpenMP)...$(RESET)"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

$(OUTDIR)/version2: version2.cpp | $(OUTDIR)
	@echo "$(AZUL)Compilando version2 (MPI 1D)...$(RESET)"
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $<



$(OUTDIR)/version3: version3.cpp | $(OUTDIR)
	@echo "$(AZUL)Compilando version3 (Benchmark)...$(RESET)"
	$(MPICXX) $(CXXFLAGS) $(OMPFLAGS) $(MPIFLAGS) -o $@ $<

### Reglas de ejecución ###
run-secuencial: $(OUTDIR)/secuencial
	@echo "$(AZUL)Ejecutando secuencial...$(RESET)"
	@./$(OUTDIR)/secuencial

run-v1: $(OUTDIR)/version1
	@echo "$(AZUL)Ejecutando version1 (OpenMP) con OMP_NUM_THREADS=$(NP)...$(RESET)"
	@OMP_NUM_THREADS=$(NP) ./$(OUTDIR)/version1

run-v1-threads: $(OUTDIR)/version1
	@echo "$(AZUL)Ejecutando version1 (OpenMP) con OMP_NUM_THREADS=$(OMP_NUM_THREADS)...$(RESET)"
	@OMP_NUM_THREADS=$(OMP_NUM_THREADS) ./$(OUTDIR)/version1

run-v2: $(OUTDIR)/version2
	@echo "$(AZUL)Ejecutando version2 (MPI) con $(NP) procesos...$(RESET)"
	@mpirun -np $(NP) ./$(OUTDIR)/version2

run-v2-quick: $(OUTDIR)/version2
	@echo "$(AZUL)Ejecutando version2 (MPI) con 4 procesos (quick)...$(RESET)"
	@mpirun -np 4 ./$(OUTDIR)/version2

## Nota: target run-v2d desactivado si no existe source version2_2d.cpp

run-v3: $(OUTDIR)/version3
	@echo "$(AZUL)Ejecutando benchmark (version3) con $(NP) procesos MPI...$(RESET)"
	@mpirun -np $(NP) ./$(OUTDIR)/version3

### Utilidades de benchmark ###
benchmark-csv: $(OUTDIR)/version3
	@echo "$(VERDE)Generando CSV: 1,2,4,8 procesos$(RESET)"
	@rm -f resultados_benchmark.csv
	@mpirun -np 1 ./$(OUTDIR)/version3
	@mpirun -np 2 ./$(OUTDIR)/version3
	@mpirun -np 4 ./$(OUTDIR)/version3
	@mpirun -np 8 ./$(OUTDIR)/version3
	@mpirun -np 16 ./$(OUTDIR)/version3
	@echo "$(VERDE)✓ CSV generado: resultados_benchmark.csv$(RESET)"

benchmark-separado: $(OUTDIR)/secuencial $(OUTDIR)/version1 $(OUTDIR)/version2
	@echo "$(VERDE)Ejecutando benchmark comparativo...$(RESET)"
	@./$(OUTDIR)/secuencial
	@OMP_NUM_THREADS=1 ./$(OUTDIR)/version1
	@OMP_NUM_THREADS=2 ./$(OUTDIR)/version1
	@OMP_NUM_THREADS=4 ./$(OUTDIR)/version1
	@mpirun -np 2 ./$(OUTDIR)/version2
	@mpirun -np 4 ./$(OUTDIR)/version2
	@mpirun -np 8 ./$(OUTDIR)/version2
	@mpirun -np 16 ./$(OUTDIR)/version2
	@echo "$(VERDE)Benchmark finalizado$(RESET)"

.PHONY: benchmark-equal
benchmark: $(OUTDIR)/version1 $(OUTDIR)/version2
	@echo "$(VERDE)Benchmark emparejado: ejecutar version3 si existe, sino OpenMP+MPI (p=2,4,8)$(RESET)"
	@for p in 2 4 8 16 32 ; do \
		echo "--- p=$$p ---"; \
		if [ -x ./$(OUTDIR)/version3 ]; then \
			echo "--- Ejecutando version3 con mpirun -np $$p y OMP_NUM_THREADS=$$p"; \
			OMP_NUM_THREADS=$$p mpirun -np $$p ./$(OUTDIR)/version3; \
		else \
			echo "--- Ejecutando OpenMP: OMP_NUM_THREADS=$$p"; \
			OMP_NUM_THREADS=$$p ./$(OUTDIR)/version1; \
			echo "--- Ejecutando MPI: mpirun -np $$p"; \
			mpirun -np $$p ./$(OUTDIR)/version2; \
		fi; \
	done
	@echo "$(VERDE)Benchmark emparejado finalizado$(RESET)"

### Limpiar ###
clean:
	@echo "$(AZUL)Limpiando ejecutables...$(RESET)"
	@rm -rf $(OUTDIR)
	@rm -f resultados_benchmark.csv

rebuild: clean all

help:
	@echo "Makefile para Ecuación de Calor 2D"
	@echo "  make           - Compila todo en ./out/"
	@echo "  make run-v2 NP=6      - Ejecuta version2 con NP procesos MPI"
	@echo "  make run-v1-threads OMP_NUM_THREADS=4 - Ejecuta version1 con X threads"
	@echo "  make benchmark  - Ejecuta comparativo básico"

.PHONY: all clean rebuild run-secuencial run-v1 run-v1-threads run-v2 run-v2-quick run-v2d run-v3 benchmark benchmark-csv help
