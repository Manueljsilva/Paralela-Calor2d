# Makefile para Ecuación de Calor 2D con colores más consistentes
# Compila los binarios en ./out y proporciona reglas para ejecutar.

# Compiladores
CXX := g++
MPICXX := mpic++

# Flags
CXXFLAGS := -Wall
OMPFLAGS := -fopenmp
MPIFLAGS :=

# Salida y ficheros
FILES := secuencial version1 version2 version3

SRCDIR := src
OUTDIR := out
TARGETS := $(addprefix $(OUTDIR)/,$(FILES))
SRCFILES := $(addprefix $(SRCDIR)/,$(addsuffix .cpp,$(FILES)))

# Parámetro por defecto para mpirun
NP ?= 4

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
	@printf "%b\n" "$(AZUL)  - version3      (Benchmark completo: Seq + MPI + OMP)$(RESET)"

$(OUTDIR):
	@printf "%b\n" "$(AZUL)Creando directorio $(OUTDIR)...$(RESET)"
	@mkdir -p $(OUTDIR)

### Reglas de compilación ###
$(OUTDIR)/secuencial: src/secuencial.cpp | $(OUTDIR)
	@printf "%b\n" "$(AZUL)Compilando secuencial...$(RESET)"
	@$(CXX) $(CXXFLAGS) -o $@ $<

$(OUTDIR)/version1: src/version1.cpp | $(OUTDIR)
	@printf "%b\n" "$(AZUL)Compilando version1 (OpenMP)...$(RESET)"
	@$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

$(OUTDIR)/version2: src/version2.cpp | $(OUTDIR)
	@printf "%b\n" "$(AZUL)Compilando version2 (MPI 1D)...$(RESET)"
	@$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $<

$(OUTDIR)/version3: src/version3.cpp | $(OUTDIR)
	@printf "%b\n" "$(AZUL)Compilando version3 (Benchmark)...$(RESET)"
	@$(MPICXX) $(CXXFLAGS) $(OMPFLAGS) $(MPIFLAGS) -o $@ $<

### Reglas de ejecución ###
run-secuencial: $(OUTDIR)/secuencial
	@printf "%b\n" "$(AZUL)Ejecutando secuencial...$(RESET)"
	@./$(OUTDIR)/secuencial

run-v1: $(OUTDIR)/version1
	@printf "%b\n" "$(AZUL)Ejecutando version1 (OpenMP) con OMP_NUM_THREADS=$(NP)...$(RESET)"
	@OMP_NUM_THREADS=$(NP) ./$(OUTDIR)/version1

run-v1-threads: $(OUTDIR)/version1
	@printf "%b\n" "$(AZUL)Ejecutando version1 (OpenMP) con OMP_NUM_THREADS=$(OMP_NUM_THREADS)...$(RESET)"
	@OMP_NUM_THREADS=$(OMP_NUM_THREADS) ./$(OUTDIR)/version1

run-v2: $(OUTDIR)/version2
	@printf "%b\n" "$(AZUL)Ejecutando version2 (MPI) con $(NP) procesos...$(RESET)"
	@mpirun -np $(NP) ./$(OUTDIR)/version2

run-v2-quick: $(OUTDIR)/version2
	@printf "%b\n" "$(AZUL)Ejecutando version2 (MPI) con 4 procesos (quick)...$(RESET)"
	@mpirun -np 4 ./$(OUTDIR)/version2

run-v3: $(OUTDIR)/version3
	@printf "%b\n" "$(AZUL)Ejecutando benchmark (version3) con $(NP) procesos MPI...$(RESET)"
	@mpirun -np $(NP) ./$(OUTDIR)/version3

### Utilidades de benchmark ###
benchmark-csv: $(OUTDIR)/version3
	@printf "%b\n" "$(VERDE)Generando CSV para 1,2,4,8,16 procesos...$(RESET)"
	@rm -f plots/resultados_benchmark.csv
	@rm -f plots/graficas_rendimiento.png
	@for np in 1 2 4 8 16; do \
		printf "%b%d%b\n" "$(AZUL)--- p=" $$np " ---$(RESET)"; \
		mpirun -np $$np ./$(OUTDIR)/version3; \
	done
	@printf "%b\n" "$(VERDE)✓ CSV generado: resultados_benchmark.csv$(RESET)"

benchmark-separado: $(OUTDIR)/secuencial $(OUTDIR)/version1 $(OUTDIR)/version2
	@printf "%b\n" "$(VERDE)Ejecutando benchmark comparativo...$(RESET)"
	@printf "%b\n" "$(AZUL)Ejecutando secuencial$(RESET)"
	@./$(OUTDIR)/secuencial
	@for threads in 1 2 4; do \
		printf "%b%d%b\n" "$(AZUL)--- threads=" $$threads " ---$(RESET)"; \
		OMP_NUM_THREADS=$$threads ./$(OUTDIR)/version1; \
	done
	@for np in 2 4 8 16; do \
		printf "%b%d%b\n" "$(AZUL)--- p=" $$np " ---$(RESET)"; \
		mpirun -np $$np ./$(OUTDIR)/version2; \
	done
	@printf "%b\n" "$(VERDE)Benchmark finalizado$(RESET)"

benchmark: $(OUTDIR)/version1 $(OUTDIR)/version2
	@printf "%b\n" "$(VERDE)Benchmark emparejado: ejecutar version3 si existe, sino OpenMP+MPI (p=2,4,8,16,32)$(RESET)"
	@for p in 2 4 8 16 32 ; do \
		printf "%b%d%b\n" "$(AZUL)--- p=" $$p " ---$(RESET)"; \
		if [ -x ./$(OUTDIR)/version3 ]; then \
			printf "%b%d%b\n" "$(AZUL)Ejecutando version3 con mpirun -np " $$p " y OMP_NUM_THREADS=" $$p "$(RESET)"; \
			OMP_NUM_THREADS=$$p mpirun -np $$p ./$(OUTDIR)/version3; \
		else \
			printf "%b%d%b\n" "$(AZUL)Ejecutando OpenMP: OMP_NUM_THREADS=" $$p "$(RESET)"; \
			OMP_NUM_THREADS=$$p ./$(OUTDIR)/version1; \
			printf "%b%d%b\n" "$(AZUL)Ejecutando MPI: mpirun -np " $$p "$(RESET)"; \
			mpirun -np $$p ./$(OUTDIR)/version2; \
		fi; \
	done
	@printf "%b\n" "$(VERDE)Benchmark emparejado finalizado$(RESET)"

### Limpiar ###
clean:
	@printf "%b\n" "$(AZUL)Limpiando ejecutables...$(RESET)"
	@rm -rf $(OUTDIR)
	@rm -f plots/resultados_benchmark.csv
	@rm -f plots/graficas_rendimiento.png

rebuild: clean all

help:
	@printf "%b\n" "$(AZUL)Makefile para Ecuación de Calor 2D$(RESET)"
	@printf "%b\n" "$(AZUL)  make                  - Compila todo en ./out/$(RESET)"
	@printf "%b\n" "$(AZUL)  make run-v2 NP=6      - Ejecuta version2 con NP procesos MPI$(RESET)"
	@printf "%b\n" "$(AZUL)  make run-v1-threads OMP_NUM_THREADS=4 - Ejecuta version1 con X threads$(RESET)"
	@printf "%b\n" "$(AZUL)  make benchmark        - Ejecuta comparativo básico$(RESET)"

.PHONY: all clean rebuild run-secuencial run-v1 run-v1-threads run-v2 run-v2-quick run-v3 benchmark benchmark-csv help
