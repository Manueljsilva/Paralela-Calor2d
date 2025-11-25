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
TARGETS = secuencial version1 version2

# Archivos fuente (en SRCDIR)
SECUENCIAL_SRC = $(SRCDIR)/secuencial.cpp
VERSION1_SRC    = $(SRCDIR)/version1.cpp
VERSION2_SRC    = $(SRCDIR)/version2.cpp

# Rutas de salida (binarios en OUTDIR)
SECUENCIAL_BIN = $(OUTDIR)/secuencial
VERSION1_BIN    = $(OUTDIR)/version1
VERSION2_BIN    = $(OUTDIR)/version2

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

# Make targets remain the same names but depend on out/<name>
secuencial: $(SECUENCIAL_BIN)

version1: $(VERSION1_BIN)

version2: $(VERSION2_BIN)

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

# Benchmark - comparar todas las versiones
benchmark: secuencial version1 version2
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
	@printf "%b\n" "$(VERDE)========================================$(RESET)"

# Limpiar ejecutables (removes binaries from out/)
clean:
	@printf "%b\n" "$(AZUL)Limpiando archivos compilados...$(RESET)"
	rm -f $(SECUENCIAL_BIN) $(VERSION1_BIN) $(VERSION2_BIN)
	@printf "%b\n" "$(VERDE)✓ Limpieza completada$(RESET)"

# Limpiar y recompilar
rebuild: clean all

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
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)Ejecución:$(RESET)"
	@printf "%b\n" "  make run-secuencial   - Ejecutar versión secuencial"
	@printf "%b\n" "  make run-v1           - Ejecutar OpenMP con 4 threads"
	@printf "%b\n" "  make run-v1-1         - Ejecutar OpenMP con 1 thread"
	@printf "%b\n" "  make run-v1-2         - Ejecutar OpenMP con 2 threads"
	@printf "%b\n" "  make run-v1-4         - Ejecutar OpenMP con 4 threads"
	@printf "%b\n" "  make run-v1-8         - Ejecutar OpenMP con 8 threads"
	@printf "%b\n" "  make run-v2           - Ejecutar MPI con 4 procesos"
	@printf "%b\n" "  make run-v2-1         - Ejecutar MPI con 1 proceso"
	@printf "%b\n" "  make run-v2-2         - Ejecutar MPI con 2 procesos"
	@printf "%b\n" "  make run-v2-4         - Ejecutar MPI con 4 procesos"
	@printf "%b\n" "  make run-v2-8         - Ejecutar MPI con 8 procesos"
	@printf "%b\n" ""
	@printf "%b\n" "$(AZUL)Utilidades:$(RESET)"
	@printf "%b\n" "  make benchmark        - Ejecutar todas las versiones y comparar"
	@printf "%b\n" "  make clean            - Eliminar ejecutables"
	@printf "%b\n" "  make rebuild          - Limpiar y recompilar"
	@printf "%b\n" "  make help             - Mostrar esta ayuda"
	@printf "%b\n" ""

.PHONY: all clean rebuild run-secuencial run-v1 run-v1-1 run-v1-2 run-v1-4 run-v1-8 run-v2 run-v2-1 run-v2-2 run-v2-4 run-v2-8 benchmark help secuencial version1 version2
