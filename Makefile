# Makefile para proyecto de Ecuación de Calor 2D
# Compilación de versiones secuencial y paralelas

# Compiladores
CXX = g++
MPICXX = mpic++

# Flags de compilación (sin optimización para pruebas de Khipu)
CXXFLAGS = -Wall
OMPFLAGS = -fopenmp
MPIFLAGS = 

# Ejecutables
TARGETS = secuencial version1 version2 version2_2d version3

# Archivos fuente
SECUENCIAL_SRC = secuencial.cpp
VERSION1_SRC = version1.cpp
VERSION2_SRC = version2.cpp
VERSION2_2D_SRC = version2_2d.cpp
VERSION3_SRC = version3.cpp

# Colores para output
VERDE = \033[0;32m
AZUL = \033[0;34m
RESET = \033[0m

# Regla por defecto
all: $(TARGETS)
	@echo "$(VERDE)✓ Compilación completada$(RESET)"
	@echo "$(AZUL)Ejecutables disponibles:$(RESET)"
	@echo "  - secuencial    (versión secuencial)"
	@echo "  - version1      (OpenMP)"
	@echo "  - version2      (MPI 1D por filas - pruebas locales)"
	@echo "  - version2_2d   (MPI 2D por bloques - Khipu)"
	@echo "  - version3      (Benchmark completo: Seq + MPI + OMP)"

# Compilar versión secuencial
secuencial: $(SECUENCIAL_SRC)
	@echo "$(AZUL)Compilando versión secuencial...$(RESET)"
	$(CXX) $(CXXFLAGS) -o $@ $<
	@echo "$(VERDE)✓ secuencial compilado$(RESET)"

# Compilar versión 1 (OpenMP)
version1: $(VERSION1_SRC)
	@echo "$(AZUL)Compilando versión 1 (OpenMP)...$(RESET)"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<
	@echo "$(VERDE)✓ version1 compilado$(RESET)"

# Compilar versión 2 (MPI 1D)
version2: $(VERSION2_SRC)
	@echo "$(AZUL)Compilando versión 2 (MPI 1D por filas)...$(RESET)"
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $<
	@echo "$(VERDE)✓ version2 compilado$(RESET)"

# Compilar versión 2D (MPI 2D por bloques - para Khipu)
version2_2d: $(VERSION2_2D_SRC)
	@echo "$(AZUL)Compilando versión 2D (MPI 2D por bloques)...$(RESET)"
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $<
	@echo "$(VERDE)✓ version2_2d compilado$(RESET)"

# Compilar versión 3 (Benchmark completo: Seq + MPI + OpenMP)
version3: $(VERSION3_SRC)
	@echo "$(AZUL)Compilando versión 3 (Benchmark completo)...$(RESET)"
	$(MPICXX) $(CXXFLAGS) $(OMPFLAGS) $(MPIFLAGS) -o $@ $<
	@echo "$(VERDE)✓ version3 compilado$(RESET)"

# Ejecutar versión secuencial
run-secuencial: secuencial
	@echo "$(AZUL)Ejecutando versión secuencial...$(RESET)"
	@./secuencial

# Ejecutar versión 1 con diferentes números de threads
run-v1: version1
	@echo "$(AZUL)Ejecutando versión 1 (OpenMP) con 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 ./version1

run-v1-1: version1
	@echo "$(AZUL)Ejecutando versión 1 con 1 thread...$(RESET)"
	@OMP_NUM_THREADS=1 ./version1

run-v1-2: version1
	@echo "$(AZUL)Ejecutando versión 1 con 2 threads...$(RESET)"
	@OMP_NUM_THREADS=2 ./version1

run-v1-4: version1
	@echo "$(AZUL)Ejecutando versión 1 con 4 threads...$(RESET)"
	@OMP_NUM_THREADS=4 ./version1

run-v1-8: version1
	@echo "$(AZUL)Ejecutando versión 1 con 8 threads...$(RESET)"
	@OMP_NUM_THREADS=8 ./version1

# Ejecutar versión 2 con diferentes números de procesos
run-v2: version2
	@echo "$(AZUL)Ejecutando versión 2 (MPI) con 4 procesos...$(RESET)"
	@mpirun -np 4 ./version2

run-v2-1: version2
	@echo "$(AZUL)Ejecutando versión 2 con 1 proceso...$(RESET)"
	@mpirun -np 1 ./version2

run-v2-2: version2
	@echo "$(AZUL)Ejecutando versión 2 con 2 procesos...$(RESET)"
	@mpirun -np 2 ./version2

run-v2-4: version2
	@echo "$(AZUL)Ejecutando versión 2 con 4 procesos...$(RESET)"
	@mpirun -np 4 ./version2

run-v2-8: version2
	@echo "$(AZUL)Ejecutando versión 2 con 8 procesos...$(RESET)"
	@mpirun -np 8 ./version2

# Ejecutar versión 2D con diferentes números de procesos (para Khipu)
run-v2d-4: version2_2d
	@echo "$(AZUL)Ejecutando versión 2D con 4 procesos...$(RESET)"
	@mpirun -np 4 ./version2_2d

run-v2d-8: version2_2d
	@echo "$(AZUL)Ejecutando versión 2D con 8 procesos...$(RESET)"
	@mpirun -np 8 ./version2_2d

run-v2d-16: version2_2d
	@echo "$(AZUL)Ejecutando versión 2D con 16 procesos...$(RESET)"
	@mpirun -np 16 ./version2_2d

run-v2d-32: version2_2d
	@echo "$(AZUL)Ejecutando versión 2D con 32 procesos...$(RESET)"
	@mpirun -np 32 ./version2_2d

# Ejecutar versión 3 (Benchmark completo) con diferentes configuraciones
run-v3-1: version3
	@echo "$(AZUL)Ejecutando benchmark con 1 proceso MPI...$(RESET)"
	@mpirun -np 1 ./version3

run-v3-2: version3
	@echo "$(AZUL)Ejecutando benchmark con 2 procesos MPI...$(RESET)"
	@mpirun -np 2 ./version3

run-v3-4: version3
	@echo "$(AZUL)Ejecutando benchmark con 4 procesos MPI...$(RESET)"
	@mpirun -np 4 ./version3

run-v3-8: version3
	@echo "$(AZUL)Ejecutando benchmark con 8 procesos MPI...$(RESET)"
	@mpirun -np 8 ./version3

# Benchmark completo - generar CSV con múltiples configuraciones
benchmark-csv: version3
	@echo "$(VERDE)========================================$(RESET)"
	@echo "$(VERDE)  GENERANDO DATOS PARA GRÁFICAS$(RESET)"
	@echo "$(VERDE)========================================$(RESET)"
	@rm -f resultados_benchmark.csv
	@echo ""
	@echo "$(AZUL)Ejecutando con 1 proceso...$(RESET)"
	@mpirun -np 1 ./version3
	@echo ""
	@echo "$(AZUL)Ejecutando con 2 procesos...$(RESET)"
	@mpirun -np 2 ./version3
	@echo ""
	@echo "$(AZUL)Ejecutando con 4 procesos...$(RESET)"
	@mpirun -np 4 ./version3
	@echo ""
	@echo "$(AZUL)Ejecutando con 8 procesos...$(RESET)"
	@mpirun -np 8 ./version3
	@echo ""
	@echo "$(VERDE)✓ CSV generado: resultados_benchmark.csv$(RESET)"

# Benchmark - comparar todas las versiones
benchmark: secuencial version1 version2
	@echo "$(VERDE)========================================$(RESET)"
	@echo "$(VERDE)  BENCHMARK - Ecuación de Calor 2D$(RESET)"
	@echo "$(VERDE)========================================$(RESET)"
	@echo ""
	@echo "$(AZUL)1. Versión Secuencial:$(RESET)"
	@./secuencial
	@echo ""
	@echo "$(AZUL)2. Versión OpenMP (1 thread):$(RESET)"
	@OMP_NUM_THREADS=1 ./version1
	@echo ""
	@echo "$(AZUL)3. Versión OpenMP (2 threads):$(RESET)"
	@OMP_NUM_THREADS=2 ./version1
	@echo ""
	@echo "$(AZUL)4. Versión OpenMP (4 threads):$(RESET)"
	@OMP_NUM_THREADS=4 ./version1
	@echo ""
	@echo "$(AZUL)5. Versión OpenMP (8 threads):$(RESET)"
	@OMP_NUM_THREADS=8 ./version1
	@echo ""
	@echo "$(AZUL)6. Versión MPI (2 procesos):$(RESET)"
	@mpirun -np 2 ./version2
	@echo ""
	@echo "$(AZUL)7. Versión MPI (4 procesos):$(RESET)"
	@mpirun -np 4 ./version2
	@echo ""
	@echo "$(AZUL)8. Versión MPI (8 procesos):$(RESET)"
	@mpirun -np 8 ./version2
	@echo ""
	@echo "$(VERDE)========================================$(RESET)"

# Limpiar ejecutables
clean:
	@echo "$(AZUL)Limpiando archivos compilados...$(RESET)"
	rm -f $(TARGETS) resultados_benchmark.csv
	@echo "$(VERDE)✓ Limpieza completada$(RESET)"

# Compilar todas las versiones incluyendo 2D
all-with-2d: secuencial version1 version2 version2_2d version3
	@echo "$(VERDE)✓ Compilación completa (incluyendo versión 2D)$(RESET)"

# Limpiar y recompilar
rebuild: clean all

# Ayuda
help:
	@echo "$(VERDE)========================================$(RESET)"
	@echo "$(VERDE)  Makefile - Ecuación de Calor 2D$(RESET)"
	@echo "$(VERDE)========================================$(RESET)"
	@echo ""
	@echo "$(AZUL)Objetivos disponibles:$(RESET)"
	@echo "  make                  - Compilar todas las versiones"
	@echo "  make all              - Compilar todas las versiones"
	@echo "  make secuencial       - Compilar solo versión secuencial"
	@echo "  make version1         - Compilar solo versión OpenMP"
	@echo "  make version2         - Compilar solo versión MPI"
	@echo "  make version3         - Compilar solo benchmark completo"
	@echo ""
	@echo "$(AZUL)Ejecución:$(RESET)"
	@echo "  make run-secuencial   - Ejecutar versión secuencial"
	@echo "  make run-v1           - Ejecutar OpenMP con 4 threads"
	@echo "  make run-v1-1         - Ejecutar OpenMP con 1 thread"
	@echo "  make run-v1-2         - Ejecutar OpenMP con 2 threads"
	@echo "  make run-v1-4         - Ejecutar OpenMP con 4 threads"
	@echo "  make run-v1-8         - Ejecutar OpenMP con 8 threads"
	@echo "  make run-v2           - Ejecutar MPI con 4 procesos"
	@echo "  make run-v2-1         - Ejecutar MPI con 1 proceso"
	@echo "  make run-v2-2         - Ejecutar MPI con 2 procesos"
	@echo "  make run-v2-4         - Ejecutar MPI con 4 procesos"
	@echo "  make run-v2-8         - Ejecutar MPI con 8 procesos"
	@echo "  make run-v3-1         - Ejecutar benchmark con 1 proceso"
	@echo "  make run-v3-2         - Ejecutar benchmark con 2 procesos"
	@echo "  make run-v3-4         - Ejecutar benchmark con 4 procesos"
	@echo "  make run-v3-8         - Ejecutar benchmark con 8 procesos"
	@echo ""
	@echo "$(AZUL)Utilidades:$(RESET)"
	@echo "  make benchmark        - Ejecutar todas las versiones y comparar"
	@echo "  make benchmark-csv    - Generar CSV con datos 1,2,4,8 procesos"
	@echo "  make clean            - Eliminar ejecutables"
	@echo "  make rebuild          - Limpiar y recompilar"
	@echo "  make help             - Mostrar esta ayuda"
	@echo ""

.PHONY: all clean rebuild run-secuencial run-v1 run-v1-1 run-v1-2 run-v1-4 run-v1-8 run-v2 run-v2-1 run-v2-2 run-v2-4 run-v2-8 run-v3-1 run-v3-2 run-v3-4 run-v3-8 benchmark benchmark-csv help
