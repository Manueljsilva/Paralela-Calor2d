# Ecuación de Calor en 2D - Paralelización con MPI y OpenMP

Proyecto de paralelización de la ecuación de calor en 2D usando diferentes paradigmas de programación paralela (OpenMP y MPI).

## Descripción del Problema

Resolución de la ecuación de calor bidimensional usando el método de diferencias finitas con esquema explícito de Euler. Se simula la transmisión de calor en una malla de 80×80 puntos con condiciones de frontera fijas.

### Ecuación

$$\frac{\partial \phi}{\partial t} = \alpha \left(\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2}\right)$$

### Parámetros
- **Malla:** 80×80 puntos (configurable)
- **Iteraciones máximas:** 20,000
- **Criterio de convergencia:** $\epsilon = 10^{-8}$
- **Esquema:** Euler explícito (diferencias finitas)

## Versiones Implementadas

### 1. **Secuencial** (`secuencial.cpp`)
- Código base sin paralelización
- Referencia para validación y cálculo de speedup
- Tiempo base: ~0.75-1.5s (sin optimización)

### 2. **V1: OpenMP** (`version1.cpp`)
- **Paradigma:** Memoria compartida (multi-threading)
- **Directivas clave:**
  ```cpp
  #pragma omp parallel for reduction(max:dphimax) schedule(static)
  ```
- **Ventajas:**
  -  Paralelización 2D: divide dominio en grid 2D de iteraciones
  -  Mejor balanceo de carga con muchos threads
  -  Comparación justa con descomposición 2D de MPI
  -  Sin comunicación explícita
  -  Buen rendimiento en sistemas multi-core

### 3. **V2: MPI No-Bloqueante** (`version2.cpp`)
- **Paradigma:** Memoria distribuida (multi-proceso)
- **Descomposición:** Dominio 1D por filas
- **Comunicación:**
  - `MPI_Isend` / `MPI_Irecv` para intercambio de halos (ghost rows)
  - `MPI_Allreduce` para convergencia global
- **Estrategia de solapamiento:**
  1. Enviar/recibir halos (no-bloqueante)
  2. Computar puntos interiores (mientras se comunica)
  3. `MPI_Waitall` (esperar comunicaciones)
  4. Computar puntos frontera
- **Ventajas:**
  -  Escalable a múltiples nodos
  -  Solapamiento comunicación/cómputo
- **Overhead:** Comunicación aumenta con más procesos (~8-15% en grids pequeños)

## Compilación y Ejecución

### **Opción 1: Usando Makefile (Recomendado)**

### Compilar todo
```bash
make all
```

Genera los ejecutables en `./out/`:
- `out/secuencial` — versión secuencial
- `out/version1` — OpenMP
- `out/version2` — MPI 1D

### Ejecutar (ejemplos)
```bash
# Ejecutar secuencial
make run-secuencial

# Ejecutar OpenMP (ajustar número de hilos)
make run-v1 OMP_THREADS=4
# o alternativamente
OMP_NUM_THREADS=4 make run-v1

# Ejecutar MPI (ajustar número de procesos)
make run-v2 NP=4
# ejemplo rápido
make run-v2-quick
```

### Benchmarks (genera CSVs en ./plots/)
```bash
# Secuencial
make benchmark-sec

# OpenMP
make benchmark-omp

# MPI
make benchmark-mpi

# Todos
make benchmark
```

### Limpieza
```bash
make clean
```

### Ayuda
```bash
make help
```

### **Opción 2: Compilación Manual**

#### **Secuencial**
```bash
g++ -Wall -o secuencial secuencial.cpp
./secuencial
```

#### **V1: OpenMP**
```bash
g++ -Wall -fopenmp -o version1 version1.cpp

# Ejecutar con 4 threads
export OMP_NUM_THREADS=4
./version1
```

#### **V2: MPI**
```bash
mpic++ -Wall -o version2 version2.cpp

# Ejecutar con 4 procesos
mpirun -np 4 ./version2
```

## Análisis de Resultados

### **Generar Gráficas**

```bash
# 1. Ejecutar benchmarks
make benchmark

# 2. Generar gráficas con Python
python3 plots/main.py
```


**Observaciones:**
-  MPI escala mejor que OpenMP en este problema
-  Overhead de comunicación crece con más procesos
-  En Khipu con grids grandes (320×320+) se espera mejor eficiencia


## Uso en Khipu (Supercomputadora)

Para ejecutar en el cluster Khipu de la universidad:

### **1. Preparar Script SLURM**

```bash
#!/bin/bash
#SBATCH --job-name=calor2d
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --partition=normal

# Compilar con optimización
module load gcc/11.2.0 openmpi/4.1.1
mpic++ -Wall -fopenmp -o version3 version3.cpp

# Ejecutar benchmark con 32 procesos (4 nodos × 8 procesos)
mpirun -np 32 out/version3
```

### **2. Enviar Trabajos**

```bash
sbatch jobs/job.slurm
```

### **3. Escalar Grid para Khipu**

Modificar en `common.hpp`:
```cpp
#ifdef CLUSTER_BENCHMARK
    constexpr int imax = 3840;
    constexpr int kmax = 3840;
    constexpr int itmax = 20000;
    constexpr double eps = 1e-8;
// ...
```

## Validación de Resultados

Todos los códigos deben converger al mismo número de iteraciones:
```
Iteraciones: 14320 (para grid 80×80)
```

Verificar que `dphimax < eps` al final de la ejecución.


##  Notas de Diseño

### **¿Por qué comunicación no-bloqueante?**
- Permite solapar comunicación con cómputo
- Mejora rendimiento 10-30% vs bloqueante
- Calcula puntos interiores mientras esperan halos

##  Autor

- Alejandro Joel Ore Garcia
- Manuel Jesus Silva Anampa
- Frans Josep Trujillo Flores

##  Licencia

Proyecto académico - Universidad
