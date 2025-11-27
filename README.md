# Ecuación de Calor en 2D - Paralelización con MPI y OpenMP

Proyecto de paralelización de la ecuación de calor en 2D usando diferentes paradigmas de programación paralela (OpenMP y MPI).

## 📋 Descripción del Problema

Resolución de la ecuación de calor bidimensional usando el método de diferencias finitas con esquema explícito de Euler. Se simula la transmisión de calor en una malla de 80×80 puntos con condiciones de frontera fijas.

### Ecuación

$$\frac{\partial \phi}{\partial t} = \alpha \left(\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2}\right)$$

### Parámetros
- **Malla:** 80×80 puntos (configurable)
- **Iteraciones máximas:** 20,000
- **Criterio de convergencia:** $\epsilon = 10^{-8}$
- **Esquema:** Euler explícito (diferencias finitas)

## 🚀 Versiones Implementadas

### 1. **Secuencial** (`secuencial.cpp`)
- Código base sin paralelización
- Referencia para validación y cálculo de speedup
- Tiempo base: ~0.75-1.5s (sin optimización)

### 2. **V1: OpenMP** (`version1.cpp`)
- **Paradigma:** Memoria compartida (multi-threading)
- **Paralelización:** Bucle espacial externo (k-loop)
- **Directivas clave:**
  ```cpp
  #pragma omp parallel for private(i, dphi) reduction(max:dphimax) schedule(static)
  ```
- **Ventajas:**
  - ✅ Fácil implementación
  - ✅ Sin comunicación explícita
  - ✅ Buen rendimiento en sistemas multi-core
- **Escalabilidad:** Limitada a un nodo

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
  - ✅ Escalable a múltiples nodos
  - ✅ Solapamiento comunicación/cómputo
- **Overhead:** Comunicación aumenta con más procesos (~8-15% en grids pequeños)

### 4. **V3: Benchmark Completo** (`version3.cpp`)
- **Propósito:** Herramienta de análisis comparativo
- **Ejecuta:** Secuencial → MPI → OpenMP (en orden)
- **Métricas generadas:**
  - Speedup: $S_p = T_{seq} / T_{par}$
  - Eficiencia: $E_p = (S_p / p) \times 100\%$
  - GFlops (rendimiento computacional)
  - Overhead de comunicación (%)
- **Salida:** CSV para análisis (`resultados_benchmark.csv`)
- **Gráficas:** Script Python incluido (`graficar.py`)

## 🛠️ Compilación y Ejecución

### **Opción 1: Usando Makefile (Recomendado)**

```bash
# Compilar todas las versiones
make all

# Compilar versiones individuales
make secuencial
make version1      # OpenMP
make version2      # MPI
make version3      # Benchmark completo

# Ejecutar con configuraciones predefinidas
make run-v1-4      # OpenMP con 4 threads
make run-v2-4      # MPI con 4 procesos
make run-v3-4      # Benchmark con 4 procesos MPI

# Generar datos para gráficas (1, 2, 4 procesos)
make benchmark-csv

# Limpiar ejecutables y datos
make clean

# Ver todas las opciones disponibles
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

#### **V3: Benchmark (MPI + OpenMP)**
```bash
# Compilador MPI con flag OpenMP
mpic++ -Wall -fopenmp -o version3 version3.cpp

# Ejecutar con 4 procesos MPI
# (OpenMP usará threads disponibles automáticamente)
mpirun -np 4 ./version3
```

**Nota:** Sin `-O3` para pruebas locales. En Khipu (producción) agregar `-O3` para optimización.

## Análisis de Resultados

### **Generar Gráficas**

```bash
# 1. Ejecutar benchmarks
make benchmark-csv

# 2. Generar gráficas con Python
python3 graficar.py
```

Esto genera `graficas_rendimiento.png` con:
- **Speedup vs Procesos** (MPI y OpenMP vs Ideal)
- **Eficiencia vs Procesos** 
- **GFlops vs Procesos**
- **Overhead de Comunicación MPI**

### **Resultados Típicos (Grid 80×80, sin -O3)**

| Procesos | Speedup MPI | Eficiencia MPI | GFlops MPI | Comunicación % |
|----------|-------------|----------------|------------|----------------|
| 1        | 0.81x       | 81.2%          | 1.83       | 0.17%          |
| 2        | 1.51x       | 75.7%          | 3.00       | 8.56%          |
| 4        | 2.80x       | 70.0%          | 4.93       | 15.40%         |

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
mpic++ -O3 -Wall -fopenmp -o version3 version3.cpp

# Ejecutar benchmark con 32 procesos (4 nodos × 8 procesos)
mpirun -np 32 ./version3
```

### **2. Enviar Trabajo**

```bash
sbatch job_calor2d.sh
```

### **3. Escalar Grid para Khipu**

Modificar en los archivos `.cpp`:
```cpp
const int imax = 320;  // Grid más grande
const int kmax = 320;
```

Recomendaciones:
- **80×80:** Pruebas locales (1-4 procesos)
- **320×320:** Khipu con 16-32 procesos
- **640×640:** Khipu con 64-128 procesos

## 🔬 Validación de Resultados

Todos los códigos deben converger al mismo número de iteraciones:
```
Iteraciones: 14320 (para grid 80×80)
```

Verificar que `dphimax < eps` al final de la ejecución.

##  Referencias y Documentación

- **OpenMP:** [openmp.org](https://www.openmp.org/)
- **MPI:** [mpi-forum.org](https://www.mpi-forum.org/)
- **Diferencias Finitas:** Método explícito de Euler para EDPs parabólicas
- **Descomposición de Dominio:** Estrategia 1D por filas (row-major en C++)

##  Notas de Diseño

### **¿Por qué comunicación no-bloqueante?**
- Permite solapar comunicación con cómputo
- Mejora rendimiento 10-30% vs bloqueante
- Calcula puntos interiores mientras esperan halos

##  Autor

- Manuel Jesus Silva Anampa 
- ...
- ...   

##  Licencia

Proyecto académico - Universidad
