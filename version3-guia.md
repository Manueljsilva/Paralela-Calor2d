Aquí tienes el archivo Markdown listo para que lo copies y lo incluyas en tu proyecto o reporte. Explica detalladamente cómo funciona el código híbrido y cómo ejecutarlo.

-----

# Guía de Uso: Benchmark Híbrido (MPI vs OpenMP)

Este documento describe el funcionamiento y uso del archivo `comparativa.cpp`. El objetivo de este programa es unificar dos paradigmas de programación paralela (Memoria Distribuida con **MPI** y Memoria Compartida con **OpenMP**) en un solo ejecutable para medir y comparar su rendimiento bajo las mismas condiciones.

## 1\. Descripción del Código

El archivo `comparativa.cpp` no ejecuta MPI y OpenMP simultáneamente para resolver el problema cooperando, sino que realiza una **"carrera"** secuencial:

1.  **Fase MPI:** Todos los procesos inician y resuelven la ecuación de calor distribuyendo el dominio.
2.  **Sincronización:** Se espera a que todos terminen.
3.  **Fase OpenMP:** Únicamente el **Proceso 0 (Master)** resuelve el mismo problema usando hilos (threads) de OpenMP.
4.  **Reporte:** Se imprimen los tiempos de ambos para comparar.

### Estrategia de Medición de Tiempos

  * **MPI (`MPI_Wtime`):**

      * Se mide el **Tiempo Total** desde el inicio hasta el fin del bucle temporal.
      * Se mide el **Tiempo de Comunicación** acumulando los segundos transcurridos durante las llamadas a `MPI_Isend`, `MPI_Irecv`, `MPI_Waitall` y `MPI_Allreduce`.
      * El **Tiempo de Cómputo** se deriva de la resta: $Total - Comunicación$.

  * **OpenMP (`omp_get_wtime`):**

      * Al ser memoria compartida, no hay latencia de red. Se mide el tiempo de pared (wall-clock) que tarda el bloque paralelo en ejecutarse.

-----

## 2\. Compilación

Para compilar este código híbrido, necesitas un compilador que entienda tanto las librerías de MPI como las directivas de OpenMP.

**Comando genérico (Linux/WSL):**

```bash
mpic++ comparativa.cpp -o benchmark -fopenmp
```

  * `mpic++`: Invoca el compilador de C++ con los wrappers de MPI.
  * `-o benchmark`: Nombre del ejecutable de salida.
  * `-fopenmp`: Habilita el soporte para `#pragma omp`.

> **Nota:** Si usas macOS con compiladores clang, la bandera podría ser `-Xpreprocessor -fopenmp -lomp`.

-----

## 3\. Ejecución

Para realizar una comparación justa, debes controlar cuántos procesos MPI lanzas y cuántos hilos OpenMP permites.

**Ejemplo: Comparar 4 Procesos MPI contra 4 Hilos OpenMP**

```bash
export OMP_NUM_THREADS=4
mpirun -np 4 ./benchmark
```

1.  `export OMP_NUM_THREADS=4`: Define que la sección de OpenMP usará 4 hilos.
2.  `mpirun -np 4`: Inicia el programa lanzando 4 procesos MPI.

-----

## 4\. Interpretación de Resultados

El programa generará una salida similar a esta:

```text
========================================
  RESULTADOS COMPARATIVOS
========================================
1. MPI (4 Procesos):
   - Tiempo Total        :     0.045000 s
   - Tiempo Computo      :     0.015000 s
   - Tiempo Comunicacion :     0.030000 s (66.67%)

2. OpenMP (4 Hilos):
   - Tiempo Total        :     0.020000 s

----------------------------------------
GANADOR: OpenMP (2.25x mas rapido)
Nota: Para grids pequeños (80x80), la latencia de red en MPI suele dominar.
========================================
```

### Puntos Clave para el Análisis:

1.  **Overhead de Comunicación (MPI):**

      * Observa el porcentaje de comunicación. En problemas pequeños (como una matriz 80x80), es probable que el tiempo de envío de mensajes sea mayor que el tiempo que tarda el CPU en calcular. Esto se llama "dominio de latencia".
      * Si el `% de comunicación` es alto (\>50%), significa que el problema es muy pequeño para la cantidad de procesos MPI usados.

2.  **Eficiencia de OpenMP:**

      * En un solo nodo (una sola computadora), OpenMP suele ser más rápido para tamaños de problema pequeños o medianos porque no tiene que empaquetar datos ni enviarlos por red; solo lee/escribe en la memoria RAM compartida.

3.  **Escalabilidad:**

      * Para demostrar la fuerza de MPI, intenta aumentar `imax` y `kmax` (por ejemplo a 1000 o 2000) en el código. Verás cómo el porcentaje de cómputo sube y MPI se vuelve más eficiente.