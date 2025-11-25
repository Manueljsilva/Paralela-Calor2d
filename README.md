# Ecuación de Calor en 2D - Paralelización

Proyecto de paralelización de la ecuación de calor en 2D usando diferentes paradigmas de programación paralela.

## Descripción del Problema

Resolución de la ecuación de calor bidimensional usando el método de diferencias finitas con esquema explícito de Euler. Se simula la transmisión de calor en una malla de 80x80 puntos con condiciones de frontera fijas.

## Versiones Implementadas

### Cuadro Comparativo de Versiones

| Característica | Secuencial | V1: OpenMP | V2: MPI No-Bloqueante | V3: Híbrida MPI+OpenMP |
|----------------|------------|------------|----------------------|------------------------|
| **Paradigma** | Secuencial | Memoria Compartida | Memoria Distribuida | Memoria Distribuida + Compartida |
| **Tecnología** | C++ estándar | OpenMP | MPI (Isend/Irecv) | MPI + OpenMP |
| **Descomposición** | - | Ninguna (threads comparten todo) | Dominio 1D (por filas) | Dominio 1D entre procesos, threads dentro |
| **Comunicación** | - | Implícita (memoria compartida) | Explícita (intercambio de fronteras) | MPI entre procesos, compartida intra-proceso |
| **Tipo de Comunicación** | - | - | No-bloqueante (solapamiento) | No-bloqueante entre procesos |
| **Sincronización** | - | Barrera implícita en directivas | MPI_Wait, MPI_Allreduce | MPI_Wait + barreras OpenMP |
| **Escalabilidad** | - | Limitada (1 nodo, multi-core) | Alta (multi-nodo) | Muy Alta (multi-nodo + multi-core) |
| **Complejidad Implementación** | Baja | Baja | Media-Alta | Alta |
| **Balanceo de Carga** | - | Automático (schedule) | Manual (distribución filas) | Híbrido |
| **Regiones Paralelizadas** | - | Bucles espaciales (k) | Descomposición de dominio por filas | Ambos niveles |
| **Variables Críticas** | - | `dphimax` (reduction) | `dphimax` (MPI_Allreduce) | Reducción en ambos niveles |
| **Comunicación/Cómputo** | - | - | Solapado (Isend/Irecv) | Solapado multi-nivel |
| **Ideal para** | Validación | Sistemas multi-core | Clusters | Supercomputadoras |

### Detalles de Implementación

#### **Versión Secuencial (Base)**
- Código original sin modificaciones
- Sirve como referencia para validación
- Medición de tiempo base para calcular speedup

#### **V1: OpenMP (Memoria Compartida)**
***Características:***
- Paralelización del bucle externo `k` con `#pragma omp parallel for`
- `reduction(max:dphimax)` para criterio de convergencia
- Variables privadas: `i`, `dphi`
- `schedule(static)` para distribución de iteraciones
- No requiere modificación de la estructura de datos
- Sincronización automática en barreras

***Ventajas:***
- Fácil implementación (pocas líneas de código)
- No requiere comunicación explícita
- Balanceo automático de carga
- Ideal para sistemas multi-core de un solo nodo
- Buen speedup en arquitecturas modernas

***Limitaciones:***
- Limitado a un solo nodo
- Escalabilidad limitada por cores disponibles

#### **V2: MPI No-Bloqueante (Memoria Distribuida)**
***Características:***
- Descomposición 1D del dominio (por filas)
- Cada proceso calcula subdominio local
- Intercambio de filas frontera (ghost rows) con vecinos
- `MPI_Isend` / `MPI_Irecv` para comunicación asíncrona
- `MPI_Allreduce` para `dphimax` global
- Solapamiento de comunicación y cómputo

***Estrategia de Comunicación:***
1. Iniciar envíos/recepciones no-bloqueantes (`MPI_Isend`/`MPI_Irecv`)
2. Calcular puntos interiores mientras se comunican halos (solapamiento)
3. Esperar que terminen comunicaciones (`MPI_Waitall`)
4. Calcular puntos frontera que usan halos
5. Actualizar valores locales
6. Reducción global de `dphimax` con `MPI_Allreduce`

***Ventajas:***
- Escalable a múltiples nodos
- Mayor rendimiento por solapamiento
- Distribución real de memoria

***Limitaciones:***
- Mayor complejidad de código
- Requiere manejo explícito de fronteras
- Overhead de comunicación

#### **V3: Híbrida MPI + OpenMP**
***Características:***
- MPI para distribución entre procesos (nodos)
- OpenMP para paralelizar dentro de cada proceso
- Dos niveles de paralelismo
- Reduce comunicaciones MPI vs V2 puro
- Aprovecha arquitectura moderna (multi-nodo + multi-core)

***Estrategia:***
- Nivel externo (MPI): Descomposición de dominio
- Nivel interno (OpenMP): Paralelización de bucles locales
- Comunicación MPI solo entre procesos
- Memoria compartida dentro de cada proceso

***Ventajas:***
- Máxima escalabilidad
- Reduce comunicaciones vs MPI puro
- Aprovecha todo el hardware disponible
- Flexible en configuración procesos/threads

***Limitaciones:**
- Mayor complejidad
- Requiere ajuste fino de parámetros
- Debugging más difícil

## Aspectos Comunes a Todas las Versiones Paralelas

### **Regiones Paralelizadas:**

**OpenMP (V1):**
- Bucle espacial externo `k` con `#pragma omp parallel for`
- Bucle de actualización de valores

**MPI (V2):**
- Descomposición de dominio: cada proceso calcula su subdominio de filas
- No hay directivas de paralelización, el paralelismo viene de la distribución

**Híbrida (V3):**
- MPI: descomposición de dominio entre procesos
- OpenMP: paralelización de bucles dentro de cada proceso

### **Regiones NO Paralelizadas:**
- Inicialización de condiciones de frontera
- Bucle temporal externo (`it = 1 to itmax`) - **todos los procesos lo ejecutan simultáneamente**, pero tiene dependencia entre iteraciones

### **Variables Críticas:**
- `dphimax`: Requiere reducción (max) para criterio de convergencia

## Compilación y Ejecución

### Usando Makefile

```bash
# Vease make help para mas informacion
make help
```

## Métricas de Evaluación

Para cada versión se debe medir:
- **Tiempo de ejecución**
- **Speedup**: $S_p = \frac{T_{secuencial}}{T_{paralelo}}$
- **Eficiencia**: $E_p = \frac{S_p}{p} \times 100\%$
- **Validación**: Comparación de resultados numéricos

## Referencias

- Método de diferencias finitas para ecuaciones parabólicas
- Descomposición de dominio para problemas 2D
- Comunicación no-bloqueante en MPI
- Programación híbrida MPI+OpenMP
