#!/usr/bin/env python3
"""
Script para generar gráficas de rendimiento desde resultados_benchmark.csv
Genera: Speedup vs Procesos, Eficiencia vs Procesos, GFlops vs Procesos
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Leer datos del CSV
df = pd.read_csv('resultados_benchmark.csv')

# Configuración de estilo
plt.style.use('seaborn-v0_8-darkgrid')
colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']

# Crear figura con 3 subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Análisis de Rendimiento - Ecuación de Calor 2D (80×80)', 
             fontsize=16, fontweight='bold')

# ============================================
# GRÁFICA 1: Speedup vs Procesos
# ============================================
ax1 = axes[0, 0]
ax1.plot(df['Procesos'], df['Speedup_MPI'], 'o-', 
         color=colors[0], linewidth=2, markersize=8, label='MPI')
ax1.plot(df['Procesos'], df['Speedup_OMP'], 's--', 
         color=colors[1], linewidth=2, markersize=8, label='OpenMP')

# Línea ideal (speedup lineal)
ideal = df['Procesos']
ax1.plot(df['Procesos'], ideal, ':', color='gray', 
         linewidth=1.5, label='Speedup Ideal')

ax1.set_xlabel('Número de Procesos', fontsize=11, fontweight='bold')
ax1.set_ylabel('Speedup (Tseq / Tpar)', fontsize=11, fontweight='bold')
ax1.set_title('Speedup vs Procesos', fontsize=12, fontweight='bold')
ax1.legend(loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_xticks(df['Procesos'])

# Anotar valores
for i, (p, s_mpi, s_omp) in enumerate(zip(df['Procesos'], 
                                           df['Speedup_MPI'], 
                                           df['Speedup_OMP'])):
    ax1.annotate(f'{s_mpi:.2f}x', 
                xy=(p, s_mpi), 
                xytext=(5, 5), 
                textcoords='offset points',
                fontsize=8, color=colors[0])

# ============================================
# GRÁFICA 2: Eficiencia vs Procesos
# ============================================
ax2 = axes[0, 1]
ax2.plot(df['Procesos'], df['Eficiencia_MPI'], 'o-', 
         color=colors[0], linewidth=2, markersize=8, label='MPI')
ax2.plot(df['Procesos'], df['Eficiencia_OMP'], 's--', 
         color=colors[1], linewidth=2, markersize=8, label='OpenMP')

# Línea de eficiencia ideal (100%)
ax2.axhline(y=100, color='gray', linestyle=':', linewidth=1.5, label='Ideal (100%)')

ax2.set_xlabel('Número de Procesos', fontsize=11, fontweight='bold')
ax2.set_ylabel('Eficiencia (%)', fontsize=11, fontweight='bold')
ax2.set_title('Eficiencia vs Procesos', fontsize=12, fontweight='bold')
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.set_xticks(df['Procesos'])

# Anotar valores
for i, (p, e_mpi) in enumerate(zip(df['Procesos'], df['Eficiencia_MPI'])):
    ax2.annotate(f'{e_mpi:.1f}%', 
                xy=(p, e_mpi), 
                xytext=(5, -15), 
                textcoords='offset points',
                fontsize=8, color=colors[0])

# ============================================
# GRÁFICA 3: GFlops vs Procesos
# ============================================
ax3 = axes[1, 0]
x = np.arange(len(df['Procesos']))
width = 0.25

bars1 = ax3.bar(x - width, df['GFlops_Seq'], width, 
                label='Secuencial', color=colors[2], alpha=0.8)
bars2 = ax3.bar(x, df['GFlops_MPI'], width, 
                label='MPI', color=colors[0], alpha=0.8)
bars3 = ax3.bar(x + width, df['GFlops_OMP'], width, 
                label='OpenMP', color=colors[1], alpha=0.8)

ax3.set_xlabel('Número de Procesos', fontsize=11, fontweight='bold')
ax3.set_ylabel('Rendimiento (GFlops)', fontsize=11, fontweight='bold')
ax3.set_title('Rendimiento (GFlops) vs Procesos', fontsize=12, fontweight='bold')
ax3.set_xticks(x)
ax3.set_xticklabels(df['Procesos'])
ax3.legend()
ax3.grid(True, alpha=0.3, axis='y')

# Anotar valores en barras
for bars in [bars1, bars2, bars3]:
    for bar in bars:
        height = bar.get_height()
        ax3.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', fontsize=8)

# ============================================
# GRÁFICA 4: Overhead de Comunicación MPI
# ============================================
ax4 = axes[1, 1]
ax4.plot(df['Procesos'], df['Comunicacion_%'], 'o-', 
         color=colors[3], linewidth=2, markersize=8)

ax4.set_xlabel('Número de Procesos', fontsize=11, fontweight='bold')
ax4.set_ylabel('Comunicación (%)', fontsize=11, fontweight='bold')
ax4.set_title('Overhead de Comunicación MPI', fontsize=12, fontweight='bold')
ax4.grid(True, alpha=0.3)
ax4.set_xticks(df['Procesos'])

# Anotar valores
for i, (p, comm) in enumerate(zip(df['Procesos'], df['Comunicacion_%'])):
    ax4.annotate(f'{comm:.2f}%', 
                xy=(p, comm), 
                xytext=(5, 5), 
                textcoords='offset points',
                fontsize=9, color=colors[3])

# Ajustar layout y guardar
plt.tight_layout()
plt.savefig('graficas_rendimiento.png', dpi=300, bbox_inches='tight')
print("✓ Gráfica guardada como: graficas_rendimiento.png")

# Mostrar gráfica
plt.show()

# ============================================
# TABLA RESUMEN
# ============================================
print("\n" + "="*80)
print("TABLA RESUMEN DE RESULTADOS")
print("="*80)
print(df[['Procesos', 'Speedup_MPI', 'Eficiencia_MPI', 
          'GFlops_MPI', 'Comunicacion_%']].to_string(index=False))
print("="*80)
print(f"\nTotal de iteraciones convergencia: {df['Iteraciones'].iloc[0]}")
print(f"Grid: 80×80 puntos")
