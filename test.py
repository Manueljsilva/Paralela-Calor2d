#!/usr/bin/env python3
"""
Estimación determinista del número de iteraciones para el esquema:
    phi_new = phi + dt * L(phi)
donde L es el operador discreto que usa tu programa C.

Construye L (sparse), forma T = I + dt*L, estima rho(T) con power method
(inicialización determinista: vector de unos), y calcula k mediante la fórmula asintótica.
"""

import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def build_L(nx, ny, dx, dy):
    """
    Construye la matriz L (NxN sparse) que representa el operador discreto
    usado en tu código para los puntos interiores:
      lap = (east+west-2*center)*dy2i + (north+south-2*center)*dx2i
    Notación: nx = imax-1 (nº puntos interiores en x), ny = kmax-1 (en y).
    Ordering: lexicographic: index = j*nx + i  (i=0..nx-1, j=0..ny-1)
    """
    N = nx * ny
    dx2i = 1.0 / (dx*dx)
    dy2i = 1.0 / (dy*dy)

    # diagonals
    main_diag = -2.0 * (dy2i + dx2i) * np.ones(N)

    # Off-diagonals: east/west (x-direction) have coefficient dy2i
    # north/south (y-direction) have coefficient dx2i

    # create three diagonals for x-direction for each row
    offsets = []
    data = []

    # main diagonal
    offsets.append(0)
    data.append(main_diag)

    # east (i+1) and west (i-1) inside same row: offset +1 and -1 but skip row boundaries
    east = dy2i * np.ones(N-1)
    west = dy2i * np.ones(N-1)
    # must zero entries that cross row boundaries (between i = nx-1 and i+1 = 0)
    for row in range(1, ny):
        pos = row*nx - 1
        east[pos] = 0.0
        west[pos] = 0.0

    offsets.append(1)
    data.append(east)
    offsets.append(-1)
    data.append(west)

    # north (j+1) offset = +nx, south (j-1) offset = -nx, coefficient dx2i
    north = dx2i * np.ones(N - nx)
    south = dx2i * np.ones(N - nx)
    offsets.append(nx)
    data.append(north)
    offsets.append(-nx)
    data.append(south)

    L = sp.diags(data, offsets, shape=(N, N), format='csr')
    return L

def power_method_spectral_radius(matvec, N, maxiter=1000, tol=1e-12):
    """
    Power method to estimate spectral radius (max |lambda|) of a linear operator
    given by matvec(v) -> vector of length N.
    Deterministic: initial vector is ones.
    """

    v = np.ones(N, dtype=float)
    v /= np.linalg.norm(v)
    prev_norm = 0.0
    for it in range(maxiter):
        w = matvec(v)
        w_norm = np.linalg.norm(w)
        if w_norm == 0:
            return 0.0
        v = w / w_norm
        if abs(w_norm - prev_norm) < tol * max(1.0, abs(prev_norm)):
            return w_norm
        prev_norm = w_norm
    return prev_norm

def estimate_iterations(rho, eps, e0_norm=1.0):
    if rho <= 0:
        return 0
    if rho >= 1.0:
        return math.inf
    k = math.log(eps / e0_norm) / math.log(rho)
    return int(math.ceil(k))

def main(imax=80, kmax=80, itmax=20000, eps=1e-8, e0_norm=1.0):
    # match your C code: dx = 1.0/kmax, dy = 1.0/imax
    dx = 1.0 / kmax
    dy = 1.0 / imax
    nx = imax - 1   # interior points in x (i=1..imax-1)
    ny = kmax - 1   # interior points in y (k=1..kmax-1)
    N = nx * ny

    # dt as in C:
    dx2 = dx*dx
    dy2 = dy*dy
    dt = min(dx2, dy2) / 4.0

    print("Parameters:")
    print(f" imax={imax}, kmax={kmax}, nx={nx}, ny={ny}, N={N}")
    print(f" dx={dx:.6g}, dy={dy:.6g}, dt={dt:.6g}, eps={eps}, e0_norm={e0_norm}")

    # Build L (sparse)
    L = build_L(nx, ny, dx, dy)

    # Define matvec for T = I + dt * L, deterministic
    def matvec_T(v):
        return v + dt * (L.dot(v))

    # Estimate spectral radius of T deterministically (initial vector ones)
    rho_T = power_method_spectral_radius(matvec_T, N, maxiter=1000, tol=1e-14)
    print(f"Estimated spectral radius rho(T) = {rho_T:.12g}")

    # If rho_T >= 1 -> no decay (or unstable); otherwise estimate iterations
    if rho_T >= 1.0:
        print("Warning: rho(T) >= 1.0 -> no convergence (or marginal).")
        est_iters = math.inf
    else:
        est_iters = estimate_iterations(rho_T, eps, e0_norm=e0_norm)

    print("Estimated iterations to reduce error from e0_norm to eps:")
    if est_iters == math.inf:
        print("  infinite (no decay / unstable)")
    else:
        print(f"  {est_iters} iterations (asymptotic estimate)")

    return {
        "imax": imax, "kmax": kmax, "nx": nx, "ny": ny, "N": N,
        "dx": dx, "dy": dy, "dt": dt,
        "rho_T": rho_T, "estimated_iterations": est_iters
    }

if __name__ == "__main__":
    res = main(
        imax=500,
        kmax=500,
        itmax=400000,
        eps=1e-10,
        e0_norm=1.0
    )
