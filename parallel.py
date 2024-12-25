"""
parallel.py

Examples of accelerating PDE loops or bacterial updates with Numba/CuPy.
"""

import numpy as np

try:
    import cupy as cp
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False

from numba import njit, prange

@njit(parallel=True)
def accelerate_bacteria_update(
    non_resistant, resistant, dt, base_growth_rate, antibiotic_eff, carrying_capacity
):
    ny, nx = non_resistant.shape
    consumption = np.zeros_like(non_resistant)
    for j in prange(ny):
        for i in prange(nx):
            nr = non_resistant[j, i]
            r  = resistant[j, i]
            total = nr + r
            if total < 1e-12:
                continue
            logistic_factor = (1 - total/carrying_capacity)
            nr_growth = base_growth_rate*(1 - antibiotic_eff)*logistic_factor
            r_growth  = base_growth_rate*logistic_factor
            d_nr = dt*nr*nr_growth
            d_r  = dt*r*r_growth
            mut_rate = 0.0001 + 0.001*antibiotic_eff
            mutated = nr*mut_rate*dt

            nr_new = nr + d_nr - mutated
            r_new  = r + d_r + mutated
            if nr_new<0: nr_new=0
            if r_new<0:  r_new=0

            non_resistant[j, i] = nr_new
            resistant[j, i] = r_new
            consumption[j, i] = d_nr + d_r

    return consumption

def gpu_diffusion_step_1d(u, D, dt, dx):
    """
    Example 1D explicit diffusion step using CuPy on GPU.
    """
    if not CUPY_AVAILABLE:
        raise ImportError("cupy not installed. GPU diffusion not available.")
    u_gpu = cp.asarray(u)
    D_gpu = cp.asarray(D)
    u_new = cp.zeros_like(u_gpu)
    i = cp.arange(1, u_gpu.size-1)
    u_new[i] = u_gpu[i] + D_gpu[i]*dt*(u_gpu[i-1] - 2*u_gpu[i] + u_gpu[i+1])/(dx*dx)
    return cp.asnumpy(u_new)
