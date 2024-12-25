import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

def apply_boundary_conditions_2d(u, bc_type, bc_value):
    """
    Applies 2D boundary conditions to array u in-place:
      - dirichlet: set boundary to bc_value
      - neumann: replicate neighbor cells (zero-flux)
    """
    ny, nx = u.shape
    if bc_type == "dirichlet":
        u[0, :] = bc_value
        u[ny - 1, :] = bc_value
        u[:, 0] = bc_value
        u[:, nx - 1] = bc_value
    elif bc_type == "neumann":
        # top/bottom
        u[0, :] = u[1, :]
        u[ny - 1, :] = u[ny - 2, :]
        # left/right
        u[:, 0] = u[:, 1]
        u[:, nx - 1] = u[:, nx - 2]

def adi_diffusion_2d_variableD(
    conc, D, dt, dx, dy, boundary_type, boundary_value, reaction=None
):
    """
    Perform one ADI step for 2D diffusion with variable or constant D.

    conc: 2D array for concentration
    D: 2D array or float for diffusion
    dt, dx, dy: floats for time step and cell sizes
    boundary_type: "dirichlet" or "neumann"
    boundary_value: boundary value (float)
    reaction: 2D array for source/sink term R (optional)
    returns updated 2D array for concentration
    """
    ny, nx = conc.shape
    if isinstance(D, float) or isinstance(D, int):
        D = np.full((ny, nx), float(D))
    if reaction is None:
        reaction = np.zeros_like(conc)

    # X sweep
    conc_half = conc.copy()
    for j in range(ny):
        row = conc[j, :]
        rxn_row = reaction[j, :]
        D_row = D[j, :]
        # Build tri-di
        A_data = np.zeros((3, nx))
        B = np.zeros(nx)
        for i in range(nx):
            B[i] = row[i] + 0.5*dt*rxn_row[i]

        for i in range(nx):
            if 0 < i < nx - 1:
                D_left = 0.5*(D_row[i] + D_row[i - 1])
                D_right= 0.5*(D_row[i] + D_row[i + 1])
            else:
                D_left = D_row[i]
                D_right= D_row[i]
            alpha_left = (D_left*dt)/(2*dx*dx)
            alpha_right= (D_right*dt)/(2*dx*dx)
            main = 1 + alpha_left + alpha_right
            A_data[1, i] = main
            if i>0:
                A_data[0, i] = -alpha_left
            if i<nx-1:
                A_data[2, i] = -alpha_right

        diag_lower = A_data[0,1:]
        diag_main  = A_data[1,:]
        diag_upper = A_data[2,:-1]
        from scipy.sparse import diags
        A = diags([diag_lower, diag_main, diag_upper], [-1,0,1], shape=(nx,nx), format="csr")
        row_new = spsolve(A, B)
        conc_half[j, :] = row_new

    apply_boundary_conditions_2d(conc_half, boundary_type, boundary_value)

    # Y sweep
    conc_new = conc_half.copy()
    for i in range(nx):
        col = conc_half[:, i]
        rxn_col = reaction[:, i]
        D_col = D[:, i]
        A_data = np.zeros((3, ny))
        B = np.zeros(ny)
        for j in range(ny):
            B[j] = col[j] + 0.5*dt*rxn_col[j]

        for j in range(ny):
            if 0< j < ny - 1:
                D_down = 0.5*(D_col[j] + D_col[j - 1])
                D_up   = 0.5*(D_col[j] + D_col[j + 1])
            else:
                D_down = D_col[j]
                D_up   = D_col[j]
            alpha_down = (D_down*dt)/(2*dy*dy)
            alpha_up   = (D_up*dt)/(2*dy*dy)
            main = 1 + alpha_down + alpha_up
            A_data[1, j] = main
            if j>0:
                A_data[0, j] = -alpha_down
            if j<ny-1:
                A_data[2, j] = -alpha_up

        diag_lower = A_data[0,1:]
        diag_main  = A_data[1,:]
        diag_upper = A_data[2,:-1]
        A = diags([diag_lower, diag_main, diag_upper], [-1,0,1], shape=(ny,ny), format="csr")
        col_new = spsolve(A, B)
        conc_new[:, i] = col_new

    apply_boundary_conditions_2d(conc_new, boundary_type, boundary_value)
    return conc_new
