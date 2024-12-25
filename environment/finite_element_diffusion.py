import dolfinx
import ufl
import math
import numpy as np

from dolfinx import fem
from petsc4py import PETSc

def solve_fem_diffusion_2d(domain, D_expr, rhs_expr, dirichlet_bcs, neumann_bcs=None):
    """
    Solve -div(D grad(u)) = rhs_expr in a 2D domain with advanced boundary conditions.
    domain: a dolfinx mesh
    D_expr: ufl expression or dolfinx.fem.Constant for diffusion
    rhs_expr: ufl expression or dolfinx.fem.Constant for source
    dirichlet_bcs: list of (locate_dofs, value)
    neumann_bcs: optional list of (marker_id, flux_expr) for BCs
    returns: (V, uh)
    """
    V = fem.FunctionSpace(domain, ("CG", 1))

    # trial/test
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Dirichlet BC
    bcs = []
    for (dof_list, val) in dirichlet_bcs:
        bc = fem.dirichletbc(val, dof_list, V)
        bcs.append(bc)

    # PDE
    a = ufl.dot(D_expr*ufl.grad(u), ufl.grad(v))*ufl.dx
    L = rhs_expr*v*ufl.dx

    # if neumann_bcs:
    #   for (marker_id, flux_expr) in neumann_bcs:
    #       L += flux_expr * v * ufl.ds(marker_id)

    a_form = fem.form(a)
    L_form = fem.form(L)
    A = fem.petsc.assemble_matrix(a_form, bcs=bcs)
    A.assemble()
    b = fem.petsc.assemble_vector(L_form)
    fem.petsc.apply_lifting(b, [a_form], bcs=[bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    fem.petsc.set_bc(b, bcs)

    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType("preonly")
    solver.getPC().setType("lu")
    x = b.duplicate()
    solver.solve(b, x)

    uh = fem.Function(V)
    uh.x.array[:] = x.array
    uh.x.scatter_forward()
    return V, uh

