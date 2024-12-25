import dolfinx
import ufl
import math
import numpy as np

def solve_fem_diffusion_2d(domain, D_expr, rhs_expr, dirichlet_bcs, neumann_bcs=None):
    """
    Solve -div(D grad(u)) = rhs_expr in domain, with advanced boundary conditions.
    domain: a dolfinx mesh
    D_expr: UFL expression or dolfin Constant
    rhs_expr: UFL expression or dolfin Constant
    dirichlet_bcs: list of (subdomain_marker, value)
    neumann_bcs: optional list of (subdomain_marker, flux_expr)
    returns: (FunctionSpace, solution Function)
    """
    V = dolfinx.fem.FunctionSpace(domain, ("CG",1))

    # define test/trial
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # define boundary conditions
    # This code is just a placeholder approach for advanced subdomain usage.
    # You can locate dofs geometrically or with subdomain markers.
    # We'll assume user does that externally.
    # For demonstration, let's assume "dirichlet_bcs" is a list of (locate_dofs, value).
    # E.g. see main HPC code for how to do it with fem.locate_dofs_geometrical.

    bcs = []
    for (dofs, val) in dirichlet_bcs:
        bc = dolfinx.fem.dirichletbc(val, dofs, V)
        bcs.append(bc)

    # define PDE: -div(D grad(u)) = rhs_expr
    a = ufl.dot(D_expr*ufl.grad(u), ufl.grad(v))*ufl.dx
    L = rhs_expr*v*ufl.dx

    # neumann BC => add boundary integrals if needed
    # e.g. L += flux_expr * v * ufl.ds(subdomain_id)

    A = dolfinx.fem.petsc.assemble_matrix(dolfinx.fem.form(a), bcs=bcs)
    A.assemble()
    b = dolfinx.fem.petsc.assemble_vector(dolfinx.fem.form(L))
    dolfinx.fem.petsc.apply_lifting(b, [dolfinx.fem.form(a)], bcs=[bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    dolfinx.fem.petsc.set_bc(b, bcs)

    # solve
    A_mat = A
    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A_mat)
    solver.setType("preonly")
    solver.getPC().setType("lu")
    solver.solve(b, b)

    # solution
    uh = dolfinx.fem.Function(V)
    uh.vector.axpy(1,b)
    uh.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
    return (V, uh)
