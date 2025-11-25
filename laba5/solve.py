from dolfin import *

# Load mesh from file
mesh = Mesh("mesh.xml")
subdomains = MeshFunction("size_t", mesh, 'mesh' + "_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, 'mesh' + "_facet_region.xml")

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.1
T = 1
mu = 0.000018
rho = 1.2
nu = mu / rho

U_inf = 10.0

# Define time-dependent pressure boundary condition
# p_in = Expression("sin(3.0*t)", t=0.0, degree=2)

# Define boundary conditions
# noslip  = DirichletBC(V, (0, 0),
#                       "on_boundary && \
#                        (x[0] < DOLFIN_EPS | x[1] < DOLFIN_EPS | \
#                        (x[0] > 0.5 - DOLFIN_EPS && x[1] > 0.5 - DOLFIN_EPS))")
# inflow  = DirichletBC(Q, p_in, "x[1] > 1.0 - DOLFIN_EPS")
# outflow = DirichletBC(Q, 0, "x[0] > 1.0 - DOLFIN_EPS")

inflow = DirichletBC(V, (U_inf, 0), boundaries, 1)
outflow = DirichletBC(V, (U_inf, 0), boundaries, 3)
bcu = [inflow, outflow]
bcp = []

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Create files for storing solution
ufile = File("/mnt/c/Users/Vasya/Downloads/velocity.pvd")
pfile = File("/mnt/c/Users/Vasya/Downloads/pressure.pvd")

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
#     p_in.t = t

    # Compute tentative velocity step
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "mumps")

    # Pressure correction
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    [bc.apply(p1.vector()) for bc in bcp]
    solve(A2, p1.vector(), b2, "mumps")

    # Velocity correction
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "mumps")

    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    t += dt