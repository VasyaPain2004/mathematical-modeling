from dolfin import *

mesh = Mesh("mesh.xml")
subdomains = MeshFunction("size_t", mesh, 'mesh' + "_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, 'mesh' + "_facet_region.xml")

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

mu = 1
rho = 1
nu = mu / rho
U_inf = 10.0
T_max = 0.5
dt = T_max / 36

inflow = DirichletBC(V, (U_inf, 0), boundaries, 1)
outflow_p = DirichletBC(Q, Constant(0.0), boundaries, 2)
bcu = [inflow]
bcp = [outflow_p]

u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

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

A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

parameters['krylov_solver']['nonzero_initial_guess'] = True

ufile = File("/mnt/c/Users/Vasya/Downloads/velocity.pvd")
pfile = File("/mnt/c/Users/Vasya/Downloads/pressure.pvd")

t = dt
while t < T_max + DOLFIN_EPS:

    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "bicgstab", "default")

    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    [bc.apply(p1.vector()) for bc in bcp]
    solve(A2, p1.vector(), b2, "bicgstab", prec)

    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "bicgstab", "default")

    ufile << u1
    pfile << p1

    u0.assign(u1)
    t += dt
    print(t)

