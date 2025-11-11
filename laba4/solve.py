from dolfin import *

meshname = "mesh1"
mesh1 = Mesh(meshname + ".xml")
subdomains1 = MeshFunction("size_t", mesh1, meshname + "_physical_region.xml")
boundaries1 = MeshFunction("size_t", mesh1, meshname + "_facet_region.xml")

meshname = "mesh2"
mesh2 = Mesh(meshname + ".xml")
subdomains2 = MeshFunction("size_t", mesh2, meshname + "_physical_region.xml")
boundaries2 = MeshFunction("size_t", mesh2, meshname + "_facet_region.xml")

meshname = "mesh3"
mesh3 = Mesh(meshname + ".xml")
subdomains3 = MeshFunction("size_t", mesh3, meshname + "_physical_region.xml")
boundaries3 = MeshFunction("size_t", mesh3, meshname + "_facet_region.xml")

dx1 = Measure('dx', domain=mesh1, subdomain_data=subdomains1)
ds1 = Measure('ds', domain=mesh1, subdomain_data=boundaries1)
dx2 = Measure('dx', domain=mesh2, subdomain_data=subdomains2)
ds2 = Measure('ds', domain=mesh2, subdomain_data=boundaries2)
dx3 = Measure('dx', domain=mesh3, subdomain_data=subdomains3)
ds3 = Measure('ds', domain=mesh3, subdomain_data=boundaries3)

mu = 40
lmbda = 80

F = Expression(("sin(pi*x[0]) * cos(pi*x[1])", "exp(-x[0]*x[0] - x[1]*x[1])"), degree=2)
u_D = Expression(("sin(pi*x[0]) * cos(pi*x[1])", "sin(pi*x[0]) * cos(pi*x[1])"), degree=2)

def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)

def sigma(u):
    return lmbda*div(u)*Identity(2) + 2*mu*epsilon(u)

V1 = VectorFunctionSpace(mesh1, "CG", 1)
V2 = VectorFunctionSpace(mesh2, "CG", 1)
V3 = VectorFunctionSpace(mesh3, "CG", 1)

bc1 = DirichletBC(V1, u_D, boundaries1, 4)
bcs1 = [bc1]
bc2 = DirichletBC(V2, u_D, boundaries2, 4)
bcs2 = [bc2]
bc3 = DirichletBC(V3, u_D, boundaries3, 4)
bcs3 = [bc3]

u1 = TrialFunction(V1)
v1 = TestFunction(V1)
u2 = TrialFunction(V2)
v2 = TestFunction(V2)
u3 = TrialFunction(V3)
v3 = TestFunction(V3)

a1 = inner(sigma(u1), epsilon(v1))*dx1
L1 = inner(F, v1)*dx1
a2 = inner(sigma(u2), epsilon(v2))*dx2
L2 = inner(F, v2)*dx2
a3 = inner(sigma(u3), epsilon(v3))*dx3
L3 = inner(F, v3)*dx3

u1 = Function(V1)
u2 = Function(V2)
u3 = Function(V3)

solve(a1 == L1, u1, bcs1)
solve(a2 == L2, u2, bcs2)
solve(a3 == L3, u3, bcs3)

file1 = File("/mnt/c/Users/Vasya/Downloads/u1.pvd")
file1 << u1
file2 = File("/mnt/c/Users/Vasya/Downloads/u2.pvd")
file2 << u2
file3 = File("/mnt/c/Users/Vasya/Downloads/u3.pvd")
file3 << u3

u1_on_V3 = interpolate(u1, V3)
E2_t_1 = inner(u3 - u1_on_V3, u3 - u1_on_V3)*dx3
E2_b_1 = inner(u3, u3)*dx3
E2_1 = sqrt(abs(assemble(E2_t_1))/abs(assemble(E2_b_1)))

u2_on_V3 = interpolate(u2, V3)
E2_t_2 = inner(u3 - u2_on_V3, u3 - u2_on_V3)*dx3
E2_b_2 = inner(u3, u3)*dx3
E2_2 = sqrt(abs(assemble(E2_t_2))/abs(assemble(E2_b_2)))

print('mesh1 pogreshnost:', E2_1)
print('mesh2 pogreshnost:', E2_2)

W1 = TensorFunctionSpace(mesh1, "DG", 0)
W2 = TensorFunctionSpace(mesh2, "DG", 0)
W3 = TensorFunctionSpace(mesh3, "DG", 0)

stress1 = project(sigma(u1), V=W1)
stress2 = project(sigma(u2), V=W2)
stress3 = project(sigma(u3), V=W3)

File("/mnt/c/Users/Vasya/Downloads/s1.pvd") << stress1
File("/mnt/c/Users/Vasya/Downloads/s2.pvd") << stress2
File("/mnt/c/Users/Vasya/Downloads/s3.pvd") << stress3