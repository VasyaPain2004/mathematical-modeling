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

V1 = FunctionSpace(mesh1, "CG", 1)
V2 = FunctionSpace(mesh2, "CG", 1)
V3 = FunctionSpace(mesh3, "CG", 1)

dx1 = Measure('dx', domain=mesh1, subdomain_data=subdomains1)
ds1 = Measure('ds', domain=mesh1, subdomain_data=boundaries1)
dx2 = Measure('dx', domain=mesh2, subdomain_data=subdomains2)
ds2 = Measure('ds', domain=mesh2, subdomain_data=boundaries2)
dx3 = Measure('dx', domain=mesh3, subdomain_data=subdomains3)
ds3 = Measure('ds', domain=mesh3, subdomain_data=boundaries3)

T_max = 100
N = 10
tau = T_max / N
dt = 0
      
g_T = Constant(10.0)
g_B = Constant(10.0)
lmd = Expression("1 + 0.1*x[0] - 0.2*x[1]", degree=1)

bc_1_V1 = DirichletBC(V1, g_T, boundaries1, 2)
bc_2_V1 = DirichletBC(V1, g_B, boundaries1, 4)
bcs1 = [bc_1_V1, bc_2_V1]

bc_1_V2 = DirichletBC(V2, g_T, boundaries2, 2)
bc_2_V2 = DirichletBC(V2, g_B, boundaries2, 4)
bcs2 = [bc_1_V2, bc_2_V2]

bc_1_V3 = DirichletBC(V3, g_T, boundaries3, 2)
bc_2_V3 = DirichletBC(V3, g_B, boundaries3, 4)
bcs3 = [bc_1_V3, bc_2_V3]

T1 = TrialFunction(V1)
v1 = TestFunction(V1)
T2 = TrialFunction(V2)
v2 = TestFunction(V2)
T3 = TrialFunction(V3)
v3 = TestFunction(V3)

T0_val = Constant(15.0)
T0_1 = interpolate(T0_val, V1)
T0_2 = interpolate(T0_val, V2)
T0_3 = interpolate(T0_val, V3)

a1 = (T1/tau)*v1*dx1 + lmd*inner(grad(T1), grad(v1)) * dx1
L1 = (T0_1/tau)*v1*dx1
a2 = (T2/tau)*v2*dx2 + lmd*inner(grad(T2), grad(v2)) * dx2
L2 = (T0_2/tau)*v2*dx2
a3 = (T3/tau)*v3*dx3 + lmd*inner(grad(T3), grad(v3)) * dx3
L3 = (T0_3/tau)*v3*dx3

T1 = Function(V1)
T2 = Function(V2)
T3 = Function(V3)
file1 = File("/mnt/c/Users/Vasya/Downloads/temperature1.pvd")
file2 = File("/mnt/c/Users/Vasya/Downloads/temperature2.pvd")
file3 = File("/mnt/c/Users/Vasya/Downloads/temperature3.pvd")

while dt < T_max:
  dt += tau
  solve(a1 == L1, T1, bcs1)
  solve(a2 == L2, T2, bcs2)
  solve(a3 == L3, T3, bcs3)
  T1.rename('T1', '0')
  T2.rename('T2', '0')
  T3.rename('T3', '0')
  file1 << T1
  file2 << T2
  file3 << T3
  T0_1.assign(T1)
  T0_2.assign(T2)
  T0_3.assign(T3)
  
T1_on_V3 = interpolate(T1, V3)
E2_t = inner(T3 - T1_on_V3, T3 - T1_on_V3)*dx3
E2_b = inner(T3, T3)*dx3
E2 = sqrt(abs(assemble(E2_t))/abs(assemble(E2_b)))

print("Error =", str(E2))