import time
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

Q1 = FunctionSpace(mesh1, "CG", 1)
Q2 = FunctionSpace(mesh2, "CG", 1)
Q3 = FunctionSpace(mesh3, "CG", 1)
DAY = 86400
Month = 30*DAY
tau = 0.5*DAY

sou = Constant(-10.0)
rhoL = Constant(6.0*1.0e7)
croS = Constant(2.0*1.0e6) #c-pho-
croL = Constant(2.5*1.0e6) #c+pho+
kS = Constant(1.8) #lambda-
kL = Constant(2.0) #lambda+
delta = Constant(1.0)

initT = Constant(2.0)
T0_1 = interpolate(initT, Q1)
T0_2 = interpolate(initT, Q2)
T0_3 = interpolate(initT, Q3)

phi_1 = conditional(gt(T0_1, delta), 1.0, conditional(lt(T0_1, -delta), 0.0,
(T0_1+delta)/(delta+delta)))
dphi_1 = conditional(gt(T0_1, delta), 0.0, conditional(lt(T0_1, -delta),
0.0, 1.0/(delta+delta)))

phi_2 = conditional(gt(T0_2, delta), 1.0, conditional(lt(T0_2, -delta), 0.0,
(T0_2+delta)/(delta+delta)))
dphi_2 = conditional(gt(T0_2, delta), 0.0, conditional(lt(T0_2, -delta),
0.0, 1.0/(delta+delta)))

phi_3 = conditional(gt(T0_3, delta), 1.0, conditional(lt(T0_3, -delta), 0.0,
(T0_3+delta)/(delta+delta)))
dphi_3 = conditional(gt(T0_3, delta), 0.0, conditional(lt(T0_3, -delta),
0.0, 1.0/(delta+delta)))

k1 = kS + phi_1*(kL - kS) #lambda(phi)
cro1 = croS + phi_1*(croL - croS) #alpha(phi)

k2 = kS + phi_2*(kL - kS) #lambda(phi)
cro2 = croS + phi_2*(croL - croS) #alpha(phi)

k3 = kS + phi_3*(kL - kS) #lambda(phi)
cro3 = croS + phi_3*(croL - croS) #alpha(phi)

w1 = Function(Q1)
w2 = Function(Q2)
w3 = Function(Q3)

bc1 = DirichletBC(Q1, sou, boundaries1, 1)
bc2 = DirichletBC(Q2, sou, boundaries2, 1)
bc3 = DirichletBC(Q3, sou, boundaries3, 1)

T1 = TrialFunction(Q1)
v1 = TestFunction(Q1)
T2 = TrialFunction(Q2)
v2 = TestFunction(Q2)
T3 = TrialFunction(Q3)
v3 = TestFunction(Q3)

dx1 = Measure('dx', domain=mesh1, subdomain_data=subdomains1)
ds1 = Measure('ds', domain=mesh1, subdomain_data=boundaries1)
dx2 = Measure('dx', domain=mesh2, subdomain_data=subdomains2)
ds2 = Measure('ds', domain=mesh2, subdomain_data=boundaries2)
dx3 = Measure('dx', domain=mesh3, subdomain_data=subdomains3)
ds3 = Measure('ds', domain=mesh3, subdomain_data=boundaries3)

a1 = (cro1 + rhoL*dphi_1) * T1/tau*v1*dx1 + inner(k1 * grad(T1), grad(v1))*dx1
L1 = (cro1 + rhoL*dphi_1) * T0_1/tau*v1*dx1

a2 = (cro2 + rhoL*dphi_2) * T2/tau*v2*dx2 + inner(k2 * grad(T2), grad(v2))*dx2
L2 = (cro2 + rhoL*dphi_2) * T0_2/tau*v2*dx2

a3 = (cro3 + rhoL*dphi_3) * T3/tau*v3*dx3 + inner(k3 * grad(T3), grad(v3))*dx3
L3 = (cro3 + rhoL*dphi_3) * T0_3/tau*v3*dx3

file1 = File("/mnt/c/Users/Vasya/Downloads/T1.pvd")
file2 = File("/mnt/c/Users/Vasya/Downloads/T2.pvd")
file3 = File("/mnt/c/Users/Vasya/Downloads/T3.pvd")

error_file1 = open("/mnt/c/Users/Vasya/Downloads/error_data1.dat", "w")
error_file2 = open("/mnt/c/Users/Vasya/Downloads/error_data2.dat", "w")

start_time = time.time()
t = tau
while t < Month:
  print(t)
  t += tau
  solve(a1 == L1, w1, bc1, solver_parameters=dict(linear_solver="cg", preconditioner="hypre_amg"))
  solve(a2 == L2, w2, bc2, solver_parameters=dict(linear_solver="cg", preconditioner="hypre_amg"))
  solve(a3 == L3, w3, bc3, solver_parameters=dict(linear_solver="cg", preconditioner="hypre_amg"))
  solve(a3 == L3, w3, bc3)
  
  w1_on_Q3 = interpolate(w1, Q3)
  E2_t_1 = inner(w3 - w1_on_Q3, w3 - w1_on_Q3)*dx3
  E2_b_1 = inner(w3, w3)*dx3
  E2_1 = sqrt(abs(assemble(E2_t_1))/abs(assemble(E2_b_1))) * 100
  
  w2_on_Q3 = interpolate(w2, Q3)
  E2_t_1 = inner(w3 - w2_on_Q3, w3 - w2_on_Q3)*dx3
  E2_b_1 = inner(w3, w3)*dx3
  E2_2 = sqrt(abs(assemble(E2_t_1))/abs(assemble(E2_b_1))) * 100
  
  error_file1.write(f"{t} {E2_1}\n")
  error_file2.write(f"{t} {E2_2}\n")
  
  T0_1.assign(w1)
  T0_2.assign(w2)
  T0_3.assign(w3)
  
  file1 << w1
  file2 << w2
  file3 << w3
  
solve_time = time.time() - start_time
print(solve_time, "seconds")

error_file1.close()
error_file2.close()