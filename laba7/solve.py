import time
from dolfin import *

meshname = "mesh2"
mesh = Mesh(meshname + ".xml")
subdomains = MeshFunction("size_t", mesh, meshname + "_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, meshname + "_facet_region.xml")

Q = FunctionSpace(mesh, "CG", 1)
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
T0 = interpolate(initT, Q)

phi = conditional(gt(T0, delta), 1.0, conditional(lt(T0, -delta), 0.0,
(T0+delta)/(delta+delta)))
dphi = conditional(gt(T0, delta), 0.0, conditional(lt(T0, -delta),
0.0, 1.0/(delta+delta)))

k = kS + phi*(kL - kS) #lambda(phi)
cro = croS + phi*(croL - croS) #alpha(phi)

w = Function(Q)

bc = DirichletBC(Q, sou, boundaries, 1)

T = TrialFunction(Q)
v = TestFunction(Q)

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

a = (cro + rhoL*dphi) * T/tau*v*dx + inner(k * grad(T), grad(v))*dx
L = (cro + rhoL*dphi) * T0/tau*v*dx

file = File("/mnt/c/Users/Vasya/Downloads/T1.pvd")

start_time = time.time()
t = tau

solver = KrylovSolver('richardson', 'amg') 
cg_prm = solver.parameters
cg_prm["absolute_tolerance"] = 1E-7 
cg_prm["relative_tolerance"] = 1E-4 
cg_prm["maximum_iterations"] = 100000

while t < Month:
  t += tau
  res = w.vector()
  A = assemble(a)
  b = assemble(L)
  bc.apply(A, b)
  iter = solver.solve(A, res, b) 
  print(iter)  
  
  T0.assign(w)
  
  file << w
  
solve_time = time.time() - start_time
print(solve_time, "seconds")
