from dolfin import *

meshname = "mesh3"
mesh = Mesh(meshname + ".xml")
subdomains = MeshFunction("size_t", mesh, meshname + "_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, meshname + "_facet_region.xml")

Q = FunctionSpace(mesh, "CG", 1)
DAY = 86400
Month = 30*DAY
tau = 0.5*DAY

sou = Constant(-30.0)
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

bc1 = DirichletBC(Q, sou, boundaries, 1)

T = TrialFunction(Q)
v = TestFunction(Q)

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

a = (cro + rhoL*dphi) * T/tau*v*dx + inner(k * grad(T), grad(v))*dx
L = (cro + rhoL*dphi) * T0/tau*v*dx

fileT = File("/mnt/c/Users/Vasya/Downloads/T1.pvd")

t = tau
while t < Month:
  print(t)
  t += tau
  solve(a == L, w, bc1)
  T0.assign(w)
  fileT << w