from dolfin import *

meshname = "mesh3"
mesh = Mesh(meshname + ".xml")
subdomains = MeshFunction("size_t", mesh, meshname + "_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, meshname + "_facet_region.xml")
      
k = Constant(200.0)
mu = Constant(0.8e9)
lmbda = Constant(1.25e9)
beta = 1.0e-5*(3*lmbda + 2*mu)
C = 1.0e6
alpAir = 100.0
TAir = 1.0
THot = 1000.0
tmax = 10000
dt = tmax/10

def epsilon(u):
  return 0.5*(grad(u) + grad(u).T)

def sigma(u):
  return lmbda*div(u)*Identity(2) + 2*mu*epsilon(u)

V = VectorElement("CG", mesh.ufl_cell(), 1)
Q = FiniteElement("CG", mesh.ufl_cell(), 1)
TH = V * Q
W = FunctionSpace(mesh, TH)

un_val = Constant((0.0, 0.0))
Tn_val = Constant(TAir)
un = interpolate(un_val, W.sub(0).collapse())
Tn = interpolate(Tn_val, W.sub(1).collapse())
w = Function(W)
g0 = Constant(0.0)
bc2 = DirichletBC(W.sub(0).sub(0), g0, boundaries, 2)
bc3 = DirichletBC(W.sub(0).sub(1), g0, boundaries, 3)
bcT1 = DirichletBC(W.sub(1), THot, boundaries, 1)
bcT3 = DirichletBC(W.sub(1), TAir, boundaries, 3)
bcT4 = DirichletBC(W.sub(1), THot, boundaries, 4)
bcs = [bc2, bc3, bcT1, bcT3, bcT4]
u, T = TrialFunctions(W)
v, q = TestFunctions(W)
dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
a = inner(sigma(u), epsilon(v))*dx + beta*inner(grad(T), v)*dx \
+ beta*Tn_val*div(u)*q*dx + C*T*q*dx \
+ dt*inner(k*grad(T), grad(q))*dx
L = C*Tn*q*dx + beta*Tn_val*div(un)*q*dx

w.rename('u', '0')
fileu = File("/mnt/c/Users/Vasya/Downloads/u.pvd")
fileT = File("/mnt/c/Users/Vasya/Downloads/T.pvd")
t = 0
while t < tmax:
  print(t)
  t += dt
  solve(a == L, w, bcs)
  (uu, TT) = w.split(deepcopy = True)
  un.assign(uu)
  Tn.assign(TT)
  print('w(%g, %g)' % (w.vector().min(), w.vector().max()))
  fileu << uu
  fileT << TT