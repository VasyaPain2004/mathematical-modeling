from dolfin import *

# Только на root-процессе лог
parameters["std_out_all_processes"] = False

# Загружаем сетку (например, вокруг тела или в канале)
mesh = Mesh("mesh3.xml")

# Пространства (скорость — P2, давление — P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Функции и тестовые функции
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# ---------------------------
# Физические параметры
rho = 1.2                # кг/м^3
mu  = 1.8e-5             # Па·с
nu  = mu / rho           # кинематическая вязкость (≈1.5e-5 м^2/с)
U_inf = 10.0             # м/с
# ---------------------------

# Параметры по времени
dt = 0.001
T = 1.0
k = Constant(dt)

# Граничные условия
# 1. На внешней границе (x[0] > ... или особая метка) задаём поток
inflow = DirichletBC(V, Constant((U_inf, 0.0)), "on_boundary && x[0] < DOLFIN_EPS")

# 2. На теле/дне — условие прилипания (no-slip)
noslip = DirichletBC(V, Constant((0.0, 0.0)),
                     "on_boundary && (x[0] > 0.999 - DOLFIN_EPS || x[1] < DOLFIN_EPS)")

# 3. На выходе — условие давления p = 0
outflow = DirichletBC(Q, Constant(0.0), "on_boundary && x[0] > 0.999 - DOLFIN_EPS")

bcu = [noslip, inflow]
bcp = [outflow]

# Начальные функции
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Правая часть
f = Constant((0, 0))

# --- Вариационные формы ---

# Предварительный шаг скорости
F1 = (rho/k)*inner(u - u0, v)*dx \
     + rho*inner(dot(grad(u0), u0), v)*dx \
     + mu*inner(grad(u), grad(v))*dx \
     - inner(f, v)*dx

a1 = lhs(F1)
L1 = rhs(F1)

# Уравнение для давления
a2 = inner(grad(p), grad(q))*dx
L2 = -(rho/k)*div(u1)*q*dx

# Коррекция скорости
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - (k/rho)*inner(grad(p1), v)*dx

# Сборка матриц
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Решатели
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Файлы для записи
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")

# --- Временной цикл ---
t = dt
while t < T + DOLFIN_EPS:
    # Шаг 1: предварительная скорость
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "bicgstab", "default")

    # Шаг 2: давление
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    [bc.apply(p1.vector()) for bc in bcp]
    solve(A2, p1.vector(), b2, "bicgstab", prec)

    # Шаг 3: коррекция скорости
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "bicgstab", "default")

    # Сохраняем
    ufile << u1
    pfile << p1

    # Следующий шаг
    u0.assign(u1)
    t += dt
