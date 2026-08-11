import sympy as sy

A = 1
B = 1
C = 1
D = 5
a_u = 1
a_l = -1
delta = 30
X = 500
c1, c2, c3, T1, T2, T, lambda_x, v = sy.symbols("c1 c2 c3 T1 T2 T lambda_x v")
A_, B_, C_, t = sy.symbols("A_ B_ C_ t")
K = sy.symbols("K")   # generic integration constant (NOT "C" -> avoids shadowing C=0.1)

# --- Fase 1: acelerar 0->T1, L_accel = delta + A v + B v^2 + C v^3 + D v a ---
lambda_v1 = K - (A + lambda_x + D*a_u)*t - B*a_u*t**2 - C*a_u**2*t**3

# --- Fase 2: crucero T1->T2, a=0, L_accel = L_brake aqui (D*a=0 de todas formas) ---
v_c = a_u*T1
lambda_v2 = K - (A + lambda_x + 2*B*v_c + 3*C*v_c**2)*(t - T1)

# --- Fase 3: frenar T2->T, L_brake = delta + A v + B v^2 + C v^3 (SIN D v a) ---
v0 = a_u*T1 - a_l*T2   # velocidad "aparente" en t=0 si se extendiera la recta hacia atras
A_ = A + lambda_x + 2*B*v0 + 3*C*v0**2
B_ = a_l*(B + 3*C*v0)
C_ = C*a_l**2
lambda_v3 = K - A_*t - B_*t**2 - C_*t**3

# --- Continuidad de lambda_v en los switches ---
# FIX 1: expr1 debe comparar ambas fases en el MISMO instante T1 (antes comparaba T1 vs T2)
expr1 = lambda_v1.subs({K: c1, t: T1}) - lambda_v2.subs({K: c2, t: T1})
expr2 = lambda_v2.subs({K: c2, t: T2}) - lambda_v3.subs({K: c3, t: T2})

# --- Distancia total (pura cinematica, no depende de L) ---
expr3 = 0.5*(T + T2 - T1)*a_u*T1 - X

# --- Hamiltoniano autonomo = 0 en los extremos ---
expr4 = delta + lambda_v3.subs({K: c3, t: T})*a_l
expr7 = delta + lambda_v1.subs({K: c1, t: 0})*a_u

# --- Cinematica de velocidad ---
expr5 = v - a_u*T1
expr6 = v + a_l*(T - T2)

# --- Condicion de conmutacion (usa L_accel, con D) ---
expr8 = D*v + lambda_v1.subs({K: c1, t: T1})

system = [expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8]
solutions = sy.nonlinsolve(system, [c1, c2, c3, T1, T2, T, lambda_x, v])

solutions_reals = [s for s in solutions if all(sy.im(sy.N(x)) == 0 for x in s )
                                        and s[3] > 0
                                        and s[4] > 0
                                        and s[5] > 0
                                        and s[-1]> 0 
                                        and s[5] > s[4]
                                        and s[4] > s[3]]

print('c1, c2, c3, T1, T2, T, lambda_x, v')
print(solutions_reals)