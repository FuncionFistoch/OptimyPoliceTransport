import sympy as sy
#Defined constants
A = 0.1 #0.1 para todo, a = 0.1, X <= 3.7 
B = 0.1
C = 0.1
D = 0.1
a_u = 0.1
a_l = -a_u
X = 20000
delta = 5

#Defined parameters of the solution optimal
c1, T1, lambda_x, v, c2, T2, T = sy.symbols('c1 T1 lambda_x v c2 T2 T', real= True)

#Defined lambda_v
c, t, a = sy.symbols('c t a')
lambda_v = c - (A + lambda_x)*t - B*a*t**2 - C*a**2*t**3
lambda_vT1 = lambda_v.subs({c: c1, t: T1, a:a_u})
lambda_vT2 = lambda_v.subs({c: c2, t: T2, a:a_l})
lambda_v0 = lambda_v.subs({c: c1, t:0, a: a_u})
lambda_vT = lambda_v.subs({c: c2, t: T, a: a_l})
#Solap into coestate_v and velocity of cruise in acelerate(T1)
expr1 = lambda_vT1 + D*v
#Solap into coestate_v and velocity of cruise in desacelerate(T2)
expr2 = lambda_vT2 #+ D*v
#Total distance in T time
expr3 = 0.5*a_u*T1**2 + (T2-T1)*v - 0.5*a_l*(T-T2)**2 - X
#Hamiltoniano autonomus
expr4 = delta + lambda_v0*a_u
expr5 = delta + lambda_vT*a_l
#Acelerate constant since 0 to T1
expr6 = v - a_u*T1
#Acelerate constante since T2 to T
expr7 = v + a_l*(T - T2)

system = [expr1, expr2, expr3, expr4, expr5, expr6, expr7]
solutions = sy.nonlinsolve(system, [c1, T1, lambda_x, v, c2, T2, T])

#Filter the Real Solutions
solutions_reals = [s for s in solutions if (all(sy.im(sy.N(v)) == 0 for v in s)) ]
print('c1, T1, lambda_x, v, c2, T2, T')
print(solutions_reals)

def cost_phase(v0, accel, duration):
    """
    Integrates L(v,a) over a phase where v(t) = v0 + accel*t,
    for t in [0, duration], with accel CONSTANT (use accel=0 for cruise).
    """
    tt = sy.symbols('tt', nonnegative=True)
    v_t = v0 + accel*tt
    if accel > 0:
        L = delta + A*v_t + B*v_t**2 + C*v_t**3 + D*v_t*accel
    else:
        L = delta + A*v_t + B*v_t**2 + C*v_t**3
    return sy.N(sy.integrate(L, (tt, 0, duration)))
 
def total_cost_bang_singular_bang(a_u, a_l, T1_val, T2_val, T_val):
    """
    Total cost for the 3-phase structure (bang-singular-bang):
    Phase 1: accelerate 0 -> T1 with a_u, starting from v=0
    Phase 2: cruise T1 -> T2 at constant v = a_u*T1 (a=0)
    Phase 3: decelerate T2 -> T with a_l, ending at v=0
    """
    v_cruise = a_u * T1_val
    c1 = cost_phase(v0=0,        accel=a_u, duration=T1_val)
    c2 = cost_phase(v0=v_cruise, accel=0,   duration=T2_val - T1_val)
    c3 = cost_phase(v0=v_cruise, accel=a_l, duration=T_val - T2_val)
    return c1 + c2 + c3
 
# Apply to each real, physically valid solution found above
for s in solutions_reals:
    c1_v, T1_v, lambda_x_v, v_v, c2_v, T2_v, T_v = [sy.N(x) for x in s]
    costo = total_cost_bang_singular_bang(a_u, a_l, T1_v, T2_v, T_v)
    print(f"T1={T1_v}, T2={T2_v}, T={T_v}  ->  Costo total = {costo}")