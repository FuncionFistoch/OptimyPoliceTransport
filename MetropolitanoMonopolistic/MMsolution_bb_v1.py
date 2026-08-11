import sympy as sy
#Defined constants
A = 0.1 #0.1 para todo, a = 0.1, X <= 3.7 
B = 0.1
C = 0.1
D = 0.1
a_u = 0.1
a_l = -a_u
X = 0.8
delta = 5

#Defined parameters of the solution optimal
c1, c2, T1, lambda_x, T = sy.symbols('c1 c2 T1 lambda_x T')

#Defined lambda_v for solap 
c, t, a = sy.symbols('c t a')
lambda_v = c - (A + lambda_x)*t - B*a*t**2 - C*a**2*t**3
lambda_v_acelerate = lambda_v.subs({c: c1, t: T1, a: a_u})
lambda_v_deacelerate = lambda_v.subs({c: c2, t: T1, a: a_l})
lambda_v0 = lambda_v.subs({c: c1, t: 0, a: a_u})
lambda_vT = lambda_v.subs({c: c2, t: T, a: a_l})

#Lambda_v Slop
expr1 = lambda_v_acelerate - lambda_v_deacelerate
#Total Length
expr2 = 0.5*a_u*T1*T - X
#Hamiltoniano autonomus
expr3 = delta + lambda_v0*a_u
expr4 = delta + lambda_vT*a_l
#Velocity Slop
expr5 = a_u*T1 + a_l*(T-T1)

system = [expr1, expr2, expr3, expr4, expr5]
solutions = sy.nonlinsolve(system, [c1, c2, T1, lambda_x, T])
#Filter solutions real
solutions_reals = [s for s in solutions if (all(sy.im(sy.N(v)) == 0 for v in s) and 
                                            (s[-1] > 0) and
                                            (s[-3] > 0) and 
                                            (s[-1] > s[-3]))]
print('c1 c2 T1 lambda_x T')
print(solutions_reals)

def cost_phase(v0, accel, duration):
    """
    Integrates L(v,a) over a phase where v(t) = v0 + accel*t,
    for t in [0, duration], with accel CONSTANT.
    """
    tt = sy.symbols('tt', nonnegative=True)
    v_t = v0 + accel*tt
    if accel > 0:
        L = delta + A*v_t + B*v_t**2 + C*v_t**3 + D*v_t*accel
    else:
        L = delta + A*v_t + B*v_t**2 + C*v_t**3
    return sy.N(sy.integrate(L, (tt, 0, duration)))
 
def total_cost_bang_bang(a_u, a_l, T1_val, T_val):
    """
    Total cost for the 2-phase structure (pure bang-bang, no cruise):
    Phase 1: accelerate 0 -> T1 with a_u, starting from v=0
    Phase 2: decelerate T1 -> T with a_l, ending at v=0
    """
    v1 = a_u * T1_val  # peak velocity at the end of the acceleration phase
    c1 = cost_phase(v0=0,  accel=a_u, duration=T1_val)
    c2 = cost_phase(v0=v1, accel=a_l, duration=T_val - T1_val)
    return c1 + c2
 
# Apply to each real, physically valid solution found above
for s in solutions_reals:
    c1_v, c2_v, T1_v, lambda_x_v, T_v = [sy.N(x) for x in s]
    costo = total_cost_bang_bang(a_u, a_l, T1_v, T_v)
    print(f"T1={T1_v}, T={T_v}  ->  Total Cost = {costo}")