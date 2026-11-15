"""
MM_conductor_optimo.py  (v4 -- solo solución numérica)
=======================================================
El Metropolitano Monopolístico: tramo entre dos paraderos.

  min J = ∫_0^T L(v,a) dt,   ẋ = v,  v̇ = a
  x(0)=0, v(0)=0,  x(T)=X, v(T)=0,  T libre
  a ∈ [a_l, a_u],  v ∈ [0, v̄]

Costos instantáneos:
  A  L = δ + A v + B v² + C v³ + D·a·v
  B  L = δ + A v + B v² + C v³ + D·max(a,0)·v
  C  L = δ + A v + B v² + C v³ + D·max(a,0)·v + E·a²

Método: control constante por tramos + simulación exacta de la dinámica
(simular) + SLSQP. Ver docstring de simular().

Salida: resultados_MM.xlsx con las hojas
  resumen, perfiles, superficie_C, elast_delta_C, elast_X_C, parametros
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize

# ============================================================================
# 1. PARÁMETROS (metros, segundos; coeficientes = precios relativos)
# ============================================================================
P = dict(
    A=0.5,        # resistencia a la rodadura (∝ v)
    B=0.01,       # efectos no captados por la dinámica (∝ v²)
    C=0.001,      # resistencia aerodinámica (∝ v³)
    D=1.0,        # costo de ganar energía cinética (∝ a·v)
    E=10.0,       # desgaste por aceleraciones bruscas (∝ a²), solo C
    delta=5.0,    # salario por unidad de tiempo: W = δ·T
    a_u=1.0,      # aceleración máxima    [m/s²]
    a_l=-1.2,     # desaceleración máxima [m/s²]
    v_bar=16.7,   # límite de velocidad   [m/s]
    X=800.0,      # distancia entre paraderos [m]
)

MODELOS = {
    "A": dict(nombre="A: CMEM simétrico", asim=False, desgaste=False),
    "B": dict(nombre="B: asimétrico", asim=True, desgaste=False),
    "C": dict(nombre="C: asimétrico + desgaste", asim=True, desgaste=True),
}

N_NODOS = 160  # número de tramos de control


# ============================================================================
# 2. COSTO INSTANTÁNEO SIN TÉRMINO DE ACELERACIÓN
# ============================================================================
def L0(v, p):
    """δ + A·v + B·v² + C·v³  (cúbico en v)."""
    return p["delta"] + p["A"] * v + p["B"] * v**2 + p["C"] * v**3


# ============================================================================
# 3. SOLUCIÓN NUMÉRICA
# ============================================================================
def simular(z, N):
    """
    Dado el vector de decisión z = [T, a⁺_0..a⁺_{N-1}, a⁻_0..a⁻_{N-1}],
    construye la trayectoria discreta que cumple la dinámica exacta.

    Discretización (N tramos, h = T/N, t_k = k·h):
      a_k = a⁺_k − a⁻_k ,  a⁺_k, a⁻_k ≥ 0      (control constante en [t_k, t_{k+1}])

    Con aceleración constante, la solución exacta de la dinámica en un tramo es
      v(t_{k+1}) = v(t_k) + a_k·h
      x(t_{k+1}) = x(t_k) + v(t_k)·h + a_k·h²/2

    Por inducción, con v_0 = 0 y x_0 = 0:
      v_k = h · Σ_{j<k} a_j
      x_k = h · Σ_{j<k} v_j + (h²/2) · Σ_{j<k} a_j

    Esto se realiza con cumsum. La dinámica se cumple exactamente
    dado el control; la única aproximación es que el control (por ello su forma escalonada) es constante por tramos.
    Además las condiciones de frontera x(0)=0, v(0)=0 quedan satisfechas por
    construcción, y las restricciones v(T)=0 y x(T)=X se imponen como restricciones
    del problema de optimización.

    Devuelve T, h, a⁺, a⁻, a, v (N+1 nodos), x (N+1 nodos).
    """
    T = z[0]
    ap, am = z[1:N + 1], z[N + 1:]
    a = ap - am
    h = T / N
    v = np.concatenate([[0.0], np.cumsum(a * h)])
    x = np.concatenate([[0.0], np.cumsum(v[:-1] * h + a * h**2 / 2)])
    return T, h, ap, am, a, v, x


def costo_numerico(z, N, p, modelo):
    """
    J = ∫_0^T L dt, evaluada exactamente para el control por tramos:
      • Término δ + A·v + B·v² + C·v³: Simpson por tramo. v es lineal en cada
        tramo, así que L0 es cúbico en t, y Simpson integra exacto polinomios
        de grado ≤ 3.
      • Término D·φ(a)·v: ∫ v dt por la regla del trapecio (exacta para v lineal).
        φ(a) = a en A, y φ(a) = a⁺ en B y C.
      • Término E·a²: a es constante por tramo, así que la integral es exacta.
    """
    T, h, ap, am, a, v, x = simular(z, N)
    v_mid = v[:-1] + a * h / 2
    J = np.sum(h / 6 * (L0(v[:-1], p) + 4 * L0(v_mid, p) + L0(v[1:], p)))
    int_v = h * (v[:-1] + v[1:]) / 2
    J += p["D"] * np.sum((ap if modelo["asim"] else a) * int_v)
    if modelo["desgaste"]:
        J += p["E"] * np.sum(a**2) * h
    """
    Se agrega para solucionar problemas con el modelo A 
    """
    J += 1e-6 * np.sum(ap * am)
    return J


def k_frenado(p):
    """Distancia que requiere acelerar hasta v y frenar: k·v²."""
    return 1 / (2 * p["a_u"]) + 1 / (2 * abs(p["a_l"]))


def semilla(p, N):
    """
    Punto inicial cinemático: bang-crucero-bang que cumple v(0)=v(T)=0 y
    x(T)=X, sin usar el costo. El crucero va a v_pk = min(v̄, √(X/k)).
    Es solo un punto de partida factible para el algoritmo SLSQP.
    """
    v_pk = min(p["v_bar"], np.sqrt(p["X"] / k_frenado(p)))
    t1 = v_pk / p["a_u"]
    t3 = v_pk / abs(p["a_l"])
    tau = (p["X"] - k_frenado(p) * v_pk**2) / v_pk
    T0 = t1 + tau + t3
    t_mid = (np.arange(N) + 0.5) * T0 / N
    a = np.where(t_mid < t1, p["a_u"], np.where(t_mid < t1 + tau, 0.0, p["a_l"]))
    return np.concatenate([[T0], np.clip(a, 0, None), np.clip(-a, 0, None)])


def solucion_numerica(p, modelo, z0, N=N_NODOS):
    """
    Problema de programación no lineal en z (2N+1 variables):

      min_z   J(z)                     (costo_numerico)
      s.a.    v_N(z) = 0               (restricción de frontera)
              x_N(z) = X               (restricción de frontera)
              v_k(z) ≥ 0,  k=1..N-1    (restricción de estado)
              v̄ − v_k(z) ≥ 0           (restricción de estado)
              0 ≤ a⁺_k ≤ a_u,  0 ≤ a⁻_k ≤ |a_l|,  T ∈ [1, 10·T0]

    Se resuelve con SLSQP. Las restricciones de frontera y de estado se
    evalúan sobre las trayectorias de simular(z).
    """
    escala = 1e-3
    restricciones = [
        dict(type="eq", fun=lambda z: simular(z, N)[5][-1]),
        dict(type="eq", fun=lambda z: simular(z, N)[6][-1] - p["X"]),
        dict(type="ineq", fun=lambda z: simular(z, N)[5][1:-1]),
        dict(type="ineq", fun=lambda z: p["v_bar"] - simular(z, N)[5]),
    ]
    cotas = [(1.0, 10 * z0[0])] + [(0, p["a_u"])] * N + [(0, abs(p["a_l"]))] * N
    res = minimize(lambda z: escala * costo_numerico(z, N, p, modelo), z0,
                   method="SLSQP", bounds=cotas, constraints=restricciones,
                   options=dict(maxiter=1000, ftol=1e-12))
    T, h, ap, am, a, v, x = simular(res.x, N)
    return dict(res=res, T=T, costo=res.fun / escala, a=a, v=v, x=x,
                t=np.linspace(0, T, N + 1), N=N)


def resolver(p, clave, N=N_NODOS):
    """Resuelve el modelo `clave` desde la semilla cinemática."""
    return solucion_numerica(p, MODELOS[clave], semilla(p, N), N)


# ============================================================================
# 4. ELASTICIDADES
# ============================================================================
def elasticidades(x, v, T):
    """ε = d ln y / d ln x por diferencias finitas sobre una malla no uniforme."""
    lx = np.log(x)
    return np.gradient(np.log(v), lx), np.gradient(np.log(T), lx)


# ============================================================================
# 5. EJECUCIÓN Y EXPORTACIÓN A EXCEL
# ============================================================================
def main():
    # ---- Perfiles y resumen de los tres modelos ----
    perfiles, resumen = [], []
    for clave in ("A", "B", "C"):
        sol = resolver(P, clave)
        a = np.append(sol["a"], np.nan)  # el control vive en tramos: último nodo sin valor
        perfiles.append(pd.DataFrame(dict(modelo=clave, t=sol["t"], v=sol["v"],
                                          x=sol["x"], a=a)))
        resumen.append(dict(modelo=clave, nombre=MODELOS[clave]["nombre"],
                            costo=sol["costo"], T=sol["T"], v_max=sol["v"].max(),
                            exito=bool(sol["res"].success),
                            iteraciones=int(sol["res"].nit),
                            mensaje=sol["res"].message))
    df_perfiles = pd.concat(perfiles, ignore_index=True)
    df_resumen = pd.DataFrame(resumen)

    # ---- Superficie del modelo C: velocidad máxima en (δ, X) ----
    DELTAS_SURF = [1, 3, 5, 8, 12]
    XS_SURF = [200, 400, 600, 800, 1000, 1200, 1600]
    filas = []
    for d in DELTAS_SURF:
        for X in XS_SURF:
            sol = resolver({**P, "delta": d, "X": X}, "C")
            filas.append(dict(delta=d, X=X, v_max=sol["v"].max(), T=sol["T"],
                              exito=bool(sol["res"].success)))
    df_superficie = pd.DataFrame(filas)

    # ---- Elasticidades del modelo C ----
    DELTAS_SW = [1, 2, 3, 4, 5, 6, 8, 10, 12]
    XS_SW = [200, 300, 400, 600, 800, 1000, 1200, 1600]

    d_arr = np.array(DELTAS_SW, dtype=float)
    v_d, T_d = np.empty(len(d_arr)), np.empty(len(d_arr))
    for i, d in enumerate(d_arr):
        sol = resolver({**P, "delta": d}, "C")
        v_d[i], T_d[i] = sol["v"].max(), sol["T"]
    ev_d, eT_d = elasticidades(d_arr, v_d, T_d)
    df_el_d = pd.DataFrame(dict(delta=d_arr, v_max=v_d, T=T_d,
                                eps_v_delta=ev_d, eps_T_delta=eT_d))

    x_arr = np.array(XS_SW, dtype=float)
    v_x, T_x = np.empty(len(x_arr)), np.empty(len(x_arr))
    for i, X in enumerate(x_arr):
        sol = resolver({**P, "X": X}, "C")
        v_x[i], T_x[i] = sol["v"].max(), sol["T"]
    ev_x, eT_x = elasticidades(x_arr, v_x, T_x)
    df_el_x = pd.DataFrame(dict(X=x_arr, v_max=v_x, T=T_x,
                                eps_v_X=ev_x, eps_T_X=eT_x))

    # ---- Parámetros ----
    df_par = pd.DataFrame([dict(parametro=k, valor=v) for k, v in P.items()] +
                          [dict(parametro="N_nodos", valor=N_NODOS)])

    # ---- Exportar ----
    salida = "resultados_MM.xlsx"
    with pd.ExcelWriter(salida, engine="openpyxl") as w:
        df_resumen.to_excel(w, sheet_name="resumen", index=False)
        df_perfiles.to_excel(w, sheet_name="perfiles", index=False)
        df_superficie.to_excel(w, sheet_name="superficie_C", index=False)
        df_el_d.to_excel(w, sheet_name="elast_delta_C", index=False)
        df_el_x.to_excel(w, sheet_name="elast_X_C", index=False)
        df_par.to_excel(w, sheet_name="parametros", index=False)

    print(df_resumen.to_string(index=False))
    print(f"\nSuperficie C: {int(df_superficie['exito'].sum())}/{len(df_superficie)} convergieron")
    print(f"\nArchivo generado: {salida}")


if __name__ == "__main__":
    main()
