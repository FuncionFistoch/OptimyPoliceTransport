# 1. El Metropolitano Monopolista

El caso más básico posible para plantear la modelización del flujo en las vías de una ciudad es suponer una larga carretera en línea recta, en la que un único transporte da el servicio en toda la ruta. Este caso encuentra ejemplos en las vías de tren y en el Metropolitano de Lima; ambos son similares, dado que se otorga a un solo operador el derecho de realizar la actividad económica.

Por lo que se plantea el problema desde los dos lados del mercado: primero la demanda del transporte, es decir, cuántas personas suben y bajan en cada paradero, y luego la oferta, es decir, cómo decide el conductor transportarlas y cuánto le cuesta hacerlo.

---

## 1.1 Demanda del transporte

Para todo el recorrido se supone una ruta con paraderos numerados $0, 1, \dots, K$. Se denomina $X_k$ a la distancia entre el paradero $k-1$ y el paradero $k$, de modo que la distancia total del recorrido es $X = \sum_{k=1}^{K} X_k$.

```
  |←—— X_1 ——→|←—— X_2 ——→|        ···        |←—— X_K ——→|
  ●———————————●———————————●—— ··· ——●—————————————————————●
Paradero 0  Paradero 1  Paradero 2          Paradero K-1  Paradero K
   ↑ u_0      ↑u_1 ↓d_1   ↑u_2 ↓d_2         ↑u_{K-1} ↓d_{K-1}  ↓d_K
```

En cada paradero $k$ hay personas que suben, $u_k$, y personas que bajan, $d_k$.

### Personas que suben

Las personas que suben en el paradero $k$ son todas las que llegaron a esperar desde el paso anterior del bus por ese paradero hasta el instante $t_k$ en que el bus llega:

$$
u_k \sim \text{Poisson}\left(\int_{t_k^{\text{ant}}}^{t_k} \lambda_k(s)\,ds\right)
$$

Es un proceso que depende de la intensidad de llegada $\lambda_k(s)$, y esta intensidad a su vez depende de la hora: en la entrada y en la salida laboral llega una enorme cantidad de personas en comparación con otras horas del día. Por lo que de esta manera se captan las diferencias temporales de la demanda de transporte (las "horas punta").

Se elige la distribución de Poisson porque cumple, en cierta medida, con sus supuestos:

1. La llegada de los pasajeros es de uno en uno.
2. La llegada de un pasajero es independiente de la llegada de otro.
3. La probabilidad de que llegue un pasajero es proporcional al intervalo de tiempo.

Además, es lo que se sigue de la literatura que modela el comportamiento de los pasajeros.

### Personas que bajan

Las personas que bajan no pueden modelarse igual, dado que dependen de la cantidad de personas que van en el bus: no pueden bajar más personas de las que hay a bordo. Si $A_k$ es la cantidad de pasajeros a bordo al llegar al paradero $k$ y cada uno baja ahí con probabilidad $q_k$:

$$
d_k \mid A_k \sim \text{Binomial}(A_k,\ q_k)
$$

### Aforo

Con este proceso se cumple además con el aforo máximo $\bar A$ del bus. Primero bajan las $d_k$ personas y luego suben solo las que caben:

$$
\bar u_k = \min\left\lbrace u_k,\ \bar A - (A_k - d_k)\right\rbrace 
$$

Este es el límite del aforo, y con él la cantidad de pasajeros evoluciona como:

$$
A_{k+1} = A_k - d_k + \bar u_k
$$

---

## 1.2 Oferta monopolista

El costo de dar el servicio se puede desglosar en tres partes:

- **El salario del conductor:** puede ser fijo o una fracción de lo que se recaude en el recorrido.
- **El consumo de combustible.**
- **El deterioro del vehículo.**

El conductor elige su aceleración $a(t)$, y con ello determina su posición $x(t)$ y su velocidad $v(t)$:

$$
\dot x(t) = v(t), \qquad \dot v(t) = a(t)
$$

Siguiendo la literatura sobre el CMEM (*Comprehensive Modal Emission Model*), el consumo de combustible es:

$$
f(v, a) = A\,v + B\,v^2 + C\,v^3 + D\,a\,v
$$

Esta función viene del diagrama de cuerpo libre del automóvil. La potencia que entrega el motor es:

$$
P = \left(m\,a + F_{\text{rod}} + \tfrac{1}{2}\rho\,C_d\,A_f\,v^2\right) v
$$

de donde salen la resistencia a la rodadura ($A\,v$), la resistencia aerodinámica ($C\,v^3$) y el cambio de energía cinética ($D\,a\,v$). El término $B\,v^2$ no sale de la dinámica: captura efectos que la dinámica no recoge, y puede ser nulo. Cabe mencionar que la función combina diferentes insumos, por lo que sus coeficientes se interpretan como precios relativos respecto del salario y de la depreciación del vehículo.

Se agrega además el desgaste por cambios bruscos de aceleración, que no está contemplado en el CMEM:

$$
g(a) = E\,a^2
$$

Su fundamento es el cambio brusco de energía cinética en un intervalo pequeño de tiempo.

### Costo total

Con $W$ la retribución al conductor y $T$ el tiempo que toma el recorrido:

$$
CT = W + \int_0^T \left[f(v, a) + g(a)\right] dt
$$

---

## 1.3 Problema del operador monopolista

Con una tarifa única $p$ (que es justamente el caso del Metropolitano y de muchas líneas de tren), el ingreso esperado del recorrido es:

$$
IT = \mathbb{E}\left[\sum_{k=1}^{n-1} p\,\bar u_k\right]
$$

Por lo que el problema del operador monopolista es maximizar su beneficio:

$$
\max_{a(t)}\ \pi = \mathbb{E}\left[\sum_{k=1}^{n-1} p\,\bar u_k\right] - W - \int_0^T \left[f(v, a) + g(a)\right] dt
$$

### El caso más simple: minimizar costos

Antes de resolver el problema estocástico, se supone que el conductor desconoce por completo la distribución de los pasajeros, por lo que prefiere minimizar el costo de dar el servicio. Suponiendo un salario que depende del tiempo, $W = \delta\,T$, el problema se separa en los tramos entre paraderos:

$$
\min_{a(t)}\ J = \sum_{k=0}^{n-1} \int_{t_k}^{t_{k+1}} \left[\delta + f(v, a) + g(a)\right] dt
$$

Como el conductor se detiene en cada paradero, cada tramo es un problema independiente con las mismas condiciones de partida y llegada. Por lo que el problema de optimización se convierte, para un solo tramo de longitud $X$, en:

$$
\min_{a(t)}\ \int_0^T \left[\delta + f(v, a) + g(a)\right] dt
$$

$$
\text{s.a.} \quad
\begin{cases}
\dot x = v, & x(0) = 0, \quad x(T) = X \\
\dot v = a, & v(0) = v(T) = 0 \\
a \in [\underline a,\ \bar a]
\end{cases}
$$

donde el tiempo final $T$ es libre: lo decide el propio conductor según qué tan rápido le convenga recorrer el tramo. Además, toda vía tiene un límite de velocidad, por lo que $v(t) \in [0,\ \bar v]$.

---

## 1.4 Hamiltoniano y los tres modelos

El hamiltoniano del problema es:

$$
H = \delta + f(v, a) + g(a) + \lambda_x\,v + \lambda_v\,a
$$

Según cómo se especifique el término de aceleración se tienen tres modelos.

### Modelo A (CMEM simétrico)

$$
H = \delta + A\,v + B\,v^2 + C\,v^3 + D\,a\,v + \lambda_x\,v + \lambda_v\,a
$$

Es la especificación original del CMEM. Tiene un problema: como $a\,v = \frac{d}{dt}\left(\frac{v^2}{2}\right)$ y el bus parte y llega detenido,

$$
\int_0^T D\,a\,v\,dt = D\left[\frac{v^2}{2}\right]_0^T = 0
$$

Es decir, todo lo que se paga al acelerar se "recupera" al frenar, y acelerar no tiene costo neto en el tramo.

### Modelo B (asimétrico)

$$
H = \delta + A\,v + B\,v^2 + C\,v^3 + D\,v\max[a, 0] + \lambda_x\,v + \lambda_v\,a
$$

En un vehículo de combustión, la energía para acelerar sale del combustible, pero al frenar esa energía no vuelve al tanque: los frenos la disipan como calor. Por lo que solo se paga cuando se acelera, y desacelerar no representa ningún ahorro. (En un vehículo eléctrico con freno regenerativo sí se recuperaría una parte, y el caso quedaría entre el modelo A y el B.)

### Modelo C (asimétrico + desgaste)

$$
H = \delta + A\,v + B\,v^2 + C\,v^3 + D\,v\max[a, 0] + E\,a^2 + \lambda_x\,v + \lambda_v\,a
$$

Al modelo B se le agrega el desgaste por cambios bruscos de aceleración. Este término sí castiga frenar de golpe, y hace que el conductor prefiera acelerar y frenar de forma suave, que es lo que se observa en la realidad.

La solución de estos tres modelos se obtiene numéricamente; el método se explica en el archivo `2OfertaMonopolista.md`.
