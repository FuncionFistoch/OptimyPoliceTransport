# 2. Solución numérica del problema del conductor

En el archivo anterior se dejó planteado el problema de minimización de costos para un tramo entre dos paraderos:

$$
\min_{a(t),\,T}\ \int_0^T L(v, a)\,dt
\qquad \text{s.a.} \quad
\begin{cases}
\dot x = v, & x(0) = 0, \quad x(T) = X \\
\dot v = a, & v(0) = v(T) = 0 \\
a \in [\underline a,\ \bar a], & v \in [0,\ \bar v]
\end{cases}
$$

donde $L(v, a)$ es el costo instantáneo de cada modelo:

| Modelo | $L(v, a)$ |
|---|---|
| A | $\delta + A\,v + B\,v^2 + C\,v^3 + D\,a\,v$ |
| B | $\delta + A\,v + B\,v^2 + C\,v^3 + D\,v\max[a, 0]$ |
| C | $\delta + A\,v + B\,v^2 + C\,v^3 + D\,v\max[a, 0] + E\,a^2$ |

Resolverlo con el principio del máximo de Pontryagin exige conocer de antemano la estructura de la solución (si el conductor acelera, va a velocidad crucero y frena, o si solo acelera y frena), y además tratar aparte la restricción $v \le \bar v$. En los primeros intentos con SymPy esto resultó bastante complicado, y para el modelo C ni siquiera hay una forma cerrada. Por lo que se opta por una estrategia más sencilla: **primero discretizar y luego optimizar**. Es decir, se convierte el problema de buscar una función $a(t)$ en un problema de buscar un vector de números, que se puede resolver con un optimizador numérico.

---

## 2.1 Normalizar el tiempo

El tiempo final $T$ es libre, por lo que no se puede fijar de antemano una malla de tiempo. Para resolver esto se normaliza el tiempo con $\tau = t/T \in [0, 1]$:

$$
\frac{dx}{d\tau} = T\,v, \qquad \frac{dv}{d\tau} = T\,a
$$

Ahora el horizonte siempre es $[0, 1]$ y $T$ pasa a ser una variable más que se optimiza. Al dividir el recorrido en $N$ tramos iguales, cada tramo dura en tiempo real:

$$
h = \frac{T}{N}, \qquad t_k = k\,h
$$

---

## 2.2 Discretizar el control

Se supone que en cada tramo el conductor mantiene una aceleración constante $a_k$. Además, para mantener suave la función $\max[a, 0]$, cada aceleración se separa en su parte positiva y su parte negativa:

$$
a_k = a_k^+ - a_k^-, \qquad a_k^+ \ge 0, \quad a_k^- \ge 0
$$

y se reemplaza $\max[a_k, 0]$ por $a_k^+$. El motivo es que el optimizador que se usa trabaja con derivadas, y $\max[a, 0]$ tiene una esquina justo en $a = 0$, que es donde el conductor va a velocidad crucero. Con la separación, la esquina desaparece del costo y se convierte en una restricción de signo, que es mucho más fácil de manejar. Como acelerar cuesta, el optimizador nunca va a elegir $a_k^+$ y $a_k^-$ positivos a la vez (eso sería pagar por acelerar para luego frenar sin cambiar nada), así que en la solución $a_k^+ = \max[a_k, 0]$.

Con ello, todo lo que hay que encontrar cabe en un solo vector:

$$
z = \big[\,T;\ \underbrace{a_0^+, \dots, a_{N-1}^+}_{\in\,\mathbb{R}^N};\ \underbrace{a_0^-, \dots, a_{N-1}^-}_{\in\,\mathbb{R}^N}\,\big] \in \mathbb{R}^{2N+1}
$$

---

## 2.3 Reconstruir la trayectoria

Con la aceleración constante en un tramo pequeño, la dinámica se resuelve de forma exacta:

$$
v_{k+1} = v_k + a_k\,h, \qquad x_{k+1} = x_k + v_k\,h + \tfrac{1}{2}\,a_k\,h^2
$$

Y como el bus parte detenido, $v_0 = 0$ y $x_0 = 0$, por inducción se tiene:

$$
v_k = h\sum_{j=0}^{k-1} a_j, \qquad x_k = h\sum_{j=0}^{k-1} v_j + \frac{h^2}{2}\sum_{j=0}^{k-1} a_j
$$

Esto es lo que hace la función `simular` del código con sumas acumuladas (`cumsum`). Cabe mencionar dos cosas:

- La única aproximación del método es que la aceleración es constante por tramos (por eso su forma escalonada en los resultados). Dada esa aceleración, la velocidad y la posición son exactas.
- Las condiciones de partida $v(0) = 0$ y $x(0) = 0$ se cumplen por construcción, así que no hace falta imponerlas. Solo quedan por imponer las de llegada: $v_N = 0$ y $x_N = X$.

---

## 2.4 El costo del recorrido

El costo $J(z)$ es la integral del costo instantáneo, resuelta tramo por tramo:

$$
J(z) = \int_0^T L(v, a)\,dt = \sum_{k=0}^{N-1} \int_{t_k}^{t_{k+1}} L\big(v(t), a_k\big)\,dt
$$

Dentro de cada tramo la velocidad es lineal, $v(t) = v_k + a_k\,(t - t_k)$, por lo que cada parte del costo se integra de forma exacta:

**Parte $L_0 = \delta + A\,v + B\,v^2 + C\,v^3$.** Como $v$ es lineal en $t$, $L_0$ es un polinomio cúbico en $t$, y la regla de Simpson integra sin error polinomios de grado menor o igual a 3:

$$
\int_{t_k}^{t_{k+1}} L_0\,dt = \frac{h}{6}\Big[L_0(v_k) + 4\,L_0(v_{k+\frac{1}{2}}) + L_0(v_{k+1})\Big], \qquad v_{k+\frac{1}{2}} = v_k + \frac{a_k\,h}{2}
$$

**Parte de aceleración.** La aceleración es constante en el tramo, y la integral de la velocidad (lineal) es exacta con el trapecio:

$$
\int_{t_k}^{t_{k+1}} D\,\varphi_k\,v\,dt = D\,\varphi_k\,\frac{h\,(v_k + v_{k+1})}{2},
\qquad
\varphi_k =
\begin{cases}
a_k & \text{Modelo A} \\
a_k^+ & \text{Modelos B y C}
\end{cases}
$$

**Parte de desgaste (solo modelo C).**

$$
\int_{t_k}^{t_{k+1}} E\,a_k^2\,dt = E\,a_k^2\,h
$$

Juntando todo:

$$
J(z) = \sum_{k=0}^{N-1}\left\lbrace \frac{h}{6}\Big[L_0(v_k) + 4\,L_0(v_{k+\frac{1}{2}}) + L_0(v_{k+1})\Big] + D\,\varphi_k\,\frac{h\,(v_k + v_{k+1})}{2} + E\,a_k^2\,h\right\rbrace  + 10^{-6}\sum_{k=0}^{N-1} a_k^+\,a_k^-
$$

donde $h = T/N$ y las velocidades $v_k$ dependen de $z$. Cabe mencionar que el término $\delta$ suma exactamente $\delta\,T$, que es el salario: así es como el tiempo entra en el costo y tardar más no sale gratis.

**Sobre el último término.** En el modelo A el costo depende solo de $a_k = a_k^+ - a_k^-$, y no de cada parte por separado. Por lo que hay infinitas combinaciones de $a_k^+$ y $a_k^-$ con el mismo costo, y el optimizador avanza muy lento. El término $10^{-6}\sum a_k^+\,a_k^-$ resuelve esa ambigüedad sin alterar la solución, porque su aporte es del orden de $10^{-4}$ frente a un costo del orden de $1000$. En los modelos B y C no hace falta, porque el término $D\,a_k^+\,v$ ya cumple ese papel. Una mejora pendiente es parametrizar el modelo A directamente con $a_k$, sin separarla, ya que ese modelo no tiene la función $\max$.

---

## 2.5 El problema de optimización

Con todo lo anterior, el problema del conductor queda como un problema de optimización en $2N+1$ variables:

$$
\begin{aligned}
\min_{z \in \mathbb{R}^{2N+1}}\quad & J(z) \\
\text{s.a.}\quad & v_N(z) = 0, \qquad x_N(z) = X \\
& v_k(z) \in [0,\ \bar v], \qquad k = 1, \dots, N-1 \\
& a_k^+ \in [0,\ \bar a], \qquad a_k^- \in [0,\ |\underline a|], \qquad T \ge 0
\end{aligned}
$$

Y con ello bastaría aplicar el algoritmo SLSQP (*Sequential Least Squares Programming*) para obtener la solución. A grandes rasgos, en cada iteración el algoritmo aproxima el problema con uno más sencillo (costo cuadrático y restricciones lineales), lo resuelve, avanza un paso en esa dirección y repite hasta que el costo deja de bajar y las restricciones se cumplen. En la práctica se usa la implementación de `scipy.optimize.minimize`.

Para empezar, el optimizador necesita un punto de partida. Se usa una semilla sencilla que ya cumple las condiciones de llegada: acelerar al máximo, ir a velocidad constante y frenar al máximo, sin tomar en cuenta los costos.

---

## 2.6 Implementación y resultados

El script `MM_conductor_optimo.py` resuelve los tres modelos con $N = 160$ tramos y guarda los resultados en `resultados_MM.xlsx`, con una hoja por modelo y las columnas $t$, $a$, $x$, $v$. La aceleración se define por tramos, así que en el último nodo queda vacía.

Con los parámetros base ($A = 0.5$, $B = 0.01$, $C = 0.001$, $D = 1$, $E = 10$, $\delta = 5$, $\bar a = 1$, $\underline a = -1.2$, $\bar v = 16.7$ m/s, $X = 800$ m):

| Modelo | Costo | $T$ [s] | $v_{\max}$ [m/s] |
|---|---|---|---|
| A: CMEM simétrico | 984.73 | 77.24 | 12.17 |
| B: asimétrico | 1047.29 | 86.21 | 10.44 |
| C: asimétrico + desgaste | 1176.04 | 104.30 | 9.49 |

Los parámetros son ilustrativos, no están calibrados. Aun así, se observa que mientras más realista es el costo de acelerar, más lento y más suave maneja el conductor: en el modelo B gana menos velocidad porque acelerar cuesta, y en el modelo C además suaviza las aceleraciones y las frenadas por el desgaste.
