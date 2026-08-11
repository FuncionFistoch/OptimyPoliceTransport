# 1. El Metropolitano Monopolístico 
El caso más básico posible, para plantear la modelización del flujo en las carreteras de una ciudad, es suponer una larga carretera en línea recta, donde exista un único transporte que de el servicio en toda la ruta. Este caso puede encontrar ejemplos en las vías de tren y en el Metropolitano, ambos casos son similares dado que otorga el derecho a un solo transporte de realizar la actividad económica. 

Por lo que se continua con el problema de optimización para el único conductor, para todo el recorrido se supondrá una distancia total $X$, con $n$ paraderos, se denominara la distancia entre el $k$ paradero y el $k+1$ como $X_k$. 

En cada paradero $k$, exceptuando el primero y último, tendrá dos variables aleatorias $u_k(t) \sim \text{Poisson}(\mu_{uk}(t))$ y $d_k(t) \sim \text{Poisson}(\mu_{dk}(t))$. Ambas variables representan la cantidad de personas que suben (o bajan) en el $k$ paradero en el instante $t$, estas variables dependen del tiempo debido a la existencia de "horas picos", dado que tanto en la entrada laboral como la salida existe una enorme cantidad de personas que llegan a diferencia de otras horas del día. Por lo que se capta las diferencias temporales de la demanda de transporte de esta manera. 

Cabe mencionar que se elige la distribución de Poisson debido a que cumple en cierta medida con los requisitos para ser una distribución de Poisson: 
1. La llegada de cada individuo es individual.
2. La llegada de cada individuo es independiente de la llegada de otro individuo.
3. La probabilidad de llegada de un individuo al paradero depende proporcionalmente del rango de tiempo.

Además que se sigue de la literatura que modela el comportamiento de los pasajeros. 

Además también deberá de tomarse en cuenta una variable de aforo máximo $\bar{A}$ y con ello una variable $A_t$ que evolucione con la subida/bajada de las pasajeros. 

$$
A_{t+1} =
\begin{cases}
A_t + u_t - d_t, & A_{t+1} < \bar{A}\\
A_t - d_t, & A_{t+1} = \bar{A}
\end{cases}
$$

Si bien esto define a la demanda de transporte, la parte de la oferta del transporte (o como mínimo la decisiones de como transportar a los pasajeros) esta definida por parte del conductor. Asumimos una función de coste a su servicio que puede ser desglosada en: El salario al conductor; que puede ser fijo en todo el trayecto o puede ser a destajo y depender de una porción de la cantidad recaudada en el recorrido, el consumo de gasolina o el combustible que permite dar el servicio y el deterioro de la carrocería. 

Se ahonda ahora el consumo de gasolina y como se modela según las variables de estado y de control que contempla el conductor: $x(t)$, $v(t)$ y $a(t)$. De esto y siguiendo la literatura sobre el modelo CMEM del consumo de gasolina de un transporte, tenemos la siguiente función $f(v, a)$: 

$$
f(v,a) = A \cdot v(t) + B \cdot v^2(t) + C \cdot v^3(t) + D \cdot a(t) \cdot v(t)  
$$

Esta función CMEM es mucho más compleja al introducir el factor de enriquecimiento del trabajo de investigación original, sin embargo, por simplicidad y dado que depende más de las condiciones del motor al arrancar, se considerara que tal factor es unitario. Los coeficientes de la función además tienen una interpretación física, al realizar el diagrama de cuerpo libre (y suponiendo que la carretera tiene pendiente constante e igual a 0) puede llegar a la misma ecuación, salvo por el factor cuadrático de la velocidad que se agrega para captar los efectos que pasan desapercibidos por la dinámica. Cabe mencionar también que estos coeficientes representan además el precio relativo que tiene con respecto del salario del conductor y el deterioro de la carrocería, debido a que nos encontramos ante una función de coste que incorpora varios insumos. 

También cabe mencionar que existe un factor importante que no esta contemplado en la función CMEM, y es el desgaste por cambios bruscos en la aceleración (variable de control para el conductor), que se incorpora mediante la función $g(a)$, su fundamento se encuentra en el cambio brusco de energía cinética en un pequeño tiempo. No se tomará en cuenta otra descomposición polinomial debido a que si bien imponer ad hoc una función polinómica haría más sencillo el trabajo, la caja negra de la función no explicaría el porque del valor de tales coeficientes, en cuanto el modelo propuesto derivado de la dinámica tiene un significado físico por lo que su calibración es más adecuada. 

$$
g(a) = E \cdot a^2(t) 
$$

Por lo que, suponiendo una función $W$ como la retribución al conductor por su trabajo, del cual se analizara con mayor ahinco después, se tiene entonces el costo total en el recorrido, que toma un tiempo $T$ en finalizar, como: 

$$
CT = W + \int_0^T  A \cdot v(t) + B \cdot v^2(t) + C \cdot v^3(t) + D \cdot a(t) \cdot v(t) + E \cdot a^2(t)dt
$$

Suponiendo que el conductor llega a cada $k$ paradero en el instante $T_k$ podemos definir una función de Ingresos esperada para el trayecto, suponiendo un precio de subida constante $p$ por simplicidad (una tarifa unica es de hecho el caso particular para el metropolitano como para muchas lineas de tren). 

$$
IT = \mathbb{E}\left[\sum_{k=0}^{n} p \cdot u_k(T_k)\right]
$$

Por lo que finalmente podemos definir una función de beneficios el cual optimizar, sujeto a las restricciones físicas.

$$
\max_{a(t)} J =\mathbb{E}\left[\sum_{k=0}^{n} p \cdot u_k(T_k)\right]- W - \int_0^T \left(f(a,v)+g(a)\right)dt
$$

# 1.1 El caso aún más simple: Minimizar costos

Supongamos, antes de resolver el problema estocástico, supongamos que el conductor desconoce completamente la distribución de los pasajeros, por lo que el problema le es irresoluble y prefiere minimizar el costo que asumirá al realizar el servicio. Debido a que el servicio se particiona entre dos paraderos, en donde el conductor tiene que detenerse para subir y bajar a pasajeros, entonces imponemos la restricción de $v(T_k) = 0$ y $x(T_k) = X_k$, además supondremos que el salario depende del tiempo, por lo que la función $W = \delta \cdot t$. Y de esto podemos tomar la siguiente partición al problema de minimización de costes: 

$$
\min_{a(t)} J = \int_{0}^{T_1} \delta+\left(f(a,v)+g(a)\right)dt + \int_{T_1}^{T_2}\delta+\left(f(a,v)+g(a)\right)dt+\cdots+\int_{T_{n-1}}^{T}\delta+\left(f(a,v)+g(a)\right)dt
$$

En esencia tendremos la estrategia de optimizar para cada partición $[T_{k}, T_{k+1}]$ e ir solapando las variables tanto de control como de estado para obtener un recorrido óptimo. Antes de presentar el problema de control óptimo se debe de aclarar las restricciones físicas que se tienen: 

1. Tanto $\dot{x}(t) = v(t)$ y $\dot{v}(t) = a(t)$ son el punto de partida para realizar el problema, se supondra que tanto la variable $x(t)$ y $v(t)$ son variables de estado que dependen de la variable de control $a(t)$, esto con el razonamiento de que el conductor toma la decisión de acelerar o frenar continuamente, y con ello elige una cierta velocidad.
2. Para la variable de control se supone la restricción de acotación de la siguiente forma: $a(t) \in [\underline{a},\overline{a}]$, que incluye al 0, por lo que el regimen negativo es donde el conductor frena y el regimen positivo es donde acelera. También en toda carretera tenemos una cota en la que se dispone que puede recorrer un vehículo, por lo que también se pone como restricción $v(t) \in [0,\overline{v}]$.
3. Como ya se había mencionado antes, el conductor debe detenerse en un lugar especifico para abordar pasajeros, por lo que, $v(T_k) = 0$ y $x(T_k) = 0$.

Por lo que finalmente, para el caso más simple donde solo se tiene un único paradero y no se tome en cuenta el gasto por cambios bruscos en la aceleración:

$$
\min_{a(t)} J = \int_{0}^{T} \delta+\left(f(a,v)\right)dt
$$

$$
s.a \quad \dot{v}(t) = a(t), \dot{x}(t) = v(t)
$$

Se define el hamiltoniano incorporando las variables de coestado:

$$
H = \delta + A \cdot v + B \cdot v^2 + C \cdot v^3 + D\cdot a \cdot v + \lambda_x \cdot v + \lambda_v \cdot a 
$$

Resolver este problema de control óptimo tiene la particularidad en la que se tiene una estructura bang-singular-bang o bang-bang del control, debido a que la linealidad del control sobre el hamiltoniano, hace depender del término $D \cdot v + \lambda_v$ para resolver el problema. Además que se tienen dos tipos de estrategia para el conductor: Acelerar, llegar a una velocidad de crucero, desacelerar, o simplemente acelerar y luego desacelerar. Estas dos estrategias dependerán de los parametros para obtener el más óptimo. 

De los archivos "MMsolution_bb.py" y "MMsolution_bSb.py", de donde se utiliza la librería SymPy para resolver el problema según el tipo de estrategia plausible, se puede comparar las soluciones. Si usted juega con los parametros, observara que para muchos casos no existe una solución para el regimen bang-singular-bang, y existe una preferencia por la estrategia bang-bang, que no parece realmente realista, por lo que se debe de afinar la función de costo, y esto no implica agregar la función de costo por el cambio abrupto de aceleración, sino que la función CMEM tiene por problematica que si se desacelera no se espera tener un ahorro (a menos que exista un mecanismo como en los carros eletricos que desacelerar si corresponde a un ahorro) por lo que la función se debe de ajustar a: 

$$
f(a,v) = \delta + A \cdot v + B \cdot v^2 + C \cdot v^3 + D\cdot max[a, 0] \cdot v
$$

Y con ello, en el archivo "MMsolution_desaceleratenocost.py" obtenemos una solución general para ambos regimenes, y del cual podemos estudiar según gráficos, dado que obtener una solución analitica es harto dificil, para comprender como el efecto del salario y la distancia total afecta la decisión de la decisión de la aceleración.
