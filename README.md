# IIQ2003_Seaweed_Drying

Proyecto del grupo 12 enfocado en simular el secado (difusión) en una lámina de algas usando un esquema de Crank–Nicolson.

Chile tiene más de 4.000 km de costa y presenta más de  800 especies de macroalgas, es por ello que nuestro país es uno de los principales productores y exportadores de macroalgas. Cuando se extraen del mar, presentan entre  90 y 95% de humedad, para evitar su degradación, deben pasar por un proceso secado, para asegurar calidad, estabilidad y valor comercial. 

Sobre todo en el norte, predomina el secado al sol; un método que es lento, genera variabilidad en la calidad, riesgo de contaminación y muy dependiente del clima. Por lo que éste método, no es el mejor para utilizar en el sur de Chile, dada la mayor presencia de humedad y frío que impide el secado natural. Para producir bioproductos de mayor valor, se necesita un método de secado más controlado y eficiente como lo es el secado por hornos.

El fenómeno de transporte predominante en el secado de algas mediante el aire (ya sea al aire libre o con estufas) es la transferencia de masa, específicamente la eliminación de agua desde el interior del alga hacia el medio o aire exterior. Este proceso depende intrínsecamente del tiempo, ya que la cantidad de humedad presente va cambiando continuamente hasta alcanzar un equilibrio con las condiciones del aire, por lo que para poder describir adecuadamente este problema se debe modelar un proceso transiente que considere las variaciones en el tiempo. 

La ecuación gobernante es la siguiente.
$$
\frac{\partial M}{\partial t} = D\,\frac{\partial^2 M}{\partial x^2}
$$

Sujeto a la condicion inicial

$$
M(x,0) = M_0
$$

Y a las condiciones de borde

- Condición de borde en el centro (x=0), flujo nulo :
$$
\left.\frac{\partial M}{\partial x}\right|_{x=0} = 0
$$

- Condición de borde en la superficie (x=L), convección:

$$
-D\left.\frac{\partial M}{\partial x}\right|_{x=L} = h\,\bigl(M(L,t) - M_{eq}\bigr)
$$

Donde 

- `M_0`: Valor inicial de la humedad en toda la lámina. (adimensional)
- `M_eq`: Humedad de equilibrio con el ambiente. (adimensional)
- `D`: Coeficiente de difusividad efectiva de la humedad en el material (m^2/s).
- `L`: Media Longitud de la lámina (m).
- `h`: Coeficiente de transferencia de masa convectiva en la superficie (aparece en la condición de Robin). Unidad: m/s (compatible con la formulación usada en el código).

Además para la discretización se va a usar
- `N_x`: Número de nodos espaciales en la discretización.
- `N_t`: Número de pasos temporales a simular.
- `dx = L / (N_x - 1)`: Paso espacial (m).
- `dt`: Paso temporal (s).
- `x_grid`: Grilla espacial
- `t_grid`: Grilla temporal
- `alpha = D*dt/dx**2`: número adimensional que mide la razón difusiva (aparece en las matrices del esquema numérico).
- `beta = 2*dx*h*M_eq/D`: término que aparece en la condición de borde en la formulación matricial.

**Discretización (esquema de diferencias finitas)**

El esquema de discretización utilizado fue diferencias finitas sobre la ecuación de difusión, donde al aplicar discretizaciones de segundo orden centradas en el espacio y de primer orden _forward_ en el tiempo. En el borde x=0 se usó una discretización de segundo orden _forward_ para la derivada espacial y en el  borde x=L se usó una discretización de segundo orden _backward_ para la derivada espacial; finalmente la discretización resulta:


- Ecuación de evolución:

$$
M_i^{\,k+1} = M_i^{\,k} + \frac{D\,\Delta t}{(\Delta x)^2}\;\bigl(M_{i-1}^{\,k} - 2\,M_i^{\,k} + M_{i+1}^{\,k}\bigr)
\qquad i=1,\dots,N_x-2,\; k=0,\dots,N_t-1
$$

- Condición inicial:

$$
M_i^{\,0} = M_0 \qquad \forall i \in \{0, \dots, N_x - 1\}
$$

- Condición de borde en $x=0$: derivada espacial nula.

$$
!-3\,M_0^{\,k} + 4\,M_1^{\,k} - M_2^{\,k} = 0
\qquad \forall k\in \{0,\dots,N_t-1\}
$$

Esta expresión corresponde a imponer $\partial M/\partial x(0,t)=0$ usando la aproximación
$$
\frac{-3M_0 + 4M_1 - M_2}{2\Delta x} \approx \left.\frac{\partial M}{\partial x}\right|_{x=0}.
$$

- Condición de borde en $x=L$ (CB2), condición de tipo Robin (intercambio convectivo):

$$
\Bigl(3 + \frac{2h\,\Delta x}{D}\Bigr)\,M_{N_x-1}^{\,k} - 4\,M_{N_x-2}^{\,k} + M_{N_x-3}^{\,k} = \frac{2h\,\Delta x}{D}\,M_{eq}
\qquad k=0,\dots,N_t-1
$$

Que se obtiene al aproximar la derivada en el extremo con una diferencia hacia atrás de segundo orden
$$
\left.\frac{\partial M}{\partial x}\right|_{x=L} \approx \frac{3M_{N_x-1} - 4M_{N_x-2} + M_{N_x-3}}{2\Delta x}
$$
y sustituir en la condición de flujo convectivo
$$
-D\,\left.\frac{\partial M}{\partial x}\right|_{x=L} = h\,\bigl(M(L,t) - M_{eq}\bigr).
$$

**Comentario sobre el esquema**
- La ecuación matricial que resume este sistema es $A @ M^{k+1} = B @ M^{k}$ con las matrices A y B definidas en `utils.py`
- Al utilizar el ezquema de Crank Nickolson resulta un paso temporal implícito con el sistema$A @ M^{k+1} = B @ M^{k}$, donde el producto del lado derecho se calcula directamente y luego se invierte el sistema lineal resultante.
- Dada la estructura _sparce_ de la matriz se utiliza un solver _sparce_ que reduce el tiempo de solución del paso temporal, lo que acelera el método.


**Contenido**
- `main.ipynb`: Notebook principal con ejemplo de uso.
- `utils.py`: Implementación de la clase `DryingModel` que resuelve la difusión y dibuja cortes temporales.

**Requisitos**
- Paquetes: `numpy`, `matplotlib`, `scipy`, `jupyter` o `ipykernel` en VS Code)


**Ejecutar el notebook**

Abrir `main.ipynb` y ejecutar las celdas en orden.


**Descripción de funcionalidades**

El repositorio sigue el paradigma de OOP y dentro del módulo utils.py se encuentra la clase DryingModel, que contiene todas las funcionalidades del proyecto. La clase DryingModel y todos sus métodos poseen documentación propia detallando los _input_ para su uso, una breve explicación  funcionalidad y su _output_ si es que lo genera o en que atributo se almacena.


