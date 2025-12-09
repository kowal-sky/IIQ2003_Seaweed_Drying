# IIQ2003_Seaweed_Drying

Proyecto del grupo 12 enfocado en simular el secado (difusión) en una lámina de algas usando un esquema de Crank–Nicolson.

Chile tiene más de 4.000 km de costa y presenta más de  800 especies de macroalgas, es por ello que nuestro país es uno de los principales productores y exportadores de macroalgas. Cuando se extraen del mar, presentan entre  90 y 95% de humedad, para evitar su degradación, deben pasar por un proceso secado, para asegurar calidad, estabilidad y valor comercial. 

Sobre todo en el norte, predomina el secado al sol; un método que es lento, genera variabilidad en la calidad, riesgo de contaminación y muy dependiente del clima. Por lo que éste método, no es el mejor para utilizar en el sur de Chile, dada la mayor presencia de humedad y frío que impide el secado natural. Para producir bioproductos de mayor valor, se necesita un método de secado más controlado y eficiente como lo es el secado por hornos.

El fenómeno de transporte predominante en el secado de algas mediante el aire (ya sea al aire libre o con estufas) es la transferencia de masa, específicamente la eliminación de agua desde el interior del alga hacia el medio o aire exterior. Este proceso depende intrínsecamente del tiempo, ya que la cantidad de humedad presente va cambiando continuamente hasta alcanzar un equilibrio con las condiciones del aire, por lo que para poder describir adecuadamente este problema se debe modelar un proceso transiente que considere las variaciones en el tiempo. 

**Contenido**
- `main.ipynb`: Notebook principal con ejemplo de uso.
- `utils.py`: Implementación de la clase `DryingModel` que resuelve la difusión y dibuja cortes temporales.

**Requisitos**
- Paquetes: `numpy`, `matplotlib`, `scipy`, `jupyter` o `ipykernel` en VS Code)


**Ejecutar el notebook**

Abrir `main.ipynb` y ejecutar las celdas en orden.


**Descripción de funcionalidades**

El repositorio sigue el paradigma de OOP y dentro del módulo utils.py se encuentra la clase DryingModel, que contiene todas las funcionalidades del proyecto. La clase DryingModel y todos sus métodos poseen documentación propia detallando los _input_ para su uso, una breve explicación  funcionalidad y su _output_ si es que lo genera o en que atributo se almacena.


