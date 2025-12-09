# IIQ2003_Seaweed_Drying

Breve proyecto para simular el secado (difusión) en una lámina usando un esquema de Crank–Nicolson.

**Contenido**
- `main.ipynb`: Notebook principal con ejemplo de uso.
- `utils.py`: Implementación de la clase `DryingModel` que resuelve la difusión y dibuja cortes temporales.

**Requisitos**
- Python 3.8+
- Paquetes: `numpy`, `matplotlib`, `scipy`, `jupyter` (u `ipykernel` si usas VS Code)

**Instalación rápida (PowerShell)**

```powershell
# (opcional) crear y activar un entorno virtual
python -m venv .venv; .\.venv\Scripts\Activate.ps1

# instalar dependencias
pip install --upgrade pip
pip install numpy matplotlib scipy jupyter
```

**Ejecutar el notebook**

- Desde VS Code: abrir `main.ipynb` y ejecutar las celdas en orden.
- Desde la terminal (PowerShell):

```powershell
jupyter notebook
# luego abrir main.ipynb en el navegador
```

Recomendación de ejecución en el notebook (celdas en orden):
- 1) Importar `DryingModel` (celda 1)
- 2) Definir parámetros y crear `model` (celda 2)
- 3) Ejecutar `model.simulate_diffusion()` y `model.plot_cuts([...])` (celda 3)
- 4) (Opcional) `model.validate_model()` si existe esa comprobación en tu versión.

**Notas importantes / Observaciones**

- Arreglo aplicado: se detectó un `IndexError: index 50 is out of bounds for axis 1 with size 50` al pedir cortes en `x = 1.0` (fracción 1.0 de la longitud). Esto ocurría por inconsistencia en el cálculo del índice de la malla: en algunos lugares se usaba `int(x0 * self.N_x)` (que para `x0=1.0` devuelve `N_x`, fuera de rango) y en otros `int(x0 * (self.N_x - 1))`. En la versión actual se normalizó el cálculo usando `N_x - 1` y se clampa el índice para garantizar que siempre esté en el rango `[0, N_x-1]`.

- Línea potencialmente incorrecta en `utils.py`: en el constructor aparece `self.D = D-9`. Es muy probable que se trate de un error tipográfico. Posibles correcciones:
  - Si `D` se pasa en unidades de m^2/s ya escaladas (por ejemplo `1.0e-9`), entonces debería ser `self.D = D`.
  - Si la intención era convertir de `1.0` a `1e-9`, debería usarse `self.D = D * 1e-9`.

  Si quieres, puedo cambiar esto por la alternativa que prefieras; por ahora no lo toco automáticamente para no romper resultados inesperadamente.

**Sugerencias**
- Revisar la línea con `self.D = D-9` y confirmar la intención (te ofrezco corregirla si me indicas qué comportamiento deseas).
- Si planeas distribuir o reproducir resultados, añade un `requirements.txt` o `pyproject.toml` con versiones fijas.

**Contacto / Próximos pasos**
- Si quieres, ejecuto el notebook para verificar que ya no aparece el IndexError y subo una celda de prueba que muestre los cortes.
- También puedo actualizar `utils.py` para corregir `self.D` si confirmas cuál es el comportamiento esperado.

---

Archivo generado automáticamente: `README.md` (contenido básico). Si deseas ajustes de idioma, formato o más secciones (por ejemplo: explicación matemática, referencia bibliográfica, ejemplos de gráficos), dime y lo amplío.