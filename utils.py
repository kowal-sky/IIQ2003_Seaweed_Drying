import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.stats import linregress



class DryingModel:

    '''
    Docstring for DryingModel
    Clase principal del proyecto, contiene todos los parametros del modelo
    y los metodos para simular la difusion de humedad en una muestra
    y utiliza el esquema de Crank-Nicolson con el paso de A @ u_next = B @ u_prev
    para resolver la ecuacion de difusión de humedad en una muestra unidimensional.

    Considera condiciones de frontera de flujo de masa en x=L (Robin) y 
    simetría en x=0 (Neumann)..
    
    '''


    def __init__(self, M_0, M_eq, D, L, N_x, N_t, h):

        '''
        Docstring for __init__
        
        :param M_0: Humedad inicial (adimensional)
        :param M_eq: Humedad de equilibrio (adimensional)
        :param D: Difusividad efectiva (m^2/s)
        :param L: Semi Grosor de la muestra (m)
                  (mitad del espesor total)
        :param N_x: Numero de puntos en el espacio
        :param N_t: Numero de pasos en el tiempo
        :param h: Coeficiente de transferencia de masa (m/s)

        '''

        self.M_0 = M_0        
        self.M_eq = M_eq
        self.N_x = N_x
        self.L = L       
        self.dx = L / (N_x - 1)
        self.D = D
        self.h = h
        self.x_grid = np.linspace(0, L, N_x)
        self.dt = 10 * (self.dx**2) / D  #Fourier con Alpha=10
        self.N_t = N_t
        self.t_grid = np.linspace(0, self.N_t * self.dt, self.N_t+1)
        self.alpha = None
        self.beta = None
        self.A = None
        self.B = None
        self.u_history = None
        self.sistem_solved = False

    def build_A(self, sparse_format = 'csc'):
        '''
        Docstring for build_A
        
        :param sparse_format: Formato de la matriz dispersa (default 'csc')
        :return: Matriz A del esquema de Crank-Nicolson
        
        La matriz A se guarda en el atributo self.A
        '''

        dx = self.dx

        alpha = self.alpha

        A1 = sp.eye(self.N_x, k = 0, format=sparse_format) * (2*(1 + alpha))
        A2 = sp.eye(self.N_x, k = -1, format=sparse_format) * -alpha
        A3 = sp.eye(self.N_x, k = 1, format=sparse_format) * -alpha
        A = A1 + A2 + A3
        
        # Boundary condition at x = 0
        A[0, 0] = -3; A[0, 1] = 4; A[0, 2] = -1

        # Boundary condition at x = L
        A[-1, -1] = 3 + 2 * dx * self.h / self.D; A[-1, -2] = -4; A[-1, -3] = 1

        self.A = A

    def build_B(self, sparse_format = 'csc'):
        '''
        Docstring for build_B
        
        :param sparse_format: Formato de la matriz dispersa (default 'csc')
        :return: Matriz B del esquema de Crank-Nicolson

        La matriz B se guarda en el atributo self.B
        '''

        alpha = self.alpha

        B1 = sp.eye(self.N_x, k = 0, format=sparse_format) * (2*(1 - alpha))
        B2 = sp.eye(self.N_x, k = -1, format=sparse_format) * alpha
        B3 = sp.eye(self.N_x, k = 1, format=sparse_format) * alpha
        B = B1 + B2 + B3
        B[0, 0] = 0; B[0, 1] = 0; B[-1, -1] = 0; B[-1, -2] = 0

        self.B = B


    def Crank_Nicolson_step(self, u_prev):

        '''
        Docstring for Crank_Nicolson_step
        
        :param u_prev: Vector de humedad en el paso de tiempo previo
        :return: Vector de humedad en el siguiente paso de tiempo
        Realiza un paso del esquema de Crank-Nicolson
        utilizando las matrices A y B ya construidas
        '''

    
        b = self.B @ u_prev; b[-1] = self.beta

        # Solve A @ u_next = B @ u_prev
        u_next = sp.linalg.spsolve(self.A, b)

        return u_next

    def simulate_diffusion(self, sparse_format = 'csc'):

        '''
        Docstring for simulate_diffusion
        
        :param sparse_format: Formato de la matriz dispersa (default 'csc')
        Realiza la simulacion de la difusion de humedad
        utilizando el esquema de Crank-Nicolson guardando
        los resultados en el atributo self.u_history
        '''

        dx = self.dx
        self.alpha = self.D * self.dt / dx**2
        self.beta = 2 * self.dx * self.h * self.M_eq / self.D

        # # Construct matrices A and B

        self.build_A(sparse_format)
        self.build_B(sparse_format)

        # Initial condition
        u = np.ones(self.N_x) * self.M_0

        # results history
        results = np.zeros((self.N_t + 1, self.N_x))
        results[0, :] = u.copy()

        for n in range(self.N_t):
            u = self.Crank_Nicolson_step(u)
            # results.append(u.copy())
            results[n + 1, :] = u.copy()

        self.u_history = np.array(results)
        self.sistem_solved = True
    
    
    def plot_cuts(self, x_points):

        '''
        Docstring for plot_cuts
        
        :param x_points: Lista de posiciones x (en fraccion de L) donde se quieren hacer los cortes
        Realiza cortes de humedad en las posiciones especificadas y grafica
        la evolucion temporal de la humedad en estos puntos
        '''


        if not self.sistem_solved:
            raise RuntimeError("The diffusion simulation has not been run yet.")
        
        plt.figure(figsize=(10, 6))
        for x0 in x_points:
            if x0 < 0 or x0 > 1:
                raise ValueError(f"Point {x0} L is out of bounds.")
            # Compute index on grid consistently using (N_x - 1) as max index
            ix = int(round(x0 * (self.N_x - 1)))
            # Clamp index to valid range [0, N_x-1]
            ix = min(max(ix, 0), self.N_x - 1)
            print(self.u_history.shape)
            print(self.u_history[:, ix])
            plt.plot(self.t_grid / 3600, self.u_history[:, ix], label=f'x={x0} L')
        
        plt.xlabel('Position (x)')
        plt.ylabel('Moisture Content (M)')
        plt.title('Moisture Content Profiles at Different Times')
        plt.legend()
        plt.grid()
        plt.show()


    def global_process(self):

        '''
        Docstring for global_process
        
        Grafica la curva de secado global (humedad promedio vs tiempo)
        y la curva logaritmica para validacion del modelo
        '''

        if not self.sistem_solved:
            raise RuntimeError("The diffusion simulation has not been run yet.")
        
        M_average = (np.sum(self.u_history, axis=1) - 0.5*self.u_history[:,0] - 0.5*self.u_history[:,-1]) / (self.N_x - 1) # Regla del trapecio
        

        plt.figure(figsize=(8, 6))
        plt.plot(self.t_grid / 3600, M_average, 'b-', linewidth=2, label='Simulación FDM')
        plt.xlabel('Tiempo (h)')
        plt.ylabel('Humedad Promedio (kg/kg)')
        plt.title('Curva de Secado Global')
        plt.grid(True)
        plt.legend()
        plt.show()

        MR = (M_average - self.M_eq) / (self.M_0 - self.M_eq) # Arranaza -->MR = (M_t - M_eq) / (M_0 - M_eq)

        valid_mask = MR > 0.001
        ln_MR = np.log(MR[valid_mask])
        time_log = self.t_grid[valid_mask]

        plt.figure(figsize=(8, 6))
        plt.plot(time_log / 3600, ln_MR, 'r-')
        plt.xlabel('Tiempo (h)')
        plt.ylabel('ln(MR)')
        plt.title('Curva Logarítmica (Validación)')
        plt.grid(True)
        plt.show()

    def validate_model(self):

        '''
        Docstring for validate_model
        
        Realiza validación de la solucion mediante regresión lineal
        de la curva logarítmica de humedad promedio vs tiempo

        Método propuesto por Arranaza et al.

        Arranza, F. J., Jiménez-Ariza, T., Diezma, B., & Correa, E. C. (2017). 
        Determination of diffusion and convective transfer coefficients in food drying revisited: 
        A new methodological approach. Biosystems Engineering, 162, 30-39. 

        https://doi.org/10.1016/j.biosystemseng.2017.07.005 
        '''
        if not self.sistem_solved:
            raise RuntimeError("The diffusion simulation has not been run yet.")

        M_average = (np.sum(self.u_history, axis=1) - 0.5*self.u_history[:,0] - 0.5*self.u_history[:,-1]) / (self.N_x - 1) # Regla del trapecio

        MR = (M_average - self.M_eq) / (self.M_0 - self.M_eq) # Arranaza -->MR = (M_t - M_eq) / (M_0 - M_eq)

        MR = (M_average - self.M_eq) / (self.M_0 - self.M_eq) # Arranaza -->MR = (M_t - M_eq) / (M_0 - M_eq)

        valid_mask = MR > 0.001

        time_log = self.t_grid[valid_mask]
        ln_MR = np.log(MR[valid_mask])

        mitad = len(time_log) // 2 #regresión solo en la parte recta

        t_cola = time_log[mitad:]
        ln_MR_cola = ln_MR[mitad:]

        slope, intercept, r_value, p_value, std_err = linregress(t_cola, ln_MR_cola)

        print("=== RESULTADOS DE LA VALIDACIÓN (ARRANZA) ===")
        print(f"Pendiente (b): {slope:.6f}")
        print(f"Intercepto (a): {intercept:.6f}")
        print(f"R-cuadrado: {r_value**2:.6f}")