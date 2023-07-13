import numpy as np
import matplotlib.pyplot as plt

def bf(alpha, chi, N_m, N_ch):
    '''
    Calculo de la energía 'interna' del gel... energía atribuida al swelling del mismo.
    '''
    term1 = (alpha**3 - 1) * np.log(1 - alpha**(-3))
    term2 = chi * (1 - alpha**(-3))
    term3 = (3/2) * N_ch * (alpha**2 - np.log(alpha) - 1)
    energy = N_m * (term1 + term2) + term3
    return energy

def interaccion(ai, aj, r, epsilon):
    '''
    Potencial de Hertz para dos particulas de radiso ai y aj a una distancia r
    '''
    if r < (ai + aj):
        return epsilon * (1 - (r / (ai + aj))**(5/2))
    else:
        return 0
    
def pair_energy(radios, posiciones, epsilon, L):
    '''
    energía total del sistema de N partículas. Usando el potencial de Hertz
    '''
    energia_total = 0.0 
    N = len(radios)
    for i in range(N):
      for j in range(i+1,N):
        dx = posiciones[i,:] - posiciones[j,:]
        dx = -L*np.round(dx/L)
        r = np.sqrt(np.sum(dx**2))
        ai = radios[i]
        aj = radios[j]
        energia_total += interaccion(ai, aj, r, epsilon)
    
    return energy_total

def simulacion(steps, L, a0, epsilon, N_m, chi, N_ch, K, T):
    N = 500
    beta = 1 / (K * T) # No necesario 
    energyt = []
    alphas_values = []
    # podría ser otro valor, por lo que sería mejor: bf0 = [bi]*N , en donde bi es la energía inicial (distinta de cero)  
    bf0 = np.zeros(N)
    a0_array = np.full(N, a0) # valores iniciales de radios
    #cambio a valores de radio aleatorio
    for i in range(N):
        alpha = np.random.uniform(1.1, 10.0)
        a0_array[i] *= alpha
        bf0[i] +=  bf(alpha, chi, N_m, N_ch)    
    
    
    posiciones = L * np.random.rand(N, 3) #posiciones iniciales en la caja
    # Copia de los radios, para guardar los cambios nuevos
    radios = a0_array.copy()
    #energias iniciales, no hay swelling. La energía de cada partícula es cero... casualmente
   
    #calculo de la energía inicial por la interacción de a pares. 
    vh0 = pair_energy(radios, posiciones, epsilon, L)
    ee = np.sum(bf0) + vh0 #energía total inicial... sin utilidad en este punto.(?)
    print(vh0)
    print(ee)
    
    paso = 0
    acc = 0
    
    for i in range(steps):
        for p in range(N):
            paso += 1
            j = np.random.randint(N) #elección partícula j al azar
            dx = L * (np.random.rand(3) - 0.5) # delta de movimiento
            alpha = np.random.uniform(1.1, 10.0) #elección de un grado de swelling al azar
            #copiado para las posibles nuevas configuraciones
            x_new = posiciones.copy()
            radios_new = radios.copy()
            bf_new = bf0.copy()
            # Modificación de las variables
            x_new[j, :] += dx
            #Condiciones de contorno:
            '''
            SINTAXIS: a = np.where(condicion, x,y)
            Si la condiciones es True: el valor de a es x, de ser False le asigna un valor de y.
            Es decir si al mover la particula j se sale de la mitad de la caja (L/2) ajusta la variable a x_new[j,:] - L/2
            de no salirse (condición False) no ajusta nada.
            Del mismo modo se compara cuando la condición sea ser menor que L/2
            '''
            x_new[j, :] = np.where(x_new[j, :] > L/2, x_new[j, :] - L/2, x_new[j, :])
            x_new[j, :] = np.where(x_new[j, :] < L/2, x_new[j, :] + L/2, x_new[j, :])
            
            radios_new[j] *= alpha
            # Calculo de las nuevas energías. 
            bf_i = bf(alpha, chi, N_m, N_ch)
            vh_i = pair_energy(radios_new, x_new, epsilon, L)
            Df = bf_i - bf0[j]
            Dvh = vh_i - vh0
            Dtotal = Df + Dvh
            print(Dtotal)
            
            if Dtotal < 0 or np.random.rand() < np.exp(-beta * Dtotal):
                acc += 1
                posiciones = x_new.copy()
                radios = radios_new.copy()
                vh0 = vh_i
                bf0 = bf_new.copy()
                ee = np.sum(bf0) + vh0
                alphas_values.append(alpha)
                energyt.append(ee)
            
            if paso % 50 == 0:
                energyt.append(ee)
                print(f'paso: {paso}, energia total: {ee:.2f}')
            
            
    
    print(paso, acc)
    
    return energyt, alphas_values, energy_values

def histo_alphas(alphas_values,energy_values):
  plt.hist(alphas_values, bins=20, density=True, weights=energy_values)
  plt.xlabel('Alpha')
  plt.ylabel('Probabilidad')
  plt.title('Distribución de Probabilidad de Alpha')
  plt.show()

# Parámetros globales
a0 = 10.0
V = 45110592.60
L = V**(1/3)
epsilon = 1.5e3
N_m = 2e5
chi = 0.0
N_ch = 200
K = 1.380649e-23
T = 300

# Ejecutar simulación
energyt, alphas_values, energy_values = simulacion(1000, L, a0, epsilon, N_m, chi, N_ch, K, T)

# histograma
#histo_alphas(alphas_values)



