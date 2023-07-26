#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 09:49:38 2023

@author: ariel
"""

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
    Potencial de Hertz para dos partículas de radios ai y aj a una distancia r
    '''
    
    mask = r < (ai + aj)
    return np.where(mask, epsilon*(1 - (r / (ai + aj))**(5/2)),0)

def pair_energy(radios, posiciones, epsilon, L):
    '''
    Energía total del sistema de N partículas. Usando el potencial de Hertz
    '''
    N = len(radios)
    distancias = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            dx = posiciones[i, :] - posiciones[j, :]
            dx = -L * np.round(dx / L)
            distancias[i, j] = np.linalg.norm(dx)

    # Utilizamos numpy.triu_indices para recorrer sólo la mitad superior de la matriz de distancias
    indices = np.triu_indices(N, k=1)
    r_values = distancias[indices]
    ai_values = radios[indices[0]]
    aj_values = radios[indices[1]]

    energia_total = np.sum(interaccion(ai_values, aj_values, r_values, epsilon))
    return energia_total


def simulacion_denton(steps, L, a0, epsilon, N_m, chi, N_ch, K, T):
    N = 500
    p_totales = steps*N
    beta = 1 / (K * T) # No necesario 
    energyt = []
    alphas_values = []
    # podría ser otro valor, por lo que sería mejor: bf0 = [bi]*N , en donde bi es la energía inicial (distinta de cero)  
    bf0 = np.zeros(N) + bf(3.5, chi, N_m, N_ch)
    a0_array = np.full(N, a0)*3.5 # valores iniciales de radios   
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
    for p in range(N)
        paso += 1
        j = np.random.randint(N) #elección partícula j al azar
        dx = L * (np.random.rand(3) - 0.5) # delta de movimiento
        rn = np.random.uniform(0.0, 1.0) #elección de un grado de swelling al azar
        alpha += (rn -0.5)*0.02 
        alpha = np.where(alpha < 1.0, 1.1, alpha)
        x_new = posiciones.copy()
        radios_new = radios.copy()
        bf_new = bf0.copy()
        x_new[j, :] += dx
        x_new[j, :] = np.where(x_new[j, :] > L/2, x_new[j, :] - L/2, x_new[j, :])
        x_new[j, :] = np.where(x_new[j, :] < L/2, x_new[j, :] + L/2, x_new[j, :])
        radios_new[j] *= alpha
        bf_i = bf(alpha, chi, N_m, N_ch)
        vh_i = pair_energy(radios_new, x_new, epsilon, L)
        Df = bf_i - bf0[j]
        Dvh = vh_i - vh0
        Dtotal = Df + Dvh
            
        if Dtotal < 0 or np.random.rand() < np.exp(-Dtotal):
            acc += 1
            posiciones[j,:] += dx
            radios[j] *= alpha
            vh0 += Dvh
            bf0[j] += Df
            
            ee += Dtotal
            
            if paso > 500:
                alphas_values.append(alpha)
                energyt.append(ee)           
            
            if paso % 500 == 0:
                porcentaje = (paso/p_totales)*100
                print(f'paso: {paso}, energia total: {ee:.2f}')
                print(f' Porcentaje de simulacion: {porcentaje:.2f}')
    
    print(paso, acc)
    
    return energyt, alphas_values


def histo_alphas(alphas_values, bins):
  plt.hist(alphas_values, bins=bins)
  plt.xlabel('Alpha')
  plt.ylabel('Probabilidad')
  plt.title('Distribución de Probabilidad de Alpha')
  plt.show()



a0 = 10.0
V = 45110592.60
L = V**(1/3) # 355.98
N_m = 2e5
chi = 0.0
N_ch = 200
K = 1.380649e-23
T = 300
epsilon = 1.5e3

enegyt, alphas_values = simulacion_denton(20, L, a0, epsilon, N_m, chi, N_ch, K, T)
histo_alphas(alphas_values,50)
