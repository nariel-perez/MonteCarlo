#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 18:52:57 2023

@author: ariel
"""
import numpy as np

class gas_ideal():
    
    def __init__(self,name):
    
        self.N = 50 # número de partículas
        self.L = 10.0 # longitud del lado del cubo
        self.V = self.L**3 # volumen del cubo
        self.T = 270.0 # temperatura
        self.beta = 1/self.T # factor de Boltzmann
        self.sigma = 1.0 # parámetro del potencial Lennard-Jones
        self.epsilon = 1.0 # parámetro del potencial Lennard-Jones
        self.rc = 2.5*self.sigma # distancia de corte del potencial
        self.name = name
        self.steps = 1000
        self.E,self.frames = self.MC_gas(self.N, self.L, self.V, self.T, self.beta, self.sigma, self.epsilon, self.rc, self.steps)
    def LJ_potential(self, r, sigma, epsilon):
        """Potencial Lennard-Jones"""
        r6 = (sigma/r)**6
        return 4*epsilon*(r6**2 - r6)
    
    def energy(self,x, L, sigma, epsilon, rc):
        """Energía total del sistema"""
        e = 0.0
        for i in range(self.N-1):
            for j in range(i+1, self.N):
                dx = x[i,:] - x[j,:]
                dx = dx - self.L*np.round(dx/self.L)
                r = np.sqrt(np.sum(dx**2))
                if r < rc:
                    e += self.LJ_potential(r, self.sigma, self.epsilon)
        return e
    
    def MC_gas(self,N, L, V, T, beta, sigma, epsilon, rc, steps):#=100):
        """Simulación de gas ideal con algoritmo Metropolis-Hastings"""
        x = self.L*np.random.rand(self.N, 3) # posiciones iniciales
        e = self.energy(x, self.L, self.sigma, self.epsilon, self.rc) # energía inicial
        E = np.zeros(steps) # arreglo para guardar la energía
        accepted = 0 # número de movimientos aceptados
        #posiciones =np.zeros((N,3))
        frames = []
        count = 0
        frames.append(x)
        for i in range(steps):
            for i in range(N):
                count += 1
                # mover una partícula aleatoria
                j = np.random.randint(N)
                dx = self.L*(np.random.rand(3)-0.5)
                x_new = x.copy()
                #aqui condiciones de contorno
                x_new[j,:] += dx
                for h in range(3):
                    if x_new[j,h] > self.L/2:
                        x_new[j,h] -= self.L/2
                    elif x_new[j,h]< self.L/2:
                        x_new[j,h] += self.L/2
                     
                # calcular la energía del sistema modificado
                e_new = self.energy(x_new, self.L, self.sigma, self.epsilon, self.rc)
                # calcular la diferencia de energía y comparar con la probabilidad de aceptación
                delta_e = e_new - e
                if np.random.rand() < np.float128(np.exp(-beta*delta_e)):
                    x = x_new.copy()
                    e = e_new
                    accepted += 1
                    E[i] = e
                    if count%100 == 0:
                        frames.append(x)
        print("Tasa de aceptación: {:.2f}%".format(100*accepted/(N*steps)))
        return E, frames
    
    
    def gro(self,iniciales,N,L):
        with open('gas.gro', 'a') as out:
            out.write(f' Gasideal \n')
            out.write(f'\t\t{N}\n')
            system ='1gas'
            tipo ='GAS'
            for i,j in enumerate(iniciales):
                out.write(f'{system:^12s}{tipo:2s}{i+1:5d}{j[0]:8.3f}{j[1]:8.3f}{j[2]:8.3f}\n')
            
            out.write(f'   {L}\t{L}\t{L}\n')
        return ()
    





