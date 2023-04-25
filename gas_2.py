#%%
import numpy as np
import matplotlib.pyplot as plt


N = 50 # número de partículas
L = 10.0 # longitud del lado del cubo
V = L**3 # volumen del cubo
T = 1.0 # temperatura
beta = 1/T # factor de Boltzmann
sigma = 1.0 # parámetro del potencial Lennard-Jones
epsilon = 1.0 # parámetro del potencial Lennard-Jones
rc = 2.5*sigma # distancia de corte del potencial


def LJ_potential(r, sigma, epsilon):
    """Potencial Lennard-Jones"""
    r6 = (sigma/r)**6
    return 4*epsilon*(r6**2 - r6)

def energy(x, L, sigma, epsilon, rc):
    """Energía total del sistema"""
    e = 0.0
    for i in range(N-1):
        for j in range(i+1, N):
            dx = x[i,:] - x[j,:]
            dx = dx - L*np.round(dx/L)
            r = np.sqrt(np.sum(dx**2))
            if r < rc:
                e += LJ_potential(r, sigma, epsilon)
    return e

def MC_gas(N, L, V, T, beta, sigma, epsilon, rc, steps):#=100):
    """Simulación de gas ideal con algoritmo Metropolis-Hastings"""
    x = L*np.random.rand(N, 3) # posiciones iniciales
    e = energy(x, L, sigma, epsilon, rc) # energía inicial
    E = np.zeros(N*steps) # arreglo para guardar la energía
    accepted = 0 # número de movimientos aceptados
    #posiciones =np.zeros((N,3))
    frames = []
    frames.append(x)
    count = 0
    for i in range(steps):
        for p in range(N):
            count += 1
            # mover una partícula aleatoria
            j = np.random.randint(N)
            dx = L*(np.random.rand(3)-0.5)
            x_new = x.copy()
            #aqui condiciones de contorno
            x_new[j,:] += dx    
            for h in range(3):
                if x_new[j,h] > L/2:
                    x_new[j,h] -= L/2
                elif x_new[j,h]< L/2:
                    x_new[j,h] += L/2
                    
            # calcular la energía del sistema modificado
            e_new = energy(x_new, L, sigma, epsilon, rc)
            # calcular la diferencia de energía y comparar con la probabilidad de aceptación
            E[count-1] = e
            delta_e = e_new - e
            if np.random.rand() < np.float128(np.exp(-beta*delta_e)):
                x = x_new.copy()
                e = e_new
                accepted += 1
                #E[i] = e
                if count%100 ==0:
                    frames.append(x)
        
    print("Tasa de aceptación: {:.2f}%".format(100*accepted/(N*steps)))
    return E, frames


E,frames = MC_gas(N, L, V, T, beta, sigma, epsilon, rc, steps=100)

#%% 2
def trayectorias(frames):
    with open('gas.xyz', 'w') as dat:
        for frame in frames:
            dat.write(f'\t{N}\n\n')
            for j in frame:
                dat.write(f'C\t{j[0]}\t{j[1]}\t{j[2]}\n')


def gro(iniciales,N):
    with open('gas.gro', 'a') as out:
        out.write(f' Gasideal \n')
        out.write(f'\t\t{N}\n')
        system ='1gas'
        tipo ='GAS'
        L =10.0
        for i,j in enumerate(iniciales):
            out.write(f'{system:^12s}{tipo:2s}{i+1:5d}{j[0]:8.3f}{j[1]:8.3f}{j[2]:8.3f}\n')
            
        out.write(f'   {L}\t{L}\t{L}\n')
    return ()
           

gro(frames[0],N)

import mdtraj as md

# Carga el archivo .gro en un objeto Trajectory
traj_gro = md.load('gas.gro')

# Guarda el archivo de topología en formato PDB
traj_gro.save_pdb('gas.pdb')

# Función para guardar las coordenadas en un archivo XTC
def guardar_coordenadas_xtc(archivo_pdb, archivo_xtc, coordenadas_montecarlo):
    # Carga la topología desde el archivo PDB
    top = md.load_pdb(archivo_pdb).topology

    # Crea un objeto Trajectory con las coordenadas y la topología
    #for i in coordenadas_montecarlo:
    #    traj = md.Trajectory(i, top)
    traj = md.Trajectory(coordenadas_montecarlo, top)
    traj.save_xtc(archivo_xtc)
    # Guarda el archivo XTC
   # traj.save_xtc(archivo_xtc)




guardar_coordenadas_xtc('gas.pdb', 'gas.xtc', frames)






plt.plot(E)
plt.xlabel("Paso")
plt.ylabel("E")
