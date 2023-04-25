# A revisar, para poder implementarlo
# verificar que las funciones mencionasa existan...
#
import numpy as np

# Obtener la matriz de posiciones de las partículas
# Cada fila es un tiempo, cada columna es una partícula
# Las posiciones están en unidades arbitrarias
positions = np.loadtxt('trayectorias.txt')

# Calcular la distancia entre todas las posibles combinaciones de pares de partículas
distances = []
N = positions.shape[1] # número de partículas
for i in range(N):
    for j in range(i+1, N):
        r_ij = np.sqrt(np.sum((positions[:,i] - positions[:,j])**2))
        distances.append(r_ij)

# Calcular la distancia media
r_mean = np.mean(distances)

print('La distancia media entre las partículas es:', r_mean)


##### alternativa..
import MDAnalysis as mda

# Cargar el archivo de trayectorias (.xtc)
u = mda.Universe('topology.pdb', 'trayectorias.xtc')

# Seleccionar todas las partículas
particulas = u.select_atoms('all')

# Calcular la distancia media
r_mean = 0.0
N = len(particulas)
for ts in u.trajectory:
    r_ij = particulas.positions[:, None, :] - particulas.positions[None, :, :]
    r_ij = np.sqrt(np.sum(r_ij**2, axis=-1))
    np.fill_diagonal(r_ij, np.inf) # Excluir la comparación de cada partícula consigo misma
    r_mean += np.sum(r_ij) / (N*(N-1))

r_mean /= len(u.trajectory)

print('La distancia media entre las partículas es:', r_mean)

###### Usando MDtraj
import mdtraj as md

# Cargar el archivo de trayectorias (.xtc)
traj = md.load('trayectorias.xtc', top='topology.pdb')

# Seleccionar todas las partículas
particulas = traj.topology.select('all')

# Calcular la distancia media
distancias = md.compute_distances(traj, particulas)
r_mean = np.mean(distancias)

print('La distancia media entre las partículas es:', r_mean)


