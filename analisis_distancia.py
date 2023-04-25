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
