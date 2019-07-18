import matplotlib.pyplot as plt
import scipy as sp
import numpy
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
import sys  # Librería de la función de salida
from scipy.fftpack import fft2  # FFT dos dimensional
from scipy.fftpack import fftshift as shift  # Corrimiento al cero
import TwoGPLAfunctions as TD


# Se asignan las variables necesarias como el número de partículas, las dimensiones
# del espacio, los puntos de malla y el paso temporal.
Npb1 = 2000 #Partículas beam 1
Npb2 = 2000 # Partículas beam 2
Npi = 2000 # Partículas beam iones
NpT = Npb1 + Npb2 + Npi
totalsteps=250 #Número total de pasos
dt=0.1 # Paso temporal
# Puntos de malla y tamaño del espacio
Nx, Ny = 64, 10
Lx, Ly = 6, 0.01
dx, dy = Lx/Nx, Ly/Ny
#Densidades (Se usan para calcular las velocidades)
n1 = Npb1/(Lx*Ly)
n2 = Npb2/(Lx*Ly)
ni = Npi /(Lx*Ly)
#Velocidad media de los haces y velocidades térmicas
vd1 = 1
vd2 = -1
vt1 = 0.3
vt2 = 0.3
vti = 1

# Se guardan los datos de que dependen los otros códigos main y plots
semilladatos = Npb1, Npb2, Npi, Nx, Ny, Lx, Ly, dx, dy, totalsteps, dt
#Se llaman funciones para asignar posiciones y velocidades
positions1 = TD.CargarRandomPosicion(Npb1, 0, 0, Lx, Ly)
positions2 = TD.CargarRandomPosicion(Npb2, 0, 0, Lx, Ly)
positionsi = TD.CargarRandomPosicion(Npi, 0, 0, Lx, Ly)
velocities1 = TD.CargarVelocidadMaxwell2D(Npb1, n1, vd1, vt1)
velocities2 = TD.CargarVelocidadMaxwell2D(Npb2, n2, vd2, vt2)
velocitiesi = TD.CargarVelocidadMaxwell2D(Npi, ni, 0, vti)
#Se acomodan todas las partículas en las condiciones de frontera
positions1 = TD.CondicionesDeFrontera(Npb1, positions1, Lx, Ly)
positions2 = TD.CondicionesDeFrontera(Npb2, positions2, Lx, Ly)
positionsi = TD.CondicionesDeFrontera(Npi, positionsi, Lx, Ly)

#Se grafica los puntos para mostrar la configuración
plt.scatter(positions1[:, 0], velocities1[:, 0])
plt.scatter(positions2[:, 0], velocities2[:, 0])
plt.scatter(positionsi[:, 0], velocitiesi[:, 0])
plt.pause(5)
plt.clf()
#Se guardan las posiciones, velocidades y variables en
#la semilla
numpy.save('Semilla/semillavelocidad1.npy', velocities1)
numpy.save('Semilla/semillaposicion2.npy', positions2)
numpy.save('Semilla/semillavelocidad2.npy', velocities2)
numpy.save('Semilla/semillaposicioni.npy', positionsi)
numpy.save('Semilla/semillavelocidadi.npy', velocitiesi)
numpy.save('Semilla/semillaposicion1.npy', positions1)
numpy.save('Semilla/semilladatos.npy', semilladatos)

