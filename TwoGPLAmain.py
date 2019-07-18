import matplotlib.pyplot as plt
import numpy
import scipy as sp
# from tqdm import tqdm
import TwoGPLAfunctions as TD
from itertools import product


# DINÁMICA = 1 grafica la evolución del espacio de fase Vx vs X
DINAMICA = 1
Npb1, Npb2, Ni,Nx, Ny, Lx, Ly, dx, dy,steps, dt = numpy.load('Semilla/semilladatos.npy')
positions1 = numpy.load("Semilla/semillaposicion1.npy")
velocities1 = numpy.load("Semilla/semillavelocidad1.npy")
positions2 = numpy.load("Semilla/semillaposicion2.npy")
velocities2 = numpy.load("Semilla/semillavelocidad2.npy")
positionsi = numpy.load("Semilla/semillaposicioni.npy")
velocitiesi = numpy.load("Semilla/semillavelocidadi.npy")
Npb1, Npb2, Ni = int(Npb1), int(Npb2), int(Ni)
steps = int(steps)
Nx, Ny = int(Nx), int(Ny)
NP = Npb1 + Npb2 + Ni
positions=sp.zeros((NP, 2))
velocities=sp.zeros((NP, 2))
for i in range(NP):
    if i<Ni:
        positions[i,0]=positionsi[i,0]
        positions[i,1]=positionsi[i,1]
        velocities[i,0]=velocitiesi[i,0]
        velocities[i,1]=velocitiesi[i,1]
    if i>=Ni:
        if i<(Ni+Npb1):
             positions[i,0]=positions1[i-Ni,0]
             positions[i,1]=positions1[i-Ni,1]
             velocities[i,0]=velocities1[i-Ni,0]
             velocities[i,1]=velocities1[i-Ni,1]
    if i >=(Ni+Npb1):
        positions[i, 0] = positions2[i -Ni-Npb1, 0]
        positions[i, 1] = positions2[i -Ni-Npb1, 1]
        velocities[i, 0] = velocities2[i -Ni-Npb1, 0]
        velocities[i, 1] = velocities2[i -Ni-Npb1, 1]
NP = len(positions)
EPS0 = 1
QoverM = numpy.append((1/100)*numpy.ones(Ni), -numpy.ones(Npb1 + Npb2))
QoverM1 = numpy.append(numpy.ones(Ni), -numpy.ones(Npb1 + Npb2))
moves = numpy.append(numpy.ones(Ni), numpy.ones(Npb1+ Npb2))
charges = EPS0* Lx * Ly * QoverM / (Npb2+Npb1)
masses = charges / QoverM
E = sp.zeros((Nx + 1, Ny + 1, 2, steps + 1))
phi = sp.zeros((Nx + 1, Ny + 1, steps + 1))
rho = sp.zeros((Nx + 1, Ny + 1, steps + 1))
print("Simulation running...\n")
for step in range(0, steps):
    print(step)
    positions = TD.CondicionesDeFrontera(NP, positions, Lx, Ly)
    if step % 10 == 0 and DINAMICA == 1:
        x1 = positions[Ni:(Ni+Npb1), 0]
        y1 = velocities[Ni:(Ni+Npb1), 0]
        x2 = positions[(Ni+Npb1):, 0]
        y2 = velocities[(Ni+Npb1):, 0]
        v = velocities[Ni:, 0]
        plt.scatter(x1, y1, s=0.5)
        plt.scatter(x2, y2, s=0.5)
        plt.pause(0.01)
        plt.clf()
    rho[:, :, step] = TD.DensidadCarga(positions, NP, charges, Nx, Ny, dx, dy)
    phi[:, :, step] = TD.potential(Nx+1, Ny+1, dx, dy, rho[:, :, step])
    E[:, :, :, step] = TD.Campo(Nx , Ny, dx, dy, phi[:, :, step])
    positions, velocities = TD.MovimientoParticulas(NP, E[:, :, :, step], positions, velocities, dx, dt, dy, QoverM, step)
    EnergiaK, EnergiaP = TD.energias(phi[:, :, step], rho[:, :, step], Nx, Ny, positions, velocities, masses)
    #Se guardan las magnitudes importantes
    if step % 20 == 0:
         numpy.save('Posiciones/positions{0}.npy'.format(step), positions)
         numpy.save('Velocidades/velocities{0}.npy'.format(step), velocities)
         numpy.save('Densidad/rhoacumulador{0}.npy'.format(step), rho)
         numpy.save('Potencial/phiacumulador{0}.npy'.format(step), phi)
         numpy.save('Campo/Eacumulador{0}.npy'.format(step), E)
         numpy.save('Energías/EnergiaK{0}.npy'.format(step), EnergiaK)
         numpy.save('Energías/EnergiaP{0}.npy'.format(step), EnergiaP)
print("\nSimulation finished!\n")




