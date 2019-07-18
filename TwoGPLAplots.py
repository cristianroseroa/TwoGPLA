import matplotlib.pyplot as plt
import scipy as sp
import numpy
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
import sys  # Librería de la función de salida
from scipy.fftpack import fft2  # FFT dos dimensional
from scipy.fftpack import fftshift as shift  # Corrimiento al cero
import TwoGPLAfunctions2 as TD



# En este código solo se grafican de momento potencial y espacio de fase
Npb1, Npb2, Ni,Nx, Ny, Lx, Ly, dx, dy, steps, dt = numpy.load('Semilla/semilladatos.npy')
totalsteps=2980
step = 30
vector_mallax = sp.arange(0, Lx + dx, dx)
vector_mallay = sp.arange(0, Ly + dy, dy)
vector_tiempo = dt * sp.arange(0, totalsteps)
Npvelx = 100  # Numero de puntos de velocidad X
Npvely = 100  # Numero de puntos de velocidad Y
v0x = 30
v0y = 30
dvx = 2 * v0x / Npvelx
dvy = 2 * v0y / Npvely
vex = sp.arange(-v0x, v0x, dvx)
vey = sp.arange(-v0y, v0y, dvy)
X=numpy.arange(0,Lx+dx,dx)
Y=numpy.arange(0,Ly+dy,dy)

positions1 = numpy.load("Semilla/semillaposicion1.npy")
velocities1 = numpy.load("Semilla/semillavelocidad1.npy")
positions2 = numpy.load("Semilla/semillaposicion2.npy")
velocities2 = numpy.load("Semilla/semillavelocidad2.npy")
positionsi = numpy.load("Semilla/semillaposicioni.npy")
velocitiesi = numpy.load("Semilla/semillavelocidadi.npy")
rho = numpy.load('Densidad/rhoacumulador{0}.npy'.format(totalsteps))
phi = numpy.load('Potencial/phiacumulador{0}.npy'.format(totalsteps))
E = numpy.load('Campo/Eacumulador{0}.npy'.format(totalsteps))
EnergiaK = numpy.load('Energías/EnergiaK{0}.npy'.format(totalsteps))
EnergiaP = numpy.load('Energías/EnergiaP{0}.npy'.format(totalsteps))

Npb1 = int(Npb1)
Npb2 = int(Npb2)
Ni = int(Ni)
NP = Npb1 + Npb2 + Ni
positions=sp.zeros((NP,2))
velocities=sp.zeros((NP,2))
for i in range (NP):
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


for time in range(2500, totalsteps):
    if  time%20==0 :
         x = numpy.load('Posiciones/positions{0}.npy'.format(time))
         v = numpy.load('Velocidades/velocities{0}.npy'.format(time))
         plt.figure('Espacio de fase Vx vs X')
         plt.plot(x[Ni:(Ni + Npb1),0],v[Ni:(Ni + Npb1),0],'.', markersize=0.05, c='darkblue')
         plt.plot(x[(Ni + Npb1):NP,0],v[(Ni + Npb1):NP,0],'.', markersize=0.05, c='r')
#         plt.ylim(-4,6)
#         plt.xlim(0,100)
#         plt.title(u'Espacio de fase',fontsize = 18)
         plt.xlabel('X', fontsize = 14)
         plt.ylabel('Vx', fontsize = 14)
         plt.savefig('EspacioDeFase/EspacioFaseXvsVx{0}.png'.format(time))
         plt.clf()
#         
#         plt.figure('Espacio de fase X vs Vx')
#         plt.plot(positions1[:,0],velocities1[:,0],'.', markersize=0.05, c='m')
#         plt.plot(positions2[:,0],velocities2[:,0],'.', markersize=0.05, c='m')
#         plt.plot(positionsi[:,0],velocitiesi[:,0],'.', markersize=0.05, c='m')
#         plt.ylim(-4,6)
#         plt.xlim(0,100)
##         plt.title(u'Espacio de fase',fontsize = 18)
#         plt.xlabel('X', fontsize = 14)
#         plt.ylabel('Vx', fontsize = 14)
#         plt.savefig('EspacioDeFase/EspacioFaseXvsVx{0}.png'.format(time))
#         plt.clf()
#         plt.figure('Espacio de velocidades Vx vs Vy')
#         plt.plot(v[Ni:(Ni + Npb1),0],v[Ni:(Ni + Npb1),1],'.', markersize=0.05, c='darkblue')
#         plt.plot(v[(Ni + Npb1):NP,0],v[(Ni + Npb1):NP,1],'.', markersize=0.05, c='r')
#         plt.ylim(-0.3,0.3)
#         plt.xlim(-4,6)
##         plt.title(u'Espacio de fase',fontsize = 18)
#         plt.xlabel('Vx', fontsize = 14)
#         plt.ylabel('Vy', fontsize = 14)
#         plt.savefig('EspacioDeFase/EspacioFaseVxvsVy{0}.png'.format(time))
#         plt.clf()
#         plt.figure('Espacio de velocidades Vx vs Vy')
#         plt.plot(velocities1[:,0],velocities1[:,1],'.', markersize=0.05, c='darkblue')
#         plt.plot(velocities2[:,0],velocities2[:,1],'.', markersize=0.05, c='r')
#         plt.ylim(-0.3,0.3)
#         plt.xlim(-4,6)
#         plt.title(u'Espacio de fase',fontsize = 18)
#         plt.xlabel('Vx', fontsize = 14)
#         plt.ylabel('Vy', fontsize = 14)
#         plt.savefig('EspacioDeFase/EspacioFaseVxvsVy{0}.png'.format(time))
#         plt.clf()
         
#         Fz=TD.FuncionDV(v, dvx, dvy, Npvelx, Npvely, v0x, v0y, Lx, Ly)
#         TD.GraficaFD(Fz, vex, vey,time)
         
         plt.figure('Potencial eléctrico')
         plt.contourf(X,Y,phi[:,:,time].T)
#         plt.title(u'Potencial Eléctrico',fontsize = 18)
#         plt.ylim(0,dt*numtime)
#         plt.xlim(0,3)
         clb=plt.colorbar()
         clb.ax.set_title('$\phi$')
         plt.tick_params(labelsize=16)
         plt.savefig('Graficas/Potencial{0}.png'.format(time))
         plt.show()
         plt.draw()
         plt.clf()



#
#TD.EspacioFase2D(positions, velocities, 1, Ni, Npb1, Npb2)

#
# TD.GraficaRho(vector_mallax,vector_mallay,rho)