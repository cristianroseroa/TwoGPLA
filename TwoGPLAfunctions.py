import matplotlib.pyplot as plt
import scipy as sp
import numpy
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
import sys  # Librería de la función de salida
from scipy.fftpack import fft2  # FFT dos dimensional
from scipy.fftpack import fftshift as shift  # Corrimiento al cero
import random

#EN ESTE CÓDIGO SE ALMACENAN TODAS LAS FUNCIONES UTILIZADAS Y
# NO UTILIZADAS PERO QUE SE PUEDEN UTILIZAR EN LOS OTROS CÓDIGOS
def Maxwell(v,vd,vt):
    Distribucion = sp.exp(-(v - vd) * (v - vd) / (vt*vt))
    return Distribucion
def CargarRandomPosicion(Np, plasma_lix, plasma_liy, Lx,Ly):
    pst= sp.zeros((Np,2)) #pst are the positions x and y of the system
    for i in range(Np):
        pst[i, 0] = random.uniform(plasma_lix, Lx)
        pst[i, 1] = random.uniform(plasma_liy, Ly)
    return pst
def CargarVelocidadMaxwell2D(Np, n, vd, vt ):
    print("Configuracion inicial de la distribución de velocidad")
    v = sp.zeros((Np, 2))
    C= n/((vt*vt*vt)*(sp.pi**(3/2))) #Constante de Normalizacion
    fmax = C
    vmin1 = vd - 5.0 * vt
    vmax1 = vd + 5.0 * vt
    vmin2 = - 5.0 * vt
    vmax2 = + 5.0 * vt
    for i in range(Np):
        while True:
            vtemp1 = vmin1 + (vmax1 - vmin1) * (sp.random.random())
            vtemp2 = vmin2 + (vmax2 - vmin2) * (sp.random.random())
            f = C * Maxwell(vtemp1, vd, vt) * Maxwell(vtemp2, 0, vt)
            xtemp1 = fmax * (sp.random.random())
            if xtemp1 < f:
                break
        v[i, 0] = vtemp1
        v[i, 1] = vtemp2
    return v
def CondicionesDeFrontera(Np,pst, Lx, Ly):
    #  Periodicos
    # pst son las posiciones
    # Lx , Ly longitudes de malla
    # Np  número de partículas
    for i in range(Np):
        if (pst[i,0] < 0.0):
            while(pst[i,0] < 0.0):
                  pst[i,0] += Lx
        elif (pst[i,0] >= Lx):
            while(pst[i,0] >=Lx):
                  pst[i,0] -= Lx
        if (pst[i,1] < 0.0):
            while (pst[i,1] < 0.0):
                pst[i,1] += Ly
        elif (pst[i,1] >= Ly):
            while (pst[i,1] >= Ly):
                pst[i,1] -= Ly
    return pst
def DensidadCarga( positions, Np, q, Np_Mallax, Np_Mallay,dx,dy,):
    rho = sp.zeros((Np_Mallax + 1, Np_Mallay + 1))  # Densidad electrones
    for i in range(Np):
        re = q[i]
        xa = (positions[i, 0] / dx)
        ya = (positions[i, 1] / dy)
        j1 = int(xa)
        j2 = int(ya)
        h1 = xa - j1
        h2 = ya - j2
        f1 = (1 - h1) * (1 - h2)
        f2 = h1 * (1 - h2)
        f3 = (1 - h1) * h2
        f4 = h1 * h2
        rho[j1][j2] += re * f1
        rho[j1 + 1][j2] += re * f2
        rho[j1][j2 + 1] += re * f3
        rho[j1 + 1][j2 + 1] += re * f4

    for i in range(1,Np_Mallax-1):
        if i==63:
            print(i)
        rho[i][0] += rho[i][Np_Mallay]
        rho[i][Np_Mallay] = rho[i][0]
    for j in range(1,Np_Mallay-1):
        rho[0][j] += rho[Np_Mallax][j]
        rho[Np_Mallax][j] = rho[0][j]
    rho[0][0]+= rho[0][Np_Mallay]+ rho[Np_Mallax][0] + rho[Np_Mallax][Np_Mallay]
    rho[0][Np_Mallay] = rho[0][0]
    rho[Np_Mallax][0] = rho[0][0]
    rho[Np_Mallax][Np_Mallay] = rho[0][0]

    # rho /= (dx * dy * dx * dy)

    return rho
def potential(Nx, Ny, dx, dy, rho):
    rho_k = numpy.fft.fftn(rho)
    Wx = numpy.exp(2 * 1j * numpy.pi / Nx)
    Wy = numpy.exp(2 * 1j * numpy.pi / Ny)
    Wn = 1.0 + 0.0j
    Wm = 1.0 + 0.0j
    dx_2 = dx * dx
    dy_2 = dy * dy

    for n in range(Nx):
        for m in range(Ny):
            denom = dy_2 * (2 - Wn - 1.0/Wn) + dx_2 * (2 - Wm - 1.0/Wm)
            if denom != 0:
                rho_k[n, m] *= dx_2 * dy_2 / denom
            Wm *= Wy
        Wn *= Wx

    phi = numpy.fft.ifftn(rho_k)
    phi = numpy.real(phi)

    return phi
def Campo(Nx, Ny, dx, dy, phi):
    E = numpy.zeros((Nx +1, Ny +1, 2))
    for j in range( Ny+1):
        for i in range(Nx+1):
            nxt_i= (i + 1)
            prv_i = (i - 1)
            if prv_i == -1:
                prv_i = Nx - 1
            if nxt_i == (Nx + 1):
                nxt_i = 1
            E[i, j, 0] = (phi[prv_i, j] - phi[nxt_i, j]) / (2* dx)

    for i in range(Nx +1):
        for j in range(Ny + 1):
            nxt_j = (j + 1)
            prv_j = (j - 1)
            if prv_j == -1:
                prv_j = Ny - 1
            if nxt_j == (Ny + 1):
                nxt_j = 1
            E[i, j, 1] = (phi[i, prv_j] - phi[i, nxt_j]) / (2 * dy)

    return E
def MovimientoParticulas(Np,E, pst, v, dx, dt, dy, carga_masa,step):
    """Se interpola el campo Ex desde la malla a la particula"""
    E /=(dx*dy)
    for i in range(Np):
        # Cloud In Cell
        # -----------------------------------
        xa = (pst[i,0] / dx)
        ya = (pst[i,1] / dy)
        j1 = int(xa)
        j2 = int(ya)
        h1 = xa - j1
        h2 = ya - j2
        f1 = (1 - h1) * (1 - h2)
        f2 = h1 * (1 - h2)
        f3 = (1 - h1) * h2
        f4 = h1 * h2
        # Repartir valores de campo al numero de particulas ingresado
        exi = f1 * E[j1][j2][0] + f2 * E[j1+1][j2][0] + f3 * E[j1][j2+1][0] + f4 * E[j1+1][j2+1][0]
        eyi = f1 * E[j1][j2][1] + f2 * E[j1+1][j2][1] + f3 * E[j1][j2+1][1] + f4 * E[j1+1][j2+1][1]
        # -----------------------------------
        # Actualizar velocidades/Fuerza de Lorentz en diferencias
        if step == 0:
             v[i, 0] -= carga_masa[i] * (dt/2) * exi
             v[i, 1] -= carga_masa[i] * (dt /2)* eyi
        v[i, 0] += carga_masa[i] * dt * exi
        v[i, 1] += carga_masa[i] * dt * eyi
        pst[i, 0] += dt * v[i,0]
        pst[i, 1] += dt * v[i,1]


    return pst, v
def FuncionDV(v, dvx, dvy, Npvelx, Npvely, v0x, v0y, Lx, Ly):
    Np = len(v)
    Fz = sp.zeros((Npvelx, Npvely))
    for i in range(Np):
        gx = (v[i,0] + v0x) / dvx
        gy = (v[i,1] + v0y) / dvy
        jx = int(gx)
        jy = int(gy)
        hx = gx - jx
        hy = gy - jy
        f1 = (1 - hx) * (1 - hy)
        f2 = hx * (1 - hy)
        f3 = (1 - hx) * hy
        f4 = hx * hy
        Fz[jx][jy] = Fz[jx][jy] + f1
        Fz[jx + 1][jy] = Fz[jx + 1][jy] + f2
        Fz[jx][jy + 1] = Fz[jx][jy + 1] + f3
        Fz[jx + 1][jy + 1] = Fz[jx + 1][jy + 1] + f4
        Fz /= Lx*Ly
    return Fz
def Grafica(x,y,vx,vy,Np):
    x1 = sp.zeros(Np//2)
    x2 = sp.zeros(Np//2)
    vx1 = sp.zeros(Np//2)
    vx2 = sp.zeros(Np//2)
    y1 = sp.zeros(Np//2)
    y2 = sp.zeros(Np//2)
    vy1 = sp.zeros(Np//2)
    vy2 = sp.zeros(Np//2)
    x1, x2 = sp.split(x, 2)
    vx1, vx2 = sp.split(vx, 2)
    [y1, y2] = sp.split(y, 2)
    [vy1, vy2] = sp.split(vy, 2)
    plt.figure('Malla')
    plt.plot(x, y, '.', markersize=1)


    plt.figure('Velocidades x')
    plt.plot(x1, vx1, '.', markersize=1, c='m')
    plt.plot(x2, vx2, '.', markersize=1, c='c')



    plt.figure('Velocidades y')
    plt.plot(y1, vy1, '.', markersize=1, c='m')
    plt.plot(y2, vy2, '.', markersize=1, c='c')


    # vex2, vey2 = sp.meshgrid(vex, vey)
    # fig = plt.figure('Función de distribución')
    # ax1 = fig.add_subplot(111,projection='3d')
    # ax1.plot_surface(vex2, vey2, Fz)

    # X, Y = sp.meshgrid(vector_mallay, vector_mallax)
    # fig = plt.figure('Densidad carga')
    # ax2 = fig.add_subplot(111, projection='3d')
    # ax2.plot_surface(X, Y, rhot)
def GraficaRho(vector_mallax,vector_mallay,rhot):
     X, Y = sp.meshgrid(vector_mallay, vector_mallax)
     fig = plt.figure('Densidad carga')
     ax2 = fig.add_subplot(111, projection='3d')
     ax2.plot_surface(X, Y, rhot)
def GraficaFD(Fz,vex,vey):
     vex2, vey2 = sp.meshgrid(vex, vey)
     fig = plt.figure('Función de distribución')
     ax1 = fig.add_subplot(111,projection='3d')
     ax1.plot_surface(vex2, vey2, Fz)
def EspacioFase2D(pst, v, dimension, Ni, Npb1, Npb2):
    N1 = Ni
    N2 = Ni + Npb1
    N3 = Ni + Npb1 + Npb2
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    # z_line = numpy.linspace(0,15,1000)
    # x_line = numpy.cos(z_line)
    # y_line = numpy.sin(z_line)
    # ax.plot3D(x_line, y_line, z_line, 'gray')
    z1_points = v[0: N1, dimension]
    x1_points = pst[0: N1, 0]
    y1_points = pst[0: N1, 1]
    ax.scatter3D(x1_points, y1_points, z1_points, c = 'b', cmap='hsv')
    z2_points = v[N1: N2, dimension]
    x2_points = pst[N1: N2, 0]
    y2_points = pst[N1: N2, 1]
    ax.scatter3D(x2_points, y2_points, z2_points, c = 'c', cmap='hsv')
    z3_points = v[N2: N3, dimension]
    x3_points = pst[N2: N3, 0]
    y3_points = pst[N2: N3, 1]
    ax.scatter3D(x3_points, y3_points, z3_points, c = 'm', cmap='hsv')
    plt.show()
def energias( phi, rho, Nx, Ny, pst, v, m):

     # Energía cinética
     np = len(pst[:,0])
     KE = 0
     PE = 0
     for p in range(np):
         v2 = numpy.dot(v[p, :], v[p, :]) # norma de v al cuadrado.
         KE += m[p] * v2
     KE = KE * 0.5
     # Energía potencial
     for i in range(Nx):
         for j in range(Ny):
              PE += rho[i, j] * phi[i, j]
     PE = PE * 0.5
     return KE, PE