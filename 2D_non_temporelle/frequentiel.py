import numpy as np
import matplotlib.pyplot as plt
from math import cos
import math
import cmath
import matplotlib.pyplot as plt2
import random
from matplotlib import cm
from matplotlib.animation import FuncAnimation

n = 30
Lx=0.003*8
Ly = 0.003*8
Lz=0.003*8
t = 10 #temps limite
dx = Lx / (n - 1)
dy = Ly / (n - 1)
H = dx  # ou H = dy, selon ce que vous souhaitez
w = 2*3.14/5
D = 0.0000178
W=0.0000593 #Perméabilité de la membrane
dt=0.1
#dt = dx**2 / (D / 4)
gamma = D/W
alpha = (complex(0.1)*w+4*dt*D/H)

def bijection(i, j):
    return i + n * j
def reciproque(k):
    i = k % n
    j = (k - i) / n
    return [i, j]

Pp=150
Pb=38.5
Ca=5.8*0.01*Pp
Cb=5.8*0.01*Pb
C1=5.8*0.01*(20.9-17)*Pp/(100*2)
C0=Ca-Cb-C1

A= np.zeros((n * n, n * n))
g = np.zeros(n * n)
b1 = np.zeros((n, n))
for i in range(0,n):
    k = bijection(i, n - 1)
    A[k][k] = 1
    g[k] = C0
for i in range(1, n - 1):
    A[i][i] =4-(1+H/gamma)
    A[i][i + 1] = -1
    A[i][i - 1] = -1
    A[i][i + n] = -1
for i in range(1, n-1):
    for j in range(1, n-1):
        k = bijection(i, j)
        A[k][k] =4
        A[k][k + 1] = -1
        A[k][k - 1] = -1
        A[k][k - n] = -1
        A[k][k + n] = -1
for j in range(1, n - 1):
    k = bijection(0, j)
    A[k][k] = 3
    A[k][k + 1] = -1
    A[k][k + n] = -1
    A[k][k - n] = -1
for j in range(1, n - 1):
    k = bijection(n - 1, j)
    A[k][k] =3
    A[k][k - 1] = -1
    A[k][k + n] = -1
    A[k][k - n] = -1
f1 = bijection(0, 0)
A[f1][f1] =2-H/gamma
A[f1][f1 + 1] =-1
A[f1][f1 + n] =-1
g1 = bijection(n - 1, 0)
A[g1][g1] = 2-H/gamma
A[g1][g1 - 1] =-1
A[g1][g1 + n] = -1
R= np.dot(np.linalg.inv(A), g)
for k in range(0, n * n):
    l = reciproque(k)
    i = l[0]
    j = int(l[1])
    b1[j, i] = R[k]
nt=int(t/dt)
T=np.linspace(0,t,nt)
X=np.zeros((nt,n*n),dtype=np.complex128)
for m in range(nt):
    c = np.zeros((n * n, n * n),dtype=np.complex64)
    h = np.zeros(n * n)
    X0 = C1 * cos(w * m)
    for i in range(0,n):
        k = bijection(i, n - 1)
        c[k][k] = 1
        h[k] = X0
    for i in range(1, n - 1):
        c[i][i] = alpha - (dt * D * (1 + H / gamma)) / H
        c[i][i + 1] = -dt * D / H
        c[i][i - 1] = -dt * D / H
        c[i][i + n] = -dt * D / H
    for i in range(1, n-1):
        for j in range(1, n-1):
            k = bijection(i, j)
            c[k][k] = alpha
            c[k][k + 1] = -dt * D / H
            c[k][k - 1] = -dt * D / H
            c[k][k - n] = -dt * D / H
            c[k][k + n] = -dt * D / H
    for j in range(1, n - 1):
        k = bijection(0, j)
        c[k][k] = alpha - dt * D / H
        c[k][k + 1] = -dt * D / H
        c[k][k + n] = -dt * D / H
        c[k][k - n] = -dt * D / H
    for j in range(1, n - 1):
        k = bijection(n - 1, j)
        c[k][k] = alpha - dt * D / H
        c[k][k - 1] = -dt * D / H
        c[k][k + n] = -dt * D / H
        c[k][k - n] = -dt * D / H
    f = bijection(0, 0)
    c[f][f] = alpha - (dt * D * (2 + H / gamma)) / H
    c[f][f + 1] = -dt * D / H
    c[f][f + n] = -dt * D / H
    v = bijection(n - 1, 0)
    c[v][v] = alpha - (dt * D * (2+ H / gamma)) / H
    c[v][v- 1] = -dt * D / H
    c[v][v + n] = -dt * D / H
    M = np.dot(np.linalg.inv(c), h)
    for k in range(n*n):
        X[m][k]=np.real(M[k]*cmath.exp(complex(0,1)*w*m)+R[k])
    #X[]=np.real(M*cmath.exp(complex(0,1)*w*i))

b = np.zeros((nt,n,n))
for i in range(nt):
    for k in range(0, n * n):
        l = reciproque(k)
        d = l[0]
        j = int(l[1])
        b[i,j,d] = X[i][k]
x = np.linspace(0, Lx, n)
y = np.linspace(0, Ly, n)
Z, Y = np.meshgrid(x, y)
fig, ax = plt.subplots()
po = b[0]
pcm = ax.pcolormesh(Z, Y, po, cmap=cm.jet,vmin=po.min(),vmax=po.max())
ax.set_title('Temporelle')
# Add colorbar
cbar = plt.colorbar(pcm, ax=ax)
flux_sortie=0
flux_entree=0
lentry=np.zeros(nt)
lsortie=np.zeros(nt)
def update(num):
    ax.clear()
    po = b[num]
    pcm=ax.pcolormesh(Z, Y, po, cmap=cm.jet,vmin=po.min(),vmax=po.max())
    ax.set_title('Temporelle')
    for i in range (n):
        k=bijection(i,n-1)
        flux_entrée=X[num][k]-X[num][k-n]
    flux_entrée=-D*flux_entrée*(dx/dy)*Lz
    for i in range (n):
         k=bijection(i,0)
         flux_sortie=X[num][k+n]-X[num][k]
    flux_sortie=-D*flux_sortie*(dx/dy)*Lz
    lentry[num]=flux_entrée*30000*0.001
    lsortie[num]=flux_sortie*30000*0.001

    #cbar = plt.colorbar(pcm, ax=ax)
    #cbar.clear()
ani = FuncAnimation(fig, update, frames=nt, interval=1000)
#plt.subplot(2, 1, 2)
# Plotting the results
#plt.pcolormesh(Z, Y, z, cmap='Viridis')  # You can also use plt.pcolormesh for a different style
'''vmin=10
vmax=C_0
step=0.1'''
'''mesh=plt.pcolormesh(Z, Y, z, cmap=cm.jet, vmin=0, vmax=C_0)  # You can also use plt.pcolormesh for a different style
colorbar=plt.colorbar(mesh)'''  # Add a colorbar for reference
plt.title('Solution en concentration')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.tight_layout()
plt.show()
pob=np.linspace(0,t,nt)
plt.subplot(2, 1, 1)
plt.plot(pob,lentry)

# Add labels and title
plt.xlabel('Temps en second')
plt.ylabel('Débit entrée en mol/s')
plt.title('Débit en fonction du temps')
# Display the plot
#plt.show()
plt.subplot(2, 1, 2)
plt.plot(pob,lsortie)

# Add labels and title
plt.xlabel('Temps en second')
plt.ylabel('Débit sortie en mol/s')
#plt.title('Debit en fonction du temps')
# Display the plot
plt.show()