import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


### la bijection et sa reciproque qui transforme carte <----------> vecteur colonne
def bij2D(i, j, n):
    return (i + (j - 1) * n)


def bij2D_1(k, n):
    j = k // n + 1
    i = k - (j - 1) * n
    return (i, j)


### Initialisation (valeurs aleatoires)
n = 10  # ligne
m = n  # colonne
M = np.eye(n * m)
D = 7 * 10 ** (-11)
Lx = 200 * 10 ** (-9)  ### profondeur de la couche
Nx = n  ## nbr des points
dx = Lx / (Nx - 1)
dt = 4.85e-07
alpha = dt * D / (2 * (dx ** 2))
beta = 2 * (1 + 2 * alpha)
gamma = 2 * (1 - 2 * alpha)
point_inter = []
### boucle qui donne la liste K des pts qui ont 4 voisins
for j in range(2, n):
    for i in range(2, m):
        point_inter.append(bij2D(i, j, n))
### boucle qui genere la matrice
for k in point_inter:
    M[k - 1, k - 1] = beta
    M[k - 1, k] = -alpha
    M[k - 1, k - 2] = -alpha
    M[k - 1, k + n - 1] = -alpha
    M[k - 1, k - n - 1] = -alpha
for k in range(n):
    M[k][k + n] = -1
for k in range(1, n):
    M[k * n][k * n + 1] = -1
    M[n - 1 + k * n][n - 2 + k * n] = -1


### implementation des CL

### la bijection ksi entre [1,n]*[1,n] --------->[1,n]*[1,n] qui verifie M*X(t+1)=ksi(X0)

def ksi(n, X):
    Ximg = np.zeros((n * n, 1))
    for k in point_inter:
        Ximg[k - 1] = alpha * (X[k - 2] + X[k] + X[k + n - 1] + X[k - n - 1]) + gamma * X[k - 1]
    return (Ximg)


X0 = np.zeros((n * n, 1))
for i in range(n):
    X0[i] = 1
col = X0
"Partie de simulation"
"initiation de la colonne du vecteur des conditions initiales"
iM = np.linalg.inv(M)
"creation des indices"
x = np.linspace(0, n - 1, n, endpoint=True).astype(int)
y = np.linspace(0, n - 1, n, endpoint=True).astype(int)
X, Y = np.meshgrid(x, y)
"creation de la figure de simulation"
fig, ax = plt.subplots()
"calcul avec la fonction ksi du premier vecteur et le rendre une matrice de taille n*n"
col = ksi(n, col)
col = np.dot(iM, col).reshape(n, n)
"affiliation de chacune des valeurs de la matrice Ã  Z"


def function_z(i):
    dif = col[i][i]
    return dif


Z = function_z(x)
"plot of values accordingly to a color"
pcm = ax.pcolormesh(X, Y, Z, cmap='plasma', shading='gouraud', vmin=0, vmax=0.005)
fig.colorbar(pcm, shrink=0.6)
plt.title("1")
# plt.pause(8)
"rendre la matrice de taille n*n Ã  un vecteur "
col = np.reshape(col, (n * n, 1))
"cette fonction fait la meme chose que la partie d'initiation mais en la repetant N fois avec dt difference de temps entre la figure et celle qui la suit "


def update(i):
    global col
    col = ksi(n, col)
    col = np.dot(iM, col).reshape(n, n)

    def function_z(i):
        dif = col[i][i]
        return dif

    Z = function_z(x)
    ax.pcolormesh(X, Y, Z, cmap='plasma', shading='gouraud', vmin=0, vmax=0.005)
    plt.title(str(i + 1))
    col = np.reshape(col, (n * n, 1))


# partie de creation de video
ani = animation.FuncAnimation(fig, update, frames=500, interval=1, repeat=False, blit=False)
plt.show()
