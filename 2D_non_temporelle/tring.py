import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib import cm

######################### Description de la solution #############
''' On essaye premièrement de construire une matrice notée A contenant des nombres de 0 à 7 suivant
la configuration ci-dessous, cette matrice nous permet essentiellement de déterminer le maillage tout en sachant
l'état de chaque points du maillage.
Deuxièmement, en parcourant la matrice A, en essaye d'affecter les coefficient de l'équation pour chaque cas
dans les matrices C et B, tel que C*X=B  ''' 

# 0 Points non actifs
# 1 Points actifs sans condition limite
# 2 Condition en entrée
# 3 Condition de Robin
# 4 Condition en Neumann
# 5 Combinaison de Robin et Neumann
# 6 Combinaison de deux conditions Neumann suivant deux axes différents
# 7 Condition de Robin en x = 0

# Les paramètres d'activation
n=101 # Nombre de points impair
Lx=0.003 # Longeur du domaine x
Ly=0.003 # Longeur du domaine y
Lz=0.003 # Longeur du domaine z
m=4 # Nombre de bronches ou membranes
a=int(n/2)
b=int((2/10) * n)
c=int((3/10) * n)
d=int((2/10) * n)
e=int((2/10) * n)


# Les paramètres de concentration
D=0.0000178 #Diffusivité de O2 dans l'air
W=0.0000593 #Perméabilité de la membrane
Pp=150 # Pression partielle de O2 à l'entrée du subacinus en mmHg
Pb=35 # Pression partielle de O2 à l'interface de la membrane (sang) en mmHg
C_0=5.8*(0.01)*Pp #Concentration en O2 en mol/m3
C_b=5.8*(0.01)*Pb #Concentration en sang en mol/m3
dx=(Lx/(n-1)) # Pas de discrétisation suivant x
dy=(Ly/(n-1)) # Pas de discrétisation suivant y
lamda=400 # Longueur d'exploration
alpha=(Lx/Ly)

# Les fonctions de bijection directe et reciproque entre (i,j) et k
def bijection(i,j):
    return(i+(j*n))

def reciproque(k):
    i=k%n
    j=int((k-i)/n)
    return[i,j]


A=np.ones((n*n)) # Matrice d'activation


l=[]
for i in range(m):
    l.append([int((i+1/2)*(n/m)),int(n/(4*m))-3,c+1])
for elt in l:
    for k in range(0,elt[1]+1):
        for q in range(0,elt[2]):
            k1=bijection(elt[0]+k,q)
            k2=bijection(elt[0]-k,q)
            A[k1]=0
            A[k2]=0

for i in range(0,a-d):
    for j in range(0,e+1):
        k1=bijection(i,n-j-1)
        k2=bijection(n-i-1,n-j-1)
        A[k1]=0
        A[k2]=0

indices = [reciproque(k) for k in range(n * n)]
for i in range(1,n-1):
    k1=bijection(i,n-1)
    k2=bijection(i,0)
    if A[k1]==1:
        A[k1]=2
    if A[k2]==1:
        A[k2]=7

for j in range(1,n-1):
    k1=bijection(0,j)
    k2=bijection(n-1,j)
    if A[k1]==1:
        A[k1]=4
    if A[k2]==1:
        A[k2]=4

for i in range(n):
    for j in range(1,c+2):
        k=bijection(i,j)
        if (A[k]==1 and (A[k+1]==0 or A[k-1]==0 or A[k-n]==0)):
            A[k]=3

for i in range(n):
    for j in range(c,n):
        k=bijection(i,j)
        if (A[k]==1 and (A[k+1]==0 or A[k-1]==0 or A[k+n]==0)):
            A[k]=4

A[bijection(0,0)]=A[bijection(n-1,0)]=5
A[bijection(0,n-e-2)]=A[bijection(n-1,n-e-2)]=6

# Séparer les coordonnées en listes i et j
i, j = zip(*indices)
plt.subplot(2, 1, 1)
# Plot en rouge (0) et en bleu (1)
plt.scatter(i, j, c=A, cmap='plasma', marker='o', s=50)
# Définir les étiquettes et le titre
plt.xlabel('Indice i')
plt.ylabel('Indice j')
plt.title('Plot de Maillage')
# Afficher la barre de couleur
plt.colorbar()
# Afficher le plot




C=np.zeros((n*n,n*n))
B=np.zeros(n*n)

# Concentration des points 5 (Combinaison de Robin et Neumann)
f=bijection(0,0)
C[f][f]=-(1+(2*alpha)-(alpha*(1-(dy/lamda))))
C[f][f+1]=1
C[f][f+n]=alpha
B[k]=-((dy*C_b)/lamda)
g=bijection(n-1,0)
C[g][g]=-(1+(2*alpha)-(alpha*(1-(dy/lamda))))
C[g][g-1]=1
C[g][g+n]=alpha
B[k]=-((dy*C_b)/lamda)

for k in range (n*n):
    # Concentration des points 4 (Condition en Neumann)
    if A[k]==4:
        if A[k+1]==0 or A[k+1]==4 or A[k+1]==6:
            C[k][k]=-(1+(2*alpha))
            C[k][k-1]=1
            C[k][k+n]=alpha
            C[k][k-n]=alpha
        if A[k-1]==0 or A[k-1]==4 or A[k-1]==5:
            C[k][k]=-(1+(2*alpha))
            C[k][k+1]=1
            C[k][k+n]=alpha
            C[k][k-n]=alpha
        if A[k+n]==0:
            C[k][k]=-(2+alpha)
            C[k][k+1]=1
            C[k][k-1]=1
            C[k][k-n]=alpha
    # Concentration des points 3 (Condition de Robin)
    if A[k]==3:
        if A[k+1]==0:
            C[k][k]=-(1+2*alpha+(dx/lamda))
            C[k][k-n]=alpha
            C[k][k-1]=1
            C[k][k+n]=alpha 
            B[k]=-((dx*C_b)/lamda)
        elif A[k-1]==0:
            C[k][k]=-(1+2*alpha+(dx/lamda))
            C[k][k-n]=alpha
            C[k][k+1]=1
            C[k][k+n]=alpha
            B[k]=-((dx*C_b)/lamda)
        elif A[k-n]==0:
            C[k][k]=-(2+alpha*(1+(dy/lamda)))
            C[k][k+1]=1
            C[k][k-1]=1
            C[k][k+n]=alpha
            B[k]=-((dy*C_b)/lamda)
    # Concentration des points 6 (Combinaison de deux conditions Neumann suivant deux axes différents)
    if A[k]==6:
        C[k][k]=1
    #Concentration des points 2 (Condition de Dirichlet)
    if A[k]==2:
        C[k][k]=1
        B[k]=C_0-C_b
    # Concentration des points 1
    if A[k]==1:
        C[k][k]=-2*(1+alpha)
        C[k][k+1]=1
        C[k][k-1]=1
        C[k][k+n]=alpha
        C[k][k-n]=alpha
    # Concentration des points non actifs (nulle)
    if A[k]==0:
        C[k][k]=1
    # Concentration des points 7 (Points ayant une condition de Robin dans x=0)
    if A[k]==7:
        C[k][k]=-(2+alpha*(1+(dy/lamda)))
        C[k][k+1]=1
        C[k][k-1]=1
        C[k][k+n]=alpha
        B[k]=-((dy*C_b)/lamda)
       
z=np.zeros((n,n))
s=np.zeros((n*n))
for i in range(n*n):
    if np.array_equal(s, C[i]):
        print(i)
X=np.linalg.solve(C,B)

# Calcul du flux
flux_entrée=0

for k in range(n*n):
    if A[k]==2:
        flux_entrée=X[k]-X[k-n]
flux_entrée=D*flux_entrée*(dx/dy)*Lz  
flux_sortie=0
flux_sortie1=0
flux_sortie2=0
flux_sortie3=0
for k in range(n*n):
    if A[k]==7 or (A[k]==3 and A[k-n]==0):
        flux_sortie1=X[k+n]-X[k]
    if A[k]==3 and A[k+1]==0:
        flux_sortie2=X[k]-X[k-1]
    if A[k]==3 and A[k-1]==0:
        flux_sortie3=X[k+1]-X[k]
flux_sortie=D*(-flux_sortie2+flux_sortie3)*(dy/dx)*Lz+D*flux_sortie1*(dx/dy)*Lz

print(flux_entrée)
print("////////")
print(flux_sortie)



for k in range(0,n*n):
    l=reciproque(k)
    i=l[0]
    j=int(l[1])
    z[j,i]=X[k]

x = np.linspace(0, Lx, n)
y = np.linspace(0, Ly, n)
Z, Y = np.meshgrid(x, y)


plt.subplot(2, 1, 2)

vmin=10
vmax=C_0
step=0.1
mesh=plt.pcolormesh(Z, Y, z, cmap='nipy_spectral', vmin=0, vmax=C_0-C_b)  # You can also use plt.pcolormesh for a different style
colorbar=plt.colorbar(mesh)  # Add a colorbar for reference
plt.title('Visualisation 2D stationnaire général')
plt.xlabel('Lx')
plt.ylabel('Ly')
plt.tight_layout()
plt.show()
