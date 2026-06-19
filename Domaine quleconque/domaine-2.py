import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib import cm

# 0 Points non actifs
# 1 Points actifs sans condition limite
# 2 Condition en entrée
# 3 Condition de Robin
# 4 Condition en Neumann
# 5 Combinaison de Robin et Neumann
# 6 Combinaison de deux conditions Neumann suivant deux axes différents


# Les paramètres d'activation
n=101 # Nombre de points impair
Lx=10 # Longeur du domaine x
Ly=10 # Longeur du domaine y
m=2 # Nombre de bronches ou membranes
a=int(n/2)
b=int((2/10) * n)
c=int((3/10) * n)
d=int((2/10) * n)
e=int((2/10) * n)


# Les paramètres de concentration
D=0.0000178 #Diffusivité de O2 dans l'air
W=0.0000593 #Perméabilité de la membrane
Pp=150 #Pression partielle en mmHg
#C_0=5.8*(10e-6)*Pp #Concentration en O2 en mol/m3
C_0=160
dx=(Lx/(n-1))
dy=(Ly/(n-1))
lamda=D/W
alpha=(Lx/Ly)

# Les fonctions de bijection directe et reciproque entre (i,j) et k
def bijection(i,j):
    return(i+(j*n))

def reciproque(k):
    i=k%n
    j=int((k-i)/n)
    return[i,j]

A=np.ones((n*n)) # Matrice d'activation

if (m==2):
    for i in range(0,b):
        for j in range(0,c+1):
            k1=bijection(a-i,j)
            k2=bijection(a+i,j)
            A[k1]=0
            A[k2]=0
else:
    l=[]
    for i in range(m):
        l.append([int((i+1/2)*(n/m)),random.randint(1,int(n/(2*m))-3),random.randint(1,c+1)])
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
        A[k2]=3

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
# Plot en rouge (0) et en bleu (1)
plt.scatter(i, j, c=A, cmap='plasma', marker='o', s=50)
# Définir les étiquettes et le titre
plt.xlabel('Indice i')
plt.ylabel('Indice j')
plt.title('Plot de Maillage')
# Afficher la barre de couleur
plt.colorbar()
# Afficher le plot
plt.show()



C=np.zeros((n*n,n*n))
B=np.zeros(n*n)

# Concentration des points 5 (Combinaison de Robin et Neumann)
f=bijection(0,0)
C[f][f]=-(1+(2*alpha)-(alpha*(1+(dy/lamda))))
C[f][f+1]=1
C[f][f+n]=alpha

g=bijection(n-1,0)
C[g][g]=-(1+(2*alpha)-(alpha*(1+(dy/lamda))))
C[g][g-1]=1
C[g][g+n]=alpha

# Concentration des points 4 (Condition en Neumann)
for k in range (n*n):
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
   
    if A[k]==3:
        for i in range(n):
            if k==bijection(i,0):
                C[k][k]=-(2+alpha*(1-(dy/lamda)))
                C[k][k+1]=1
                C[k][k-1]=1
                C[k][k+n]=alpha
            else:
                if A[k+1]==0:
                    C[k][k]=-(1+2*alpha+(dx/lamda))
                    C[k][k-n]=alpha
                    C[k][k-1]=1
                    C[k][k+n]=alpha 
                elif A[k-1]==0:
                    C[k][k]=-(1+2*alpha-(dx/lamda))
                    C[k][k-n]=alpha
                    C[k][k+1]=1
                    C[k][k+n]=alpha
                elif A[k-n]==0:
                    C[k][k]=-(2+alpha*(1-(dy/lamda)))
                    C[k][k+1]=1
                    C[k][k-1]=1
                    C[k][k+n]=alpha
    if A[k]==6:
        C[k][k]=1
    if A[k]==2:
        C[k][k]=1
        B[k]=C_0
    if A[k]==1:
        C[k][k]=-2*(1+alpha)
        C[k][k+1]=1
        C[k][k-1]=1
        C[k][k+n]=alpha
        C[k][k-n]=alpha
    if A[k]==0:
        C[k][k]=1

       
z=np.zeros((n,n))
s=np.zeros((n*n))
for i in range(n*n):
    if np.array_equal(s, C[i]):
        print(i)
X=np.linalg.solve(C,B)

for i in range(n*n):
    print(X[i])

for k in range(0,n*n):
    l=reciproque(k)
    i=l[0]
    j=int(l[1])
    z[j,i]=X[k]
x = np.linspace(0, Lx, n)
y = np.linspace(0, Ly, n)
Z, Y = np.meshgrid(x, y)



# Plotting the results
#plt.pcolormesh(Z, Y, z, cmap='Viridis')  # You can also use plt.pcolormesh for a different style
plt.contourf(Z, Y, z, cmap='plasma', vmin=-50, vmax=160)  # You can also use plt.pcolormesh for a different style
plt.colorbar()  # Add a colorbar for reference
plt.title('Solution en concentration')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()