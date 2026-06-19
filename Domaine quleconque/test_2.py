import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import random

# 2 Condition en entrée
# 3 Condition de Robin
# 4 Condition en flux nul (Neumann)
# 5 Combinaison de Robin et flux nul (Dirichlet)

def bijection(i,j):
    return(i+(j*n))
def reciproque(k):
    i=k%n
    j=(k-i)/n
    return [i,j]

n=21 # Nombre de points
Lx=10 # Longeur du domaine x
Ly=10 # Longeur du domaine y

Lx=100
Ly=100

D=0.0000178 #Diffusivité de O2 dans l'air
W=0.0000593 #Perméabilité de la membrane
Pi=150 #Pression partielle en mmHg
#C_0=5.8*(10e-6)*150 #Concentration en O2 en mol/m3
C_0=160

dx=(Lx/(n-1))
dy=(Ly/(n-1))

lamda=D/W
a=(Lx/Ly)
m=2 # Nombre de bronches ou membranes 

a=int(n/2)

b=int((3/10) * n)
c=int((3/10) * n)

d=int((2/10) * n)
e=int((2/10) * n)

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
        l.append([int((i+1/2)*(n/m)),random.randint(1,int((n/(2*m)))-1),random.randint(1,c+1)])
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


c=np.zeros((n*n,n*n))
h=np.zeros(n*n)

f=bijection(0,0)
c[f][f]=-(1+(2*a)-(a*(1+(dy/lamda))))
c[f][f+1]=1
c[f][f+n]=a
g=bijection(n-1,0)
c[g][g]=-(1+(2*a)-(a*(1+(dy/lamda))))
c[g][g-1]=1
c[g][g+n]=a

p=np.zeros((n*n))
for k in range(n*n):
    if A[k]==1:
        p[k]=1
        c[k][k]=-(2+(2*a))
        c[k][k+1]=1
        c[k][k-1]=1
        c[k][k+n]=a
        c[k][k-n]=a
    if A[k]==4 and A[k+n]==0:
        p[k]=2
        c[k][k]=-(2+a)
        c[k][k+1]=1
        c[k][k-1]=1
        c[k][k-n]=a
    if A[k]==4 and A[k+1]==0:
        p[k]=3
        c[k][k]=-(1+(2*a))
        c[k][k-1]=1
        c[k][k+n]=a
        c[k][k-n]=a
    if A[k]==4 and A[k-1]==0:
        p[k]=4
        c[k][k]=-(2+a)
        c[k][k+1]=1
        c[k][k+n]=a
        c[k][k-n]=a
    if A[k]==4 and A[k+1]==4 and A[k+n]!=0:
        p[k]=5
        c[k][k]=-(1+(2*a))
        c[k][k+1]=1
        c[k][k+n]=a
        c[k][k-n]=a
    if A[k]==4 and A[k+1]!=4 and A[k+1]!=0:
        p[k]=6
        c[k][k]=-(1+(2*a))
        c[k][k-1]=1
        c[k][k+n]=a
        c[k][k-n]=a
    if A[k]==0:
        p[k]=7
        c[k][k]=1
    if A[k]==2:
        p[k]=8
        c[k][k]=1
        h[k]=C_0
    if A[k]==3 and A[k-n]==0:
        p[k]=9
        c[k][k]=-(2+a-((a*dy)/lamda))
        c[k][k+1]=1
        c[k][k-1]=1
        c[k][k+n]=a
    if A[k]==3 and A[k-1]==0:
        p[k]=10
        c[k][k]=-(1+(2*a)-(dx/lamda))
        c[k][k+1]=1
        c[k][k-n]=a
        c[k][k+n]=a
    if A[k]==3 and A[k+1]==0:
        p[k]=11
        c[k][k]=-(1+(2*a)+(dx/lamda))
        c[k][k-1]=1
        c[k][k-n]=a
        c[k][k+n]=a
for i in range(0,n):
    k=bijection(i,0)
    if A[k]==3:
        p[k]=12
        c[k][k]=-(2+a-((a*dy)/lamda))
        c[k][k+1]=1
        c[k][k-1]=1
        c[k][k+n]=a

z=np.zeros((n,n))
s=np.zeros((n*n))
for i in range(n*n):
    if np.array_equal(s, c[i]):
        print(i)
# Séparer les coordonnées en listes i et j
i, j = zip(*indices)
# Plot en rouge (0) et en bleu (1)
plt.scatter(i, j, c=p, cmap='plasma', marker='o', s=50)
# Définir les étiquettes et le titre
plt.xlabel('Indice i')
plt.ylabel('Indice j')
plt.title('Plot de Maillage')
# Afficher la barre de couleur
plt.colorbar()
# Afficher le plot
plt.show()