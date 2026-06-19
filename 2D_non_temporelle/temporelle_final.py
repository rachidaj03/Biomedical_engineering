import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
import numpy as np
import random
from matplotlib import cm
from matplotlib.animation import FuncAnimation


######################### Description de la solution #############
''' On essaye premierement de construire une matrice notee A contenant des nombres de 0 a 7 suivant
la configuration dessous, cette matrice nous permet essentiellement de determiner le maillage tout en sachant l
etat de chaque points du maillage.
Deuxiement, en parcourant la matrice A, en essaye d'affecter les coefficient de l'equation pour chaque cas
dans les matrices C et B, tel que C*X=B  '''

# 0 Points non actifs
# 1 Points actifs sans condition limite
# 2 Condition en entrée
# 3 Condition de Robin
# 4 Condition en Neumann
# 5 Combinaison de Robin et Neumann
# 6 Combinaison de deux conditions Neumann suivant deux axes différents


# Les paramètres d'activation
n=101 # Nombre de points impair
Lx=0.003*8 # Longeur du domaine x
Ly=0.003*8 # Longeur du domaine y
Lz=0.003*8
m=4 # Nombre de bronches ou membranes
a=int(n/2)
b=int((2/10) * n)
c=int((5/10) * n)
d=int((2/10) * n)
e=int((2/10) * n)


# Les paramètres de concentration
D=0.0000178 #Diffusivité de O2 dans l'air
W=0.0000593 #Perméabilité de la membrane
Pp=150 #Pression partielle en mmHg
Pb=80 # Pression partielle en MmHg du sang
C_0=5.8*(0.01)*Pp #Concentration en O2 en mol/m3
#C_0=160
C_b=5.8*(0.01)*Pb #Concentration en sang en mol/m3
C_1=5.8*(0.01)*(21-17)*Pp/(2*100)
dx=(Lx/(n-1))
dy=(Ly/(n-1))
lamda=D/W
alpha=(Lx/Ly)
omega=2*3.14/5
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
#plt.subplot(2, 1, 1)
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
B[f]=-((dy*C_b)/lamda)

g=bijection(n-1,0)
C[g][g]=-(1+(2*alpha)-(alpha*(1-(dy/lamda))))
C[g][g-1]=1
C[g][g+n]=alpha
B[g]=-((dy*C_b)/lamda)


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
    # Concentration des points 4 (Condition en Robin)
    if A[k]==3:
        if A[k+1]==0:
            C[k][k]=-(1+2*alpha+(dx/lamda))
            C[k][k-n]=alpha
            C[k][k-1]=1
            C[k][k+n]=alpha 
            B[k]=-((dx*C_b)/lamda)
        elif A[k-1]==0:
            C[k][k]=-(1+2*alpha+(dx/lamda))###-(1+2*alpha-(dx/lamda))
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
    # Concentration des points 2 (Condition de Dirichlet)
    if A[k]==2:
        C[k][k]=1
        B[k]=1
    # Concentration des points 1 (Les points actifs sans condition limite)
    if A[k]==1:
        C[k][k]=-2*(1+alpha)
        C[k][k+1]=1
        C[k][k-1]=1
        C[k][k+n]=alpha
        C[k][k-n]=alpha
    # Concentration des points 0 (Les points non actifs)
    if A[k]==0:
        C[k][k]=1
    # Concentration des points 7 (Points de j=0 ayant une condition de Robin)
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

# Solution de l'équation linéaire
X=np.linalg.solve(C,B)

for i in range(n*n):
    if X[i]<0:
        print(i,'=',X[i])
t_values = np.linspace(0, 10, 100)  

# Reformulation de la matrice pour le plotage
for k in range(0,n*n):
    l=reciproque(k)
    i=l[0]
    j=int(l[1])
    z[j,i]=X[k]



x = np.linspace(0, Lx, n)
y = np.linspace(0, Ly, n)
Z, Y = np.meshgrid(x, y)

t = np.linspace(0,10, 10)
fig, ax = plt.subplots()
f = z * (C_0-C_b+C_1*np.cos(omega * t[0]))
pcm = ax.pcolormesh(Z, Y, f, cmap=cm.jet,vmin=0,vmax=C_0-C_b+C_1)
ax.set_title(f'f = X * (C_0-C_b+C_1cos({t[0]:.2f}))')
ax.set_xlabel('Lx')
ax.set_ylabel('Ly')
# Add colorbar
cbar = plt.colorbar(pcm, ax=ax)
lentry=np.zeros(10)
lsortie=np.zeros(10)

def update(num):
    ax.clear()
    f = z*(C_0-C_b+C_1*np.cos(omega*t[num]))
    pcm=ax.pcolormesh(Z, Y, f, cmap=cm.jet,vmin=0,vmax=C_0-C_b+C_1)
    flux_entrée=0

    for k in range(n*n):
        if A[k]==2:
            flux_entrée=(X[k]-X[k-n])*(C_0-C_b+C_1*np.cos(omega*t[num]))
    flux_entrée=D*flux_entrée*(dx/dy)*(8*Lz)  
    lentry[num]=flux_entrée*30000*0.001

    flux_sortie=0
    flux_sortie1=0
    flux_sortie2=0
    flux_sortie3=0
    for k in range(n*n):
        if A[k]==7 or (A[k]==3 and A[k-n]==0):
            flux_sortie1=(X[k+n]-X[k])*(C_0-C_b+C_1*np.cos(omega*t[num]))
        if A[k]==3 and A[k+1]==0:
            flux_sortie2=(X[k]-X[k-1])*(C_0-C_b+C_1*np.cos(omega*t[num]))
        if A[k]==3 and A[k-1]==0:
            flux_sortie3=(X[k+1]-X[k])*(C_0-C_b+C_1*np.cos(omega*t[num]))
    flux_sortie=(D*(-flux_sortie2+flux_sortie3)*(dy/dx)*(8*Lz)+D*flux_sortie1*(dx/dy)*(Lz))*30000*0.001
    lsortie[num]=flux_sortie

   
    ax.set_title(f'f = X * (C_0-C_b+C_1cos({t[num]:.2f}))')
    #plt2.plot(num,flux_entrée,flux_sortie)
    #cbar = plt.colorbar(pcm, ax=ax)
    #print(f)
    #cbar.clear()
   
'''for i in range(n):
    print(f[i])  ''' 

ani = FuncAnimation(fig, update, frames=range(len(t)), interval=1000)
ani.save('myAnimation.gif', writer='imagemagick', fps=30)

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


plt.subplot(2, 1, 1)
plt.plot(t,lentry)

# Add labels and title
plt.xlabel('Temps en second')
plt.ylabel('Débit entrée en mol/s')
plt.title('Débit en fonction du temps')
# Display the plot
#plt.show()
plt.subplot(2, 1, 2)
plt.plot(t,lsortie)

# Add labels and title
plt.xlabel('Temps en second')
plt.ylabel('Débit sortie en mol/s')
#plt.title('Debit en fonction du temps')
# Display the plot
plt.show()



 
