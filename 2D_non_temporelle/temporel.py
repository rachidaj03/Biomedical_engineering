import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

n=10
Lx=0.003*8
Ly=0.003*8
Lz=0.003*8
u0=1#valeur temporraire pour le moment


D=0.0000178 #Diffusivité de O2 dans l'air
W=0.0000593 #Perméabilité de la membrane
Pm=150 #Pression partielle en mmHg
C_0=5.8*(0.01)*Pm #Concentration en O2 en mol/m3
Pb=80
C_b=5.8*(0.01)*Pb #Concentration en O2 en mol/m3
#C_0=160

dx=(Lx/(n-1))
dy=(Ly/(n-1))
dt=1/1500

lamda=D/W
#a=(Lx/Ly)
alpha=dt*D/(2*(dx**2))
beta=1+dy/lamda
gama=u0*dt/dy

def bijection(i,j):
    return i+(n*j)

def reciproque(k):
    i=k%n
    j=(k-i)/n
    return [i,j]

 

c=np.zeros((n*n,n*n))
d=np.zeros((n*n,n*n))
Xn=np.zeros(n*n)
X0=np.zeros(n**2)+1

for i in range(0,n):
    K=bijection(i,n-1)
    X0[K]=C_0


# le remplissage de la matrice avec les coefficients du milieu du carré ou seule l'équation de diffusion est appliquée/////
for i in range(1,n-1):
    for j in range(1,n-1):
        k=bijection(i,j)
        c[k][k]=1+4*alpha
        c[k][k+1]=-alpha
        c[k][k-1]=-alpha
        c[k][k+n]=-alpha
        c[k][k-n]=-alpha
        d[k][k]=1-4*alpha
        d[k][k+1]=alpha
        d[k][k-1]=alpha
        d[k][k+n]=alpha
        d[k][k-n]=alpha
#Condition de Robin//////////
for i in range(1,n-1):
    k=bijection(i,0)
    c[k][k]=1+(4-beta)*alpha
    c[k][k+1]=-alpha
    c[k][k-1]=-alpha
    c[k][k+n]=-alpha
    d[k][k]=1+alpha*(beta-4)
    d[k][k+1]=alpha
    d[k][k-1]=alpha
    d[k][k+n]=alpha
   
   
# coefficient de concentration dans la limite Newman gauche///////////    
for j in range(1,n-1):
    k=bijection(0,j)
    c[k][k]=1+3*alpha
    c[k][k+1]=-alpha
    c[k][k+n]=-alpha
    c[k][k-n]=-alpha
    d[k][k]=1-3*alpha
    d[k][k+1]=alpha
    d[k][k+n]=alpha
    d[k][k-n]=alpha
# coefficient de concentration dans la limite Newman droite///////////
for j in range(1,n-1):
    k=bijection(n-1,j)
    c[k][k]=1+3*alpha
    c[k][k-1]=-alpha
    c[k][k+n]=-alpha
    c[k][k-n]=-alpha
    d[k][k]=1-3*alpha
    d[k][k-1]=alpha
    d[k][k+n]=alpha
    d[k][k-n]=alpha
   
#points particuliers////////////
    #point particulier à droite
k=bijection(0,0)
c[k][k]=1+(3-beta)*alpha
c[k][k+1]=-alpha
c[k][k+n]=-alpha
d[k][k]=1-3*alpha
d[k][k+1]=alpha
d[k][k+n]=alpha
d[k][k-n]=alpha
    #point particulier à gauche
k=bijection(n-1,0)
c[k][k]=1+(3-beta)*alpha
c[k][k-1]=-alpha
c[k][k+n]=-alpha
d[k][k]=1-3*alpha
d[k][k-1]=alpha
d[k][k+n]=alpha
d[k][k-n]=alpha

#Condition de dirichletXXXXXXXXXXXXXXX
def ConditionDirichlet(X,t):
    if 0<=t and t<=2:
        k=bijection(0,n-1)
        c[k][k]=1
        X[k]=C_0
        k=bijection(n-1,n-1)
        c[k][k]=1
        X[k]=C_0
    elif 2<t and t<=5:
        k=bijection(0,n-1)
        c[k][k]=1+2*alpha
        c[k][k+1]=-alpha
        c[k][k-n]=-alpha
        d[k][k]=1+2*((1/gama)-1)*alpha
        d[k][k+1]=alpha
        d[k][k-n]=alpha
        k=bijection(n-1,n-1)
        c[k][k]=1+2*alpha
        c[k][k-1]=-alpha
        c[k][k-n]=-alpha
        d[k][k]=1+2*((1/gama)-1)*alpha
        d[k][k-1]=alpha
        d[k][k-n]=alpha
   
    for i in range(1,n-1):
        k=bijection(i,n-1)
        if 0<=t and t<=2:
            c[k][k]=1
            X[k]=C_0
        elif 2<t and t<=5:
            c[k][k]=1+3*alpha
            c[k][k+1]=-alpha
            c[k][k-1]=-alpha
            c[k][k-n]=-alpha
            d[k][k]=1+((2/gama)-3)*alpha
            d[k][k+1]=alpha
            d[k][k-1]=alpha
            d[k][k-n]=alpha
for i in range(0,n):
    k=bijection(i,n-1)
    c[k][k]=1
    d[k][k]=1
    #Xn[k]=C_0
print(c)
print(d)
Xn=X0
def definitionXn1(C,X,t,c,X0,i,d):
    ConditionDirichlet(X,t)
    ic=np.linalg.inv(c)
    produit=np.dot(ic, d)
    puissance=np.linalg.matrix_power(produit, i)
    X=np.dot(puissance,X0)
    for k in range(n**2):
        l=reciproque(k)
        i=l[0]
        j=int(l[1])
        C[j,i]=X[k]
t=0 
time=10
i=0
Q=np.zeros((n,n))
# Visualizing with a plot

fig, axis = plt.subplots()
ti = np.linspace(0,time, 1500)
for k in range(n**2):
    l=reciproque(k)
    i=l[0]
    j=int(l[1])
    Q[j,i]=X0[k]
pcm = axis.pcolormesh(Q, cmap=plt.cm.jet)
flux_sortie=0
flux_entree=0
lentry=np.zeros(1500)
lsortie=np.zeros(1500)
def update(num):
    axis.clear()
    ic=np.linalg.inv(c)
    produit=np.dot(ic, d)
    puissance=np.linalg.matrix_power(produit, num)
    Xn=np.dot(puissance,X0)
    #puissance=np.linalg.matrix_power(produit, num)
    for k in range(n**2):
        l=reciproque(k)
        i=l[0]
        j=int(l[1])
        Q[j,i]=Xn[k]
    pcm = axis.pcolormesh(Q, cmap=plt.cm.jet)
    for i in range (n):
        k=bijection(i,n-1)
        flux_entrée=Xn[k]-Xn[k-n]
    flux_entrée=-D*flux_entrée*(dx/dy)*Lz
    for i in range (n):
         k=bijection(i,0)
         flux_sortie=Xn[k+n]-Xn[k]
    flux_sortie=D*flux_sortie*(dx/dy)*Lz
    lentry[num]=flux_entrée*30000*0.01
    lsortie[num]=flux_sortie*30000*0.01

   
plt.colorbar(pcm, ax=axis)  

ani = FuncAnimation(fig,update, frames=range(len(ti)), interval=10)

   



'''while t < time:
    #definitionXn1(Q,Xt1,t,c,X0,i,d)
    ic=np.linalg.inv(c)
    produit=np.dot(ic, d)
    puissance=np.linalg.matrix_power(produit, i)
    Xt1=np.dot(puissance,X0)
    for k in range(n**2):
        l=reciproque(k)
        i=l[0]
        j=int(l[1])
        Q[j,i]=Xt1[k]
   
    i+=1
    t+=dt
    pcm.set_array(Q)
    axis.set_title("Distribution at t: {:.3f} [s].".format(t))
    plt.pause(0.01)'''
plt.show()


plt.subplot(2, 1, 1)
plt.plot(ti,lentry)

# Add labels and title
plt.xlabel('Temps en second')
plt.ylabel('Débit entrée en mol/s')
plt.title('Débit en fonction du temps')
# Display the plot
#plt.show()
plt.subplot(2, 1, 2)
plt.plot(ti,lsortie)

# Add labels and title
plt.xlabel('Temps en second')
plt.ylabel('Débit sortie en mol/s')
#plt.title('Debit en fonction du temps')
# Display the plot
plt.show()