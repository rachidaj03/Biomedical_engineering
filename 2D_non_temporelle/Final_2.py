import numpy as np
import matplotlib.pyplot as plt

n=100
Lx=0.003*8
Ly=0.003*8
Lz=0.003*8

D=0.0000178 #Diffusivité de O2 dans l'air
W=0.0000593 #Perméabilité de la membrane
Pm=150 #Pression partielle en mmHg
C_0=5.8*Pm*0.01 #Concentration en O2 en mol/m3
Pb=38.5
C_b=5.8*0.01*Pb #Concentration en O2 en mol/m3
dx=(Lx/(n-1))
dy=(Ly/(n-1))


lamda=D/W
a=(Lx/Ly)

def bijection(i,j):
    return i+(n*j)
def reciproque(k):
    i=k%n
    j=(k-i)/n
    return [i,j]

c=np.zeros((n*n,n*n))
h=np.zeros(n*n)
b=np.zeros((n,n))

for i in range(0,n):
    k=bijection(i,n-1)
    c[k][k]=1
    h[k]=C_0

for i in range(1,n-1):
    c[i][i]=-(2+a*(1+(dy/lamda)))
    c[i][i+1]=1
    c[i][i-1]=1
    c[i][i+n]=a
    h[i]=-((dy*C_b)/lamda)
for i in range(1,n-1):
    for j in range(1,n-1):
        k=bijection(i,j)
        c[k][k]=-(2+(2*a))
        c[k][k+1]=1
        c[k][k-1]=1
        c[k][k+n]=a
        c[k][k-n]=a
for j in range(1,n-1):
    k=bijection(0,j)
    c[k][k]=-(1+(2*a))
    c[k][k+1]=1
    c[k][k+n]=a
    c[k][k-n]=a
for j in range(1,n-1):
    k=bijection(n-1,j)
    c[k][k]=-(1+(2*a))
    c[k][k-1]=1
    c[k][k+n]=a
    c[k][k-n]=a

f=bijection(0,0)
c[f][f]=-(1+(2*a)-(a*(1-(dy/lamda))))
c[f][f+1]=1
c[f][f+n]=a
h[f]=-((dy*C_b)/lamda)
g=bijection(n-1,0)
c[g][g]=-(1+(2*a)-(a*(1-(dy/lamda))))
c[g][g-1]=1
c[g][g+n]=a
h[g]=-((dy*C_b)/lamda)

X=np.linalg.solve(c,h)

for k in range(0,n*n):
    l=reciproque(k)
    i=l[0]
    j=int(l[1])
    b[j,i]=X[k]

flux_entrée=0
for i in range (n):
    k=bijection(i,n-1)
    flux_entrée=X[k]-X[k-n]
flux_entrée=-D*flux_entrée*(dx/dy)*Lz

flux_sortie=0
for i in range (n):
    k=bijection(i,0)
    flux_sortie=X[k+n]-X[k]
flux_sortie=-D*flux_sortie*(dx/dy)*Lz*30000

print(flux_entrée)
print("////////")
print(flux_sortie)

x = np.linspace(0, Lx, n)
y = np.linspace(0, Ly, n)
Z, Y = np.meshgrid(x, y)
# Plotting the results
#plt.pcolormesh(Z, Y, b, cmap='plasma')  # You can also use plt.pcolormesh for a different style
plt.contourf(Z, Y, b, cmap='plasma')  # You can also use plt.pcolormesh for a different style
plt.colorbar()  # Add a colorbar for reference
plt.title('Solution en concentration')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()