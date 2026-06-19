import numpy as np
import matplotlib.pyplot as plt

n=200
Lx=10
Ly=10

dx=(Lx/(n-1))
dy=(Ly/(n-1))

C_0=160
lamda=10000000000
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
    c[i][i]=-(2+a-((a*dy)/lamda))
    c[i][i+1]=1
    c[i][i-1]=1
    c[i][i+n]=a
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
    c[k][k+n]=1
    c[k][k-n]=a
for j in range(1,n-1):
    k=bijection(n-1,j)
    c[k][k]=-(1+(2*a))
    c[k][k-1]=1
    c[k][k+n]=a
    c[k][k-n]=a

f=bijection(0,0)
c[f][f]=-(1+(2*a)-(a*(1+(dy/lamda))))
c[f][f+1]=1
c[f][f+n]=a

g=bijection(n-1,0)
c[g][g]=-(1+(2*a)-(a*(1+(dy/lamda))))
c[g][g-1]=1
c[g][g+n]=a

X=np.dot(np.linalg.inv(c),h)

for k in range(0,n*n):
    l=reciproque(k)
    i=l[0]
    j=int(l[1])
    b[j,i]=X[k]

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
