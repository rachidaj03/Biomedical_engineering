import numpy as np
import matplotlib.pyplot as plt

n=100
Lx=0.003
Ly=0.003
Lz=0.003

D=0.0000178 #Diffusivité de O2 dans l'air
W=0.0000593 #Perméabilité de la membrane
Pm=150 #Pression partielle en mmHg
C_0=5.8*(0.01)*Pm #Concentration en O2 en mol/m3

dx=(Lx/(n-1))
dy=(Ly/(n-1))
lamda=D/W
b=np.zeros((n,n))
x = np.linspace(0, Lx, n)
y = np.linspace(0, Ly, n)
Z, Y = np.meshgrid(x, y)
for i in range(n):
    for j in range(n):
        b[j,i]=C_0*(lamda+y[j])/(lamda+Lx)
        
# Plotting the results
#plt.pcolormesh(Z, Y, b, cmap='plasma')  # You can also use plt.pcolormesh for a different style
plt.contourf(Z, Y, b, cmap='plasma')  # You can also use plt.pcolormesh for a different style
plt.colorbar()  # Add a colorbar for reference
plt.title('Solution anayltique en 2D stationnaire')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()
