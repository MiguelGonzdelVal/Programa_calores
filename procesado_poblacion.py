import pandas as pd
import numpy as np
import fem_velocidades  
import matplotlib.pyplot as plt
from matplotlib import cm
import time
import seaborn as sns

#point = Point(0.5, 0.5)
#polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
Elementos = pd.read_csv('elementos.txt')
Velocidades = pd.read_csv('flujo_velocidades.axdt')
Superficie = pd.read_csv('perfil.csv')

Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie)
Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie)
Modelo.set_T_remanso(273.15-5)
Modelo.set_presion_remanso(1e5)
N_gotas =100
mu_D = 20e-6
sd_D = 5e-6
Diameters = np.abs(np.random.normal(mu_D, sd_D, N_gotas))
ns, bins, patches = plt.hist(x=Diameters*10**6, bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.title('Diameters of drops')
plt.xlabel(r'$\mu m$')
plt.show()
x_0 = min(Modelo.x_nodo) + 1e-3
x_f = 0.08
plt.figure(figsize =(15,10))
plt.plot(Modelo.x_superficie,Modelo.y_superficie,'o')
beta_x = []
U_d0 = 70
V_d0 = 0
#D =20e-6
N=1
print('calculando tiempo...')
tiempo_inicial=time.time()

for D in Diameters:
    #if (N*100//N_gotas)==(N*100/N_gotas):print(N*100/N_gotas)
    #print(N*100/N_gotas)
    n_nodos = 500
    for y_0 in np.linspace(-0.01,0.01,30):
        (x,y,U_d,V_d,t) = Modelo.trayectoria_gota(x_0,x_f,y_0,U_d0,V_d0,D,n_nodos)
        #trayectoria_gota(self,x_0,x_f,y_0,U_d0,V_d0,D,n_nodos)
        x_proyecccion =Modelo.proyeccion_gota(x_0,y_0,x_f,n_nodos)
        ds = np.sqrt((y[-1]-y[0])**2+(x[-1]-x_proyecccion)**2)
        dy = np.abs(y[-1]-y[0])
        betha = dy/ds
        if y[-1] >=0:beta_x.append([D,x[-1],betha,'extrados'])
        else:beta_x.append([D,x[-1],betha,'intrados'])
        #plt.plot(x,y,'k')
    print(str((time.time()-tiempo_inicial)/N*N_gotas/60)+' min')
    N=N+1


bethas = pd.DataFrame(data=beta_x,columns=['x','betha','zona'])
#bethas[bethas.x<0.1].plot.scatter(x='x',y='betha',c='zona',colormap='viridis')
g =sns.scatterplot(x="x", y="betha",
              hue="zona",
              data=bethas[bethas.x<0.04])  
plt.grid()