import pandas as pd
import numpy as np
import fem_velocidades  
import matplotlib.pyplot as plt
from matplotlib import cm
#point = Point(0.5, 0.5)
#polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
Elementos = pd.read_csv('elementos.txt')
Velocidades = pd.read_csv('flujo_velocidades.axdt')
Superficie = pd.read_csv('perfil.csv')
Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie)
Modelo.set_T_remanso(273.15-5)
Modelo.set_presion_remanso(1e5)
Y = np.linspace(-0.04,0.04,12)
x_0 = min(Modelo.x_nodo) + 1e-3
x_f = max(Modelo.x_nodo) - 1e-3
plt.figure()
plt.plot(Modelo.x_superficie,Modelo.y_superficie,'o')
for y in Y:
    (x_stream,y_stream)= Modelo.stream_line(x_0,y,x_f,1000)
    plt.plot(x_stream,y_stream,'k')
plt.show()
