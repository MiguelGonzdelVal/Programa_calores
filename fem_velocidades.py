import numpy as np
from shapely.geometry import Point,LineString
from shapely.geometry.polygon import Polygon
import pickle
import pandas as pd
class modelo_fem(object):

    def __init__(self,Elementos,Velocidades,Superficie):
        #https://www.sciencedirect.com/topics/engineering/linear-triangular-element#:~:text=A%20linear%20triangular%20element%20is%20a%20two%2Ddimensional%20finite%20element,three%20sides%20shown%20in%20Fig.&text=Each%20linear%20triangular%20element%20has,X%20and%20Y%20axes%2C%20respectively.
        self.x_nodo = np.array(Velocidades['X[m](XCoordinate)'])
        self.y_nodo = np.array(Velocidades['Y[m](YCoordinate)'])
        self.u_nodo = np.array(Velocidades['Velocityu[ms^-1]'])
        self.v_nodo = np.array(Velocidades['Velocityv[ms^-1]'])
        self.Elementos = Elementos
        self.x_superficie = np.array(Superficie['X[m]'])
        self.y_superficie = np.array(Superficie['Y[m]'])
        
        self.rho_w = 1000 
        self.Superficie = Superficie
        self.perfil()
    

    def perfil(self):
        x = np.array(self.Superficie['X[m]'])
        y = np.array(self.Superficie['Y[m]'])
        coordenadas = []
        for i in range(len(x)):
            coordenadas.append((x[i],y[i]))
        self.polygon = Polygon(coordenadas)
        coordenadas = np.array(self.polygon.exterior.coords)
        self.exterior = LineString(coordenadas)
        normal = np.zeros((len(x),2))
        for i in range(1,len(x)):
            mod = np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)
            if mod > 0:
                if y[i]<0:
                    
                    normal[i,1] = -(x[i]-x[i-1])/mod
                    normal[i,0] = (y[i]-y[i-1])/mod
                else:
                    mod = np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)
                    normal[i,1] = (x[i]-x[i-1])/mod
                    normal[i,0] = -(y[i]-y[i-1])/mod

        self.normal=normal
        V_infinito = []
        for i in range(len(x)):
            if (normal[i,0]!=0) or (normal[i,1]!=0):
                V=self.Velocidad(x[i]+normal[i,0]*0.02,y[i]+normal[i,1]*0.02)
                V = np.sqrt(V[0]**2+V[1]**2)
                if V !=0:V_infinito.append([x[i],V])
                
        self.V_infinito = np.array(V_infinito)


    def Velocidad(self,x,y):
        """[Calcula la velocidad en las coordenadas x,y]

        Arguments:
            x {[float]} -- [coordenada x (m)]
            y {[float]} -- [coordenada y (m)]

        Returns:
            [tuple] -- [(u,v) del flujo]
        """        
        d=1e6
        U=0
        V=0
        for i in range(len(self.x_nodo)):
            dis =((self.x_nodo[i]-x)**2+(self.y_nodo[i]-y)**2)**.5
            if dis<d:
                U = self.u_nodo[i]
                V = self.v_nodo[i]
                d = dis
        return U,V

    def stream_line(self,x_0,y_0,x_f,n):
        """[Calcula las líneas de corriente del flujo]

        Arguments:
            x_0 {[float]} -- [coordenada x inicial en metros]
            y_0 {[float]} -- [coordenada y inicial en metros]
            x_f {[float]} -- [coordenada x final en metros]
            n {[int]} -- [numero de nodos]

        Returns:  (x,y)
            [tuple] -- [(x,y) del streamline]
        """        
        x = np.linspace(x_0,x_f,n)
        y = np.zeros(n)
        y[0] = y_0
        for i in range(1,n):
            (u,v) = self.Velocidad(x[i-1],y[i-1])
            if u!=0:
                y[i] = y[i-1] + v/u*(x[i]-x[i-1])
            else:   y[1] = y[0]
        return (x,y)

    def viscosity(self,T):
        """[viscosidad en función de la temperatura]

        Arguments:
            T {[float]} -- [Temperatura en K]

        Returns:
            [mu] -- [viscosidad en el SI]
        """        
        return 5e-8*T+3.46e-5

    def set_T_remanso(self,T_remanso):
        """[Fija la temperatura de remanso del problema]

        Arguments:
            T_remanso {[float]} -- [Temperatura (K)]
        """        
        self.T_remanso = T_remanso
    
    def Temperatura(self,u,v):
        """[Temperatura estatica]

        Arguments:
            u {[float]} -- [velocidad horizontal (m/s)]
            v {[float]} -- [velocidad horizontal (m/s)]

        Returns:
            [float] -- [Temperatura estatica (K)]
        """        
        V = np.sqrt(u**2+v**2)
        return self.T_remanso - 0.5*V**2/(2*1004.5)

    def set_presion_remanso(self,P_0):
        self.P_0 = P_0

    def presion_estatica(self,T,V):
        M =V/np.sqrt(1.4*287*T)
        gamma = 1.4
        return self.P_0/(1+(gamma -1)/2*M**2)**(gamma/(gamma-1))
    
    def densidad(self,T,P):
        return P/(287*T)

    def Reynolds(self,V,rho_a,mu_a,D):
        return rho_a*V*D/mu_a
    
    def CD(self,Re):
        if (Re<=1000)and (Re >1): 
            c=  24./Re+6./(1+np.sqrt(Re))+.27
        elif ((Re<=3600.) and (Re>1000)):
            c=0.6649+0.2712e-3*Re+1.22e-7*Re**2.-10.919e-12*Re**3.  
        elif (Re <=1): c=0
        return c 

    def trayectoria_gota(self,x_0,x_f,y_0,U_d0,V_d0,D,n_nodos):
        y = [y_0]
        U_d = [U_d0]
        V_d = [V_d0]
        t=[0]
        delta_x=(x_f-x_0)/(n_nodos-1.)
        x=[x_0]
        poly = self.polygon
        for i in range(1,n_nodos):
            x.append(x[-1]+delta_x)
            (u_a,v_a) = self.Velocidad(x[i-1],y[i-1])
            V_a = np.sqrt(u_a**2 + v_a**2)
            T_a=self.Temperatura(u_a,v_a)
            P_a = self.presion_estatica(T_a,V_a)
            mu_a=self.viscosity(T_a)
            rho_a = self.densidad(T_a,P_a)
            fg=9.81*(rho_a/self.rho_w-1)
            uf_x = u_a - U_d[i-1]
            uf_y = v_a - V_d[i-1]
            V_f = np.sqrt(uf_x**2 + uf_y**2)
            Re = self.Reynolds(V_a,rho_a,mu_a,D)
            f_d=3./4.*self.CD(Re)/D*rho_a/self.rho_w*V_f    
            #u_x(i)=delta_x/u_x(i-1)*f_d*uf_x+u_x(i-1)
            U_d.append(delta_x/U_d[i-1]*f_d*uf_x+U_d[i-1])
            #u_y(i)=delta_x/u_x(i-1)*(f_d*uf_y+fg)+u_y(i-1)
            V_d.append(delta_x/U_d[i-1]*(f_d*uf_y+fg)+V_d[i-1])
            t.append((x[i]-x[i-1])/U_d[i-1]+t[i-1])
            y.append(V_d[i-1]*(t[i]-t[i-1])+y[i-1])
            p1 = Point(x[-1],y[-1])
            if p1.within(poly):
                linea = LineString([(x[-2],y[-2]),(x[-1],y[-1])])
                Interseccion=self.exterior.intersection(linea)
                try:
                    x[-1] = Interseccion.x
                    y[-1] = Interseccion.y
                except:
                    
                    x[-1]=Interseccion[0].x
                    y[-1]=Interseccion[0].y
                break
        return (x,y,U_d,V_d,t)

    def proyeccion_gota(self,x_0,y_0,x_f,n_nodos):
        x= np.linspace(x_0,x_f,n_nodos)
        poly = self.polygon
        for i in range(n_nodos):
            p1 = Point(x[i],y_0)
            if p1.within(poly):
                break
        return x[i]

class analisis_termico(object):
    def __init__(self):       
        self.cp_a = 0.24 #cal/g K
        self.tabla_presiones = pickle.load( open( "presiones_vapor.pickle", "rb" ) )
        self.t_f = 273.15 #punto de congelacion del agua
        self.tabla_hielo = pickle.load( open( "propiedades_hielo.pickle", "rb" ) )
        self.T_remanso=None
        self.T_superficie= None
        self.recovery_factor=None 
        self.p_remanso= None
        self.LWC=None
        self.cuerda= None
        self.d=None
        self.V= None
        self.n_0=None
        self.betha=None
        self.m_in=None
        self.T_s_anterior=None
        self.cp_ws_anterior=None

    def set_T_remanso(self,T_0):
        self.T_remanso=T_0

    def set_recovery_factor(self,r):
        self.recovery_factor=r 

    def set_T_superficie(self,T_sur):
        self.T_superficie=T_sur

    def set_presion_remanso(self,P_0):
        self.p_remanso = P_0

    def set_LWC(self,LWC):
        self.LWC=LWC

    def set_cuerda(self,c):
        self.cuerda = c

    def set_diametro_caracteristico(self,d):
        self.d =d

    def set_velocidad_flujo(self,V):
        self.V=V

    def set_coeficiente_convectivo(self,h_c):
        self.h_c = h_c

    def set_freezing_fraction(self,n):
        self.n_0 = n

    def set_local_collection_efficiency(self,betha):
        self.betha = betha

    def set_flujo_masico_entrada(self,m_in):
        self.m_in =m_in

    def set_T_superficie_anterior(self,T_sur_i_menos1):
        self.T_s_anterior=T_sur_i_menos1

    def set_cp_ws_anterior(self,cp_ws_i_menos1):
        self.cp_ws_anterior = cp_ws_i_menos1

    def set_V_exterior(self,Ve):
        self.V_e
    
    def __str__(self):
        texto ='Copyright INTA'
        texto = texto+'\nArea de Materiales compuestos, Departamento de Materiales 2020'
        texto = texto+'\nDesarrollador Miguel Gonzalez del Val'
        texto = texto+'\nReferencias:'
        texto = texto+'\n[1] Zarling 1981 Heat and mass transfer from freely falling drops at low temperatures'
        texto = texto+'\n[2] NASA, Manual of Scaling Methods'
        texto = texto+'\n[3] BERNARD L. MESSINGER, Equilibrium Temperature of an Unheated Icing Surface as a Function of Air Speed'
        return texto

    def set_calor_conductivo(self,calor_conductivo):
	    self.calor_conductivo = calor_conductivo
    

    def calor_convectivo(self,h_c,T_superficie,T_estatica,V):
        '''
        Calcula el calor convectivo en la linea de remanso (manual
        scaling methods eq. 3.53)
        variables:
            h_c: film coefficient (cal/(hr m^2 K)
            T_superficie: K
            T_estatica: K
            V:Velocidad del flujo
            cp_a: coeficiente calorifico a presion cte
            q_h (cal/hr m^2)
        '''
        
        return h_c*(T_superficie-T_estatica-V**2/(2*1004.5))/3600
    def film_coefficient(self,k_a,Nu_a,d):
        '''
        Calcula el coeficiente convectivo
        Variables:
            ka: coeficiente de conductivo del aire (cal/(m K hr)
            Nu_a: Nº de Nusselt del aire
            d: diametro caracteristico del perfil (borde de ataque) (metros)
        '''
        return k_a*Nu_a/d

    def Nusselt(self,Pr,Re):
        '''
        Ecuación (3.33) Manual de Scaling Methods (NASA)
        Pr: Nº de Prandtl
        Re: Nº De Reynolds
        '''
        return 1.14*Pr**0.4*Re**0.5
    def Pr(self,cp_a,mu_a,k_a):
        '''
        Numero de Prandtl (ecuacion 3.34)
        c_p: coeficiente calorifico a presion cte 
        mu_a  viscosidad
        k_a: conductividad termica
        '''
        return cp_a*mu_a/k_a

    def Reynolds(self,V,c,rho,mu_a):
        '''
        Numero de Reynolds (ecuacion 3.39)
        Se hace respecto a la cuerda
        rho_a: densidad
        c: cuerda del perfil
        mu_a viscosidad del aire
        V_velocidad del flujo
        '''
        return V*c*rho/mu_a
    def thermal_conductivity(self,T_film):
        '''
        Conductividad térmica del aire. 
        es una funcion de la temperatura (ver ecuacion A.3 de Manual Scaling Methods NASA)
        '''
        k_a =-12.69 + 2.029 *np.sqrt(T_film) # cal/(hr m K)
        return k_a

    def Temperatura_film(self,T_st,T_s):
        '''
        Ecuacion 3.35 de Manual Scaling methods
        '''
        return 0.5*(T_s+T_st)
    def convective_mass_coeff(self,h_c,Sc,Pr):
        '''
        Calcula el coeficiente de transmisión de masa
        Ecuacion A35 de manual scaling methods
        '''
        
        return h_c/self.cp_a*(Pr/Sc)**.67

    def difusividad(self,T_film,p_st):
        '''
        Calcula la difusividad en funcion de la temperatura. 
        Ecuacion A.4 Manual de Scaling Methods
        
        '''
        
        return 0.211*(T_film/273.15)**1.94*(101325/p_st)

    def Schmidt(self,mu_a,rho_a,Difusividad):
        '''
        Numero de Schmidt (Sc)
        Difusividad del aire, mu_a viscosidad y rho_a es la densidad
        Ecuacion  3.48 Manual Scaling NASA
        '''
        return mu_a/(rho_a*Difusividad)


    def presion_vapor(self,T):
        '''
        Presion de vapor del agua en funcion del ratio de humedad
        Datos de Beltramino, G., Rosso, L., Cuccaro, R., Tabandeh, S., Smorgon, D., & Fernicola, V. (2019).
        Accurate vapour pressure measurements of super cooled water in the temperature range between 252 K and 273 K. 
        The Journal of Chemical Thermodynamics, 105944. doi:10.1016/j.jct.2019.105944 
        
        '''
        Temperatura =self.tabla_presiones[:,0]
        P_vapor = self.tabla_presiones[:,2]
        p = P_vapor[-1]
        if T > Temperatura[-1]: p=P_vapor[-1]
        for i in range(1,len(Temperatura)):
            if (T>=Temperatura[i-1]) and (T<=Temperatura[i]):
                p = P_vapor[i-1]+ (P_vapor[i]-P_vapor[i-1])/(Temperatura[i]-Temperatura[i-1])*(T-Temperatura[i-1])
                break
        return p

    def mass_water_evaporation(self,p_ww,p_w,h_g,p_estatica,T_surface,T_remanso,p_remanso,n):
        '''
        p_ww,p_w,h_g,p_estatica,T_estatica,T_remanso,p_remanso,0.5
        Calcula la masa de aire que se evapora (g/(hr m^2)).
        Ecuacion 3.46 de Manual Scaling
        h_g constante de transmision de masa por evaporacion (#g/(m^2 hr) )
        p_ww Presion parcial de vapor del agua en la superficie
        p_w Presion parcial de vapor del agua en el flujo externo
        p_st presion estatica
        n: freezing factor
        Messiger  dice que si n=1 no hay evaporacion (m_e=0). 
        '''
       
        m = h_g*((p_ww/T_surface-p_remanso/T_remanso*p_w/p_estatica)/(p_remanso/T_remanso/.622-p_ww/T_surface))
        if n>=1 or m<0:m=0

        return m
    def heat_evaporation(self,m_e,Lambda_v):
        '''
        Calcula el flujo de calor por evaporacion (cal/(m^2 hr)).
        m_e : la masa de aire que se evapora (g/(hr m^2))
        Lambda_v = calor latente de evaporacion  (cal/g)
        Se toma como incompresible
        '''
        return m_e*Lambda_v

    def viscosity(self,T):
        '''
        #g/(cm s)
        Viscosidad del aire en funcion de la temperatura
        '''    
        return 10**(-4)/(.12764+124.38/T)#g/(cm s)

    def Latent_heat_vaporisation(self,T):
            '''
            calor latente de vaporizacion (cal/g) en funcion de la Temperatura
            de la superficie (glaze T=273.15 K) 
            T en Kelvin
            '''
            E =  0.197 + 3.670e-4*T
            return 597.3*(273.15/T)**E
        

    def cp_ws(self,T_superficie):
        '''
        Calcula el calor especifico del agua en la superficie 
        unidades cal/(g K)
        cp = 8.29e-5*(T_superficie.273.15)**2
        '''
        cp = 1.0074+8.29e-5*(T_superficie-273.15)**2
        return cp

    def impinging_mass_water(self,LWC,V,betha_0):
        '''
        impinging_mass_water Masa de agua que choca
        Ecuacion 3.50 Manual scaling methods
        LWC: liquid water content
        V: Velocidad del flujo
        betha_0 = catch efficiency at stagnation (eficiencia de adhesion en remanso)
        '''
        return LWC*V*betha_0

    def modified_inertia_parameter(self,LAMBDA_LAMBDAstokes,K):
        '''
        Langmuir and Blodgett’s expression for modified inertia
        parameter (eq. 3.8 Scaling methods)
        '''
        if K>1/8:return 1/8+ LAMBDA_LAMBDAstokes*(K-1/8)
        else:print('K no valida')
    def inertia_parameter(self,rho_w,delta,V,d,mu_a):
        '''
        K=rho_w*delta**2*V/(18*d*mu_a)
        K is the non-dimensional inertia parameter
        defined by Langmuir and Blodgett  Eq. 3.5
        rho_w, densidad del agua
        delta: diametro gota
        d: diametro curvatura del borde de ataque del perfil
        mu_a: viscosidad
        '''
        return rho_w*delta**2*V/(18*d*mu_a)
        
    def dimensionless_range_parameter(self,V,delta,rho,mu_a):
        '''
        define el parametro lambda/lambda_stokes (equation 3.9 Manual scaling methods)
        delta (m): tamaño de la gota
        V (m/s): velocidad
        mu_a (g/m s) viscosidad
        rho (g/m^3)
        '''
        Re = self.Reynolds(V,delta,rho,mu_a)
        return 1/(.8388+0.001483 *Re +.1847*Re**0.5)

    def catch_efficiency_stagnation(self,K_0):
        '''
        Ecuacion 3.13 del Manual of Scaling Methods
        (betha_0)
        K_0 es el parametro de inercia (ver inertia_parameter)
        '''
        return 1.4*(K_0-1/8)**.84/(1.4*(K_0-1/8)**.84+1)

    def Sensible_Heat_Water(self,cp_ws,dot_m,t_st):
        ''' 
        Ecuacion 5 del apartado 3.5 Manual Scaling methods
        cp,ws: calor especifico del agua a la temperatura de la superficie,
        dot_m: flujo masico del agua que impregna la superficie
        t_st temperatura estatica
        '''
        return dot_m*cp_ws*(self.t_f-t_st)


    def heat_kinetic(self,V,dm_dt):
        '''
        calor debido a Energia cinetica del flujo (cal/(m^2 s))S
        ver termino (10) del apartado 3.5
        dm_dt (kg/(s m^2):impinging_mass_water
        '''
        return dm_dt*V**2/2/4.1868 #cal/m^2 s

    def Temperatura_pared(self,r,V,T_inf):
        '''
        FUente: http://www.thermopedia.com/content/291/
        Se supone flujo laminar
            Pr:Numero de Prandtl
            T_inf: temperatura estatica del flujo exterior
            V Velocidad del flujo
             For some simple cases, its value can be estimated as follows: 
                At the front stagnation point of bodies in the flow, r = 1;
                in a laminar boundary layer on a plane plate, r =Pr^0.5 for Prandtl numbers 0.5 < Pr <10; 
                in a turbulent boundary layer on a plate, r =Pr^0.33 for Prandtl numbers close to 1
        '''
        return T_inf+r*V**2/(2*self.cp_a)
    
    def latent_heat_freezing(self,T_superficie):
        '''
        Ecuacion A19 Manual Scaling Methods
        Unidades: cal/g
        T_superficie (K): en glaze es 0 K
        79.7 + .485*(T_superficie-273.15)-2.5e-3*(T_superficie-273.15)**2 #cal/g
        '''
        return 79.7 + .485*(T_superficie-273.15)-2.5e-3*(T_superficie-273.15)**2 #cal/g
    def relative_heat_factor(self,LWC,cp_ws,betha_0,V,h_c):
        '''
        Ecuacion 3.55 del Manual Scaling Methods
        LWC (g/m^3)
        V (m/s)
        cp_ws (calor especifico del agua en la superficie) (cal/g)
        h_c: cal/(s m^2) coeficiente de transmision de calor convectivo
        b=LWC*V*betha_0*cp_ws/h_c
        '''
        return LWC*V*betha_0*cp_ws/h_c

    def drop_energy_transfer(self,T_estatica,V,cp_ws):
        '''
        ecuacion 3.56
        phi = T_f - T_estatica - V**2/(2*cp_ws)
        V(m/s)
        Temperaturas en K
        cp_ws (calor especigico del agua en la superficie J/(kg K))
        '''
        return self.t_f - T_estatica - V**2/(2*cp_ws)
        
    def air_energy_transfer(self,T_superficie,T_remanso,V,h_g,h_c,p_ww,p_w,p_estatica,p_remanso,T_estatica,Lambda_v):
        '''
        Ecuacion A52 
        Air energy transfer
        pww: Pa presion de vapor del agua en la superficie
        P_st: Pa presion estatica
        p_w: Pa p vapor del agua en el flujo
        p_remanso: Pa
        T_remanso:
        T_st
        T_superficie
        P_remanso
        h_g:coeficiente de transmision de masa por evaporacion (cal/(hr m^2 hr))
        h_c:coeficiente de transmision de calor por conveccion (cal/(hr m^2 hr))
        cp_a : J/Kg K
        Lambda_v: calor por evaporacion (cal/g)
        '''
        termino_1 = T_superficie - T_estatica-V**2/(2*1004.5) #K
        termino_2 = h_g/h_c*Lambda_v*((p_ww/T_estatica-p_remanso/T_remanso*p_w/p_estatica)/(
            p_remanso/T_remanso/.622-p_ww/T_estatica))
        return termino_1+termino_2

    def freezing_fraction_stagnation(self,cp_ws,Lamdbda_f,phi,theta,b):
        '''
        ecuacion 3.59
        n_0 es un numero adimensional que describe la cantidad de hielo que congela al llegar a la superficie
        Lambda_f es el calor latente de congelacion (cal/g)
        cp_ws calor especifico del agua en la superficie (cal/(K g))
        phi: drop_energy_transfer en Kelvin
        theta: air_energy_transfer en Kelvin
        b:relative_heat_factor adimensional
        '''
        return cp_ws/Lamdbda_f*(phi+theta/b)


    def calor_flujo_agua(self,n,dot_m,dot_me,cp_ws,T_superficie):
        '''
        Flujo de calor debido al flujo de agua que se escapa del volumen de control 
        unidades: cal/(m^2 s)
        Ecuacion 6 del apartado 3.5
        dot_m flujo de agua que impregna la superficie (g/(s m^2))
        dot_me flujo de agua que se evapora (g/(s m^2))
        cp_we: calor especifico del agua en la superficie (cal/(g K))
        n : ratio de congelacion
        T_superficie: K
        '''
        return ((1-n)*dot_m-dot_me)*cp_ws*(self.t_f-T_superficie)

    def cp_is(self,T_superficie):
        '''
        calor especifico del hielo (cal/g K)
        '''
        Temperatura =self.tabla_hielo[:,0]
        Cp = self.tabla_hielo[:,3]
        T=T_superficie-273.15
        cp = Cp[0]
        if T >0: cp = Cp[0]
        for i in range(1,len(Temperatura)):
            if (T>=Temperatura[i]) and (T<=Temperatura[i-1]):
                cp = Cp[i-1]+ (Cp[i]-Cp[i-1])/(Temperatura[i]-Temperatura[i-1])*(T-Temperatura[i-1])
                #cp = Cp[i-1]
                break
        return cp/4.1848

    def conductividad_hielo(self,T_superficie):
        '''
        calor especifico del hielo (cal/m s K)
        '''
        
        Temperatura =self.tabla_hielo[:,0]
        K = self.tabla_hielo[:,2]
        T=(T_superficie)/2-273.15
        k= K[0]
        for i in range(1,len(Temperatura)):
            if (T>=Temperatura[i]) and (T<=Temperatura[i-1]):
                k = K[i-1]+ (K[i]-K[i-1])/(Temperatura[i]-Temperatura[i-1])*(T-Temperatura[i-1])
                #cp = Cp[i-1]
                break
        return k/4.1848
    
    def heat_from_ice(self,dot_m,n,T_superficie,cp_is):
        '''
        flujo de calor del hielo
        q_i = dot_m*n*cp_is*(self.t_f-T_superficie)
        q_i: unidades 
        cal/(m^2 s)
        dot_m (g/s m^2)
        n adimensional
        T_superficie K
        cp_is:cal/(g K)
        '''
        return dot_m*cp_is*n*(self.t_f-T_superficie)

    def heat_freezing (self,dot_m,Lambda_f,n):
        '''
        ecuacion 8 de la seccion 3.5 manual scaling methods
        cal/(m^2 s)
        '''
        return n*dot_m*Lambda_f

    def heat_conduction(self,Delta,chi,l,k_i,T_surface,r,V,T_estatica):
        """[Calor perdido por el modelo por conducción]

        Arguments:
            Delta {[Float]} -- [Espesor de la capa de hielo (m)]
            chi {[Float]} -- [arco característico de la superficie de control (m)]
            l {[Float]} -- [longitud característica de la superficie de control (m)]
            k_i {[Float]} -- [Conductividad térmica hielo (cal /(m K s)]
            T_surface {[Float]} -- [Temperatura de la superficie (K)]
            r {[Float]} -- [factor de relajación]
            V {[Float]} -- [Temperatura de la superficie (K)]
            T_estatica {[Float]} -- [Temperatura de la superficie (K)]

        Returns:
            [q_cond] -- [calor de conduccion (cal/(m2 s))]
        """         
        if self.calor_conductivo:
            calor = k_i*Delta/chi/l*(T_surface-T_estatica-r*V**2/(1004.5*2))
        else: calor=0
        return calor

    def coeficiente_convectivo_lineal(self,a,b,T_film):
        return a*T_film+b
    def balance_calores(self,q_c,q_e,q_w,q_rb,q_f,q_i,q_k,q_in):
        '''
        Calculo de calores del Manual de Scaling methods
        ver seccion 3.5
        los calores deben ir en las mismas unidades
        predeterminadas: cal/(m^2 s)
        esta ecuacion deberia ser cero
        
        '''
        calores_perdidos = q_c +q_e +q_w+q_rb
        calores_ganados = q_f+q_i+q_k+q_in
        return calores_ganados-calores_perdidos

    def heat_water_into(self,dot_mi,T_superficie_anterior,cp_ws_anterior):
        return dot_mi*cp_ws_anterior*(T_superficie_anterior-273.15)

    def set_tamano_gota(self,D_gota):
        self.D_gota =D_gota
    def calculo_todos_calores(self):
        """[Calcula todos los calores del problema. Ver la siguiente bibliografia:
            NASA, Manual of Scaling Methods
            BERNARD L. MESSINGER, Equilibrium Temperature of an Unheated Icing Surface as a Function of Air Speed]

        Arguments:
            T_pared_inicial {[float]} -- [Temperatura de la pared antes de la nebulizacion (K)]
            T_modelo {[float]} -- [Temperatura de la pared durante la nebulizacion (K)]
            T_superficie {[float]} -- [Temperatura de la superficie del hielo durante la nebulizacion (K)]
            recovery_factor {[float]} -- [factor de recuperacion ]
            p_remanso {[float]} -- [Presión de remanso del problema (Pa) suele ser 10^5]
            LWC {[float]} -- [contenido de agua líquida (g/m^3)]
            cuerda {[float]} -- [cuerda del perfil (cm)]
            D_gota {[float]} -- [diametro de la gota (m)]
            d {[float]} -- [diametro caracteristico del perfil (cm)]
            Grosor_hielo {[float]} -- [Maximo espesor que se va a computar (m)]
            V {[float]} -- [Velocidad del flujo (m/s)]
            n_0 {[float]} -- [freezing fraction (adimensional)]

        Returns:
            [tuple] -- [(q_convectivo,q_evaporacion,q_w,q_rb,q_cond,q_f,q_i,q_k) calores (cal/(m^2 s))]
        """        
        T_remanso = self.T_remanso
        T_superficie = self.T_superficie
        recovery_factor=self.recovery_factor
        p_remanso=self.p_remanso
        LWC=self.LWC
        d=self.d
        V=self.V
        n_0=self.n_0
        betha=self.betha
        m_in=self.m_in
        T_s_anterior=self.T_s_anterior
        cp_ws_anterior = self.cp_ws_anterior
        V_e =self.V_e
        T_estatica =self.T_estatica
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film) # cal/(hr m K)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        self.Prandtl = Pr
        Nu = self.Nusselt(Pr, Re_film)
        #h_c = self.film_coefficient(k_a, Nu, d/100) #cal/(m K hr)
        #h_c = h_c/3600 #cal/(m K s)
        h_c = self.h_c
        difusividad = self.difusividad(T_film,p_estatica) #cm^2/s
        difusividad=difusividad*10**(-4)
        Sc = self.Schmidt(mu_a*100,rho_a,difusividad)
        h_g = self.convective_mass_coeff(h_c, Sc, Pr) #g/(m^2 s)
        Lambda_v =  self.Latent_heat_vaporisation(T_superficie) #cal/g
        p_w = self.presion_vapor(T_estatica)
        p_ww = self.presion_vapor(T_superficie)
        dot_me=self.mass_water_evaporation(p_ww,p_w,h_g,p_estatica,T_superficie,T_remanso,p_remanso,n_0) #g/(m^2 s)
        self.m_e =dot_me
        dot_m = self.impinging_mass_water(LWC, V, betha) # g m^-2 s^-1
        alpha_a = k_a/3600/rho_a/self.cp_a
        self.m_c = dot_m
        cp_ws =self.cp_ws(T_superficie) ##cal/g K
        Lambda_f = self.latent_heat_freezing(T_superficie) #cal/g
        cp_is=self.cp_is(T_superficie) #(cal/g K)
        q_1 = dot_m*(cp_ws*(T_estatica-273.15)+V**2/2/4.1868/1000)
        q_2 = m_in*cp_ws_anterior*(T_s_anterior-273.15)
        #q_3 = dot_me*(cp_ws*(T_superficie-273.15)+Lambda_v)
        q_3=0
        m_out = (1-n_0)*(dot_m+m_in)-dot_me
        if m_out < 0: 
            m_out=0
        self.m_out = m_out
        
        q_4 = m_out*cp_ws*(T_superficie-273.15)
        q_5 = n_0*(dot_m+m_in)*(cp_is*(T_superficie-273.15)-Lambda_f)
        q_6 = h_c*(T_superficie-T_estatica-recovery_factor*V_e**2/(2*1004.5))
        return q_1+q_2-q_3-q_4-q_5-q_6


    def calculo_todos_calores_d(self,d):
        self.set_diametro_caracteristico(d)
        T_remanso = self.T_remanso
        T_superficie = self.T_superficie
        p_remanso=self.p_remanso
        d=self.d
        V=self.V
        T_estatica =self.T_estatica
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        self.Prandtl = Pr
        Nu = self.Nusselt(Pr, Re_film)
        h_c = self.film_coefficient(k_a, Nu, d/100) #cal/(m K hr)
        self.set_coeficiente_convectivo(h_c)
        return self.calculo_todos_calores()


    def calculo_factor_coleccion(self,D_gota):
        T_remanso = self.T_remanso
        T_superficie = self.T_superficie
        p_remanso=self.p_remanso
        LWC=self.LWC
        d=self.d
        V=self.V
        T_estatica =self.T_estatica
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        self.Prandtl = Pr
        K=self.inertia_parameter(1e6,D_gota, V, d, mu_a)
        LAMBDA_LAMBDAstokes= self.dimensionless_range_parameter(V,self.D_gota,rho_a,mu_a*100)
        K_0 = self.modified_inertia_parameter(LAMBDA_LAMBDAstokes, K)
        betha_0 =self.catch_efficiency_stagnation(K_0)
        self.betha_stagnation=betha_0
        return betha_0

    def calculo_todos_calores_HCT(self):
        """[Calcula todos los calores del problema. Ver la siguiente bibliografia:
            NASA, Manual of Scaling Methods
            BERNARD L. MESSINGER, Equilibrium Temperature of an Unheated Icing Surface as a Function of Air Speed]

        Arguments:
            T_pared_inicial {[float]} -- [Temperatura de la pared antes de la nebulizacion (K)]
            T_modelo {[float]} -- [Temperatura de la pared durante la nebulizacion (K)]
            T_superficie {[float]} -- [Temperatura de la superficie del hielo durante la nebulizacion (K)]
            recovery_factor {[float]} -- [factor de recuperacion ]
            p_remanso {[float]} -- [Presión de remanso del problema (Pa) suele ser 10^5]
            LWC {[float]} -- [contenido de agua líquida (g/m^3)]
            cuerda {[float]} -- [cuerda del perfil (cm)]
            D_gota {[float]} -- [diametro de la gota (m)]
            d {[float]} -- [diametro caracteristico del perfil (cm)]
            Grosor_hielo {[float]} -- [Maximo espesor que se va a computar (m)]
            V {[float]} -- [Velocidad del flujo (m/s)]
            n_0 {[float]} -- [freezing fraction (adimensional)]

        Returns:
            [tuple] -- [(q_convectivo,q_evaporacion,q_w,q_rb,q_cond,q_f,q_i,q_k) calores (cal/(m^2 s))]
        """  
        a=self.a
        b=self.b
        
        T_remanso = self.T_remanso
        T_superficie = self.T_superficie
        recovery_factor=self.recovery_factor
        p_remanso=self.p_remanso
        LWC=self.LWC
        d=self.d
        V=self.V
        n_0=self.n_0
        betha=self.betha
        m_in=self.m_in
        T_s_anterior=self.T_s_anterior
        cp_ws_anterior = self.cp_ws_anterior
        V_e =self.V_e
        T_estatica =self.T_estatica
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film) # cal/(hr m K)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        self.Prandtl = Pr
        Nu = self.Nusselt(Pr, Re_film)
        h_c=self.coeficiente_convectivo_lineal(a,b,T_superficie)
        difusividad = self.difusividad(T_film,p_estatica) #cm^2/s
        difusividad=difusividad*10**(-4)
        Sc = self.Schmidt(mu_a*100,rho_a,difusividad)
        h_g = self.convective_mass_coeff(h_c, Sc, Pr) #g/(m^2 s)
        Lambda_v =  self.Latent_heat_vaporisation(T_superficie) #cal/g
        p_w = self.presion_vapor(T_estatica)
        p_ww = self.presion_vapor(T_superficie)
        dot_me=self.mass_water_evaporation(p_ww,p_w,h_g,p_estatica,T_superficie,T_remanso,p_remanso,n_0) #g/(m^2 s)
        self.m_e =dot_me
        dot_m = self.impinging_mass_water(LWC, V, betha) # g m^-2 s^-1
        alpha_a = k_a/3600/rho_a/self.cp_a
        self.m_c = dot_m
        cp_ws =self.cp_ws(T_superficie) ##cal/g K
        Lambda_f = self.latent_heat_freezing(T_superficie) #cal/g
        cp_is=self.cp_is(T_superficie) #(cal/g K)
        q_1 = dot_m*(cp_ws*(T_estatica-273.15)+V**2/2/4.1868/1000)
        q_2 = m_in*cp_ws_anterior*(T_s_anterior-273.15)
        q_3 = dot_me*(cp_ws*(T_superficie-273.15)+Lambda_v)
        m_out = (1-n_0)*(dot_m+m_in)-dot_me
        if m_out < 0: 
            m_out=0
        self.m_out = m_out
        
        q_4 = m_out*cp_ws*(T_superficie-273.15)
        q_5 = n_0*(dot_m+m_in)*(cp_is*(T_superficie-273.15)-Lambda_f)
        q_6 = h_c*(T_superficie-T_estatica-recovery_factor*V_e**2/(2*1004.5))
        return q_1+q_2-q_3-q_4-q_5-q_6
