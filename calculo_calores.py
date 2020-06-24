##Copyright INTA 
##Area de Materiales compuestos, Departamento de Materiales 2020
#Desarrollador Miguel Gonzalez del Val
# Referencias:
#             [1] Zarling 1981 Heat and mass transfer from freely falling drops at low temperatures
#             [2] NASA, Manual of Scaling Methods
#             [3] BERNARD L. MESSINGER, Equilibrium Temperature of an Unheated Icing Surface as a Function of Air Speed


import numpy as np
import pickle
import pandas as pd
class analisis_termico_flujo(object):
      
    def __init__(self):       
        self.cp_a = 0.24 #cal/g K
        self.tabla_presiones = pickle.load( open( "presiones_vapor.pickle", "rb" ) )
        self.t_f = 273.15 #punto de congelacion del agua
        self.tabla_hielo = pickle.load( open( "propiedades_hielo.pickle", "rb" ) )

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
        if (n>=1) or m<0: return 0
        else:return m
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

    def calculo_todos_calores(self,T_remanso,T_superficie,recovery_factor,p_remanso,LWC,cuerda,d,V,n_0,betha,m_in,T_s_anterior,cp_ws_anterior):
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

        
        T_estatica=T_remanso - V**2/1004.5 # K temperatura del flujo (estática)
        
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        Nu = self.Nusselt(Pr, Re_film)
        h_c = self.film_coefficient(k_a, Nu, d/100) #se pasa el diametro del perfil a metros
        q_convectivo = self.calor_convectivo(h_c, T_superficie, T_estatica, V)
        difusividad = self.difusividad(T_film,p_estatica) #cm^2/s
        difusividad=difusividad*10**(-4)
        Sc = self.Schmidt(mu_a*100,rho_a,difusividad)
        h_g = self.convective_mass_coeff(h_c, Sc, Pr) #g/(m^2 hr)
        Lambda_v =  self.Latent_heat_vaporisation(T_superficie) #cal/g
        p_w = self.presion_vapor(T_estatica)
        p_ww = self.presion_vapor(T_superficie)
        dot_me=self.mass_water_evaporation(p_ww,p_w,h_g,p_estatica,T_superficie,T_remanso,p_remanso,n_0)/3600
        
        q_evaporacion = self.heat_evaporation(dot_me,Lambda_v)
        dot_m = self.impinging_mass_water(LWC, V, betha) # g m^-2 s^-1
        cp_ws= self.cp_ws(T_superficie)
        q_w = self.Sensible_Heat_Water(cp_ws,dot_m+m_in,T_estatica)
        q_k = self.heat_kinetic(V,dot_m*10**-3) #cal/m^2 s
        cp_ws =self.cp_ws(T_superficie) ##cal/g K
        Lambda_f = self.latent_heat_freezing(T_superficie) #cal/g
        cp_is=self.cp_is(T_superficie)
        q_rb = self.calor_flujo_agua(n_0,dot_m+m_in,dot_me,cp_ws,T_superficie)
        q_i = self.heat_from_ice(dot_m+m_in,n_0,T_superficie,cp_is)
        q_f=self.heat_freezing(dot_m,Lambda_f,n_0)
        k_hielo = self.conductividad_hielo(T_superficie)
        q_in =self.heat_water_into(m_in,T_s_anterior,cp_ws_anterior)
        
        return (q_convectivo,q_evaporacion,q_w,q_rb,q_f,q_i,q_k,q_in)

    
    def calculo_todas_masas(self,T_remanso,T_superficie,recovery_factor,p_remanso,LWC,cuerda,d,V,n_0,betha,m_in,T_s_anterior,cp_ws_anterior):
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

        
        T_estatica=T_remanso - V**2/1004.5 # K temperatura del flujo (estática)
        
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        Nu = self.Nusselt(Pr, Re_film)
        h_c = self.film_coefficient(k_a, Nu, d/100) #se pasa el diametro del perfil a metros
        q_convectivo = self.calor_convectivo(h_c, T_superficie, T_estatica, V)
        difusividad = self.difusividad(T_film,p_estatica) #cm^2/s
        difusividad=difusividad*10**(-4)
        Sc = self.Schmidt(mu_a*100,rho_a,difusividad)
        h_g = self.convective_mass_coeff(h_c, Sc, Pr) #g/(m^2 hr)
        Lambda_v =  self.Latent_heat_vaporisation(T_superficie) #cal/g
        p_w = self.presion_vapor(T_estatica)
        p_ww = self.presion_vapor(T_superficie)
        dot_me=self.mass_water_evaporation(p_ww,p_w,h_g,p_estatica,T_superficie,T_remanso,p_remanso,n_0)/3600
        q_evaporacion = self.heat_evaporation(dot_me,Lambda_v)
        dot_m = self.impinging_mass_water(LWC, V, betha) # g m^-2 s^-1
        cp_ws= self.cp_ws(T_superficie)
        q_w = self.Sensible_Heat_Water(cp_ws,dot_m,T_estatica)
        q_k = self.heat_kinetic(V,dot_m*10**-3) #cal/m^2 s
        cp_ws =self.cp_ws(T_superficie) ##cal/g K
        Lambda_f = self.latent_heat_freezing(T_superficie) #cal/g
        cp_is=self.cp_is(T_superficie)
        q_rb = self.calor_flujo_agua(n_0,dot_m,dot_me,cp_ws,T_superficie)
        q_i = self.heat_from_ice(dot_m,n_0,T_superficie,cp_is)
        q_f=self.heat_freezing(dot_m,Lambda_f,n_0)
        k_hielo = self.conductividad_hielo(T_superficie)
        q_in =self.heat_water_into(m_in,T_s_anterior,cp_ws_anterior)
        m_out = (1-n_0)*(dot_m+m_in)-dot_me
        return (m_in,m_out,dot_m,dot_me)


    def balance_calores_reales(self,T_remanso,T_superficie,recovery_factor,p_remanso,LWC,cuerda,d,V,n_0,betha,m_in,T_s_anterior,cp_ws_anterior):
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

        
        T_estatica=T_remanso - V**2/1004.5 # K temperatura del flujo (estática)
        
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        Nu = self.Nusselt(Pr, Re_film)
        h_c = self.film_coefficient(k_a, Nu, d/100) #se pasa el diametro del perfil a metros
        
        difusividad = self.difusividad(T_film,p_estatica) #cm^2/s
        difusividad=difusividad*10**(-4)
        Sc = self.Schmidt(mu_a*100,rho_a,difusividad)
        h_g = self.convective_mass_coeff(h_c, Sc, Pr) #g/(m^2 hr)
        Lambda_v =  self.Latent_heat_vaporisation(T_superficie) #cal/g
        p_w = self.presion_vapor(T_estatica)
        p_ww = self.presion_vapor(T_superficie)
        dot_me=self.mass_water_evaporation(p_ww,p_w,h_g,p_estatica,T_superficie,T_remanso,p_remanso,n_0)/3600
        
        dot_m = self.impinging_mass_water(LWC, V, betha) # g m^-2 s^-1
        cp_ws= self.cp_ws(T_superficie)
        q_w = self.Sensible_Heat_Water(cp_ws,dot_m,T_estatica)
        q_k = self.heat_kinetic(V,dot_m*10**-3) #cal/m^2 s
        cp_ws =self.cp_ws(T_superficie) ##cal/g K
        Lambda_f = self.latent_heat_freezing(T_superficie) #cal/g
        cp_is=self.cp_is(T_superficie)
    
        k_hielo = self.conductividad_hielo(T_superficie)
        m_out = (1-n_0)*(dot_m+m_in)-dot_me
        return (m_in,m_out,dot_m,dot_me)

    def tiempo(self,Delta,LWC,V,betha_0):
        """[summary]

        Arguments:
            
            Delta {[type]} -- [capa hielo formada en metros]
            LWC {[type]} -- [contenido de agua liquida (g/m^3)]
            V {[type]} -- [velocidad del flujo (m/s)]
            betha_0 {[type]} -- [catch efficiency]

        Returns:
            [tao] -- [Tiempo de acreción de hielo en segundos]
        """ 
        rho_i =916.8e3       #g/m^3 (incompresible M<0.3)
        return rho_i*Delta/(LWC*V*betha_0)      


        

