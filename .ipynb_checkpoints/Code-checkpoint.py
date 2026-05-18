import CoolProp.CoolProp as CP


# Q1.1

def iteration_T_cond(ecart_min,n_iterations):
    Tk = 273.15 + 25
    p_atm= 101300
    cp_air = CP.PropsSI('C', 'T', Tk-5, 'P', p_atm, 'Air')
    cp_co2_liq = CP.PropsSI('CP0MASS', 'T', Tk, 'Q', 0, 'CO2')
    L_v_co2 = CP.PropsSI('L', 'T', Tk, 'Q', 0, 'CO2')
    Tk1= 15 + 273.15 + 5 + ((L_v_co2 + 2*cp_co2_liq)/cp_air)
    k=0
    while(k<n_iterations and Tk1-Tk>ecart_min):
        Tk=Tk1
        cp_air = CP.PropsSI('C', 'T', Tk-5, 'P', p_atm, 'Air')
        cp_co2_liq = CP.PropsSI('CP0MASS', 'T', Tk, 'Q', 0, 'CO2')
        L_v_co2 = CP.PropsSI('L', 'T', Tk, 'Q', 0, 'CO2')
        Tk1= 15 + 273.15 + 5 + ((L_v_co2 + 2*cp_co2_liq)/cp_air)
    p_sat = CP.PropsSI('P', 'T', Tk1, 'Q', 1, 'CO2')
    return Tk1, p_sat

print(iteration_T_cond(1e-3,1000))    


#Q1.2
#Calcul de P_5s : 
print(CP.PropsSI('P', 'T', -8+273.15, 'Q', 1, 'CO2'))

#Calcul de T_2f :
print(294.82761632902583 + (CP.PropsSI('C', 'T', 294.82761632902583, 'P', 101300, 'Air')/CP.PropsSI('CP0MASS', 'T', 294.82761632902583, 'Q', 1, 'CO2')) * (20+273.15 - (294.82761632902583-5)))

#Calcul de q_c
print(-CP.PropsSI('C', 'T', 294.82761632902583, 'P', 101300, 'Air')*5)


# +
def recherche_P1(ecart_min,n_iterations): #servira pour le COP
    
def COP():
    Tcond,P3=iteration_T_cond(1e-3,1000)
    h3sv=CP.PropsSI('H', 'T', Tcond, 'Q', 1, 'CO2')
    h3sl=CP.PropsSI('H', 'T', Tcond, 'Q', 0, 'CO2')
    h4=CP.PropsSI('H', 'T', Tcond-2, 'Q', 0, 'CO2')
    s4=CP.PropsSI('S', 'T', Tcond-2, 'Q', 0, 'CO2')
    #h5=h4 #détente isenthalpique
    #h5sv=CP.PropsSI('H', 'T', -8+273.15 , 'Q', 1, 'CO2')
    cp_airA2 = CP.PropsSI('C', 'T', 20+273.15, 'P', 101300, 'Air')
    s2 = CP.PropsSI('S', 'T', 20 + 273.15, 'P', 101300, 'Air') #s2=sA2
    P2 = P3 #Pas de pertes de charges
    T2 = CP.PropsSI('T', 'P', P2, 'S', s2, 'CO2')
    T1 = -3 + 273.15
    P1 = P2 * (T1/T2) #Loi des gaz parfaits
    rho = 0.5*(CP.PropsSI('D', 'P', P1, 'T', T1, 'CO2') + CP.PropsSI('D', 'P', P2, 'T', T2, 'CO2')) #On calcule la masse volumique du CO2
    w_is = (1/rho)*(P2-P1) #car dh = Tds - vdP et on est dans le cas isentropique
    w = w_is/(1-0.121*(P2/P1))
    return -1 + cp_airA2*5/w

print(COP())


# +
#Q2.1
def calcul_COP(P3):
    cp_air = CP.PropsSI('C', 'T', 30+273.15, 'P', 101300, 'Air')
    h3 = CP.PropsSI('H', 'P', P3, 'T', 35+273.15, 'CO2')
    h4 = h3 #détente isenthalpique
    h2 = h3 + cp_air * 5
    s2 = CP.PropsSI('S', 'P', 101300, 'T', 35+273.15, 'Air')
    P2 = CP.PropsSI('P', 'S', s2, 'H', h2, 'CO2')
    T2 = CP.PropsSI('T', 'S', s2, 'H', h2, 'CO2')
    T1 = -8 + 5 + 273.15
    P1 = P2 * (T1/T2) #On applique la loi des gaz parfaits
    rho = 0.5*(CP.PropsSI('D', 'P', P1, 'T', T1, 'CO2') + CP.PropsSI('D', 'P', P2, 'T', T2, 'CO2')) #On calcule la masse volumique du CO2
    w_is = (1/rho)*(P2-P1) #on utilise pour cela que dh = Tds - vdP mais dans le cas isentropique
    w = w_is/(1-0.121*(P2/P1))
    return -1 + cp_air * 5/w

print(calcul_COP(120e5))
# -


