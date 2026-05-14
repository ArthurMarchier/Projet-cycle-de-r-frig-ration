import CoolProp.CoolProp as CP

# +
#Q1.1


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
# +
#Q1.2
#Calcul de P_5s : 
print(CP.PropsSI('P', 'T', -8+273.15, 'Q', 1, 'CO2'))

#Calcul de T_2f :
print(294.82761632902583 + (CP.PropsSI('C', 'T', 294.82761632902583, 'P', 101300, 'Air')/CP.PropsSI('CP0MASS', 'T', 294.82761632902583, 'Q', 1, 'CO2')) * (20+273.15 - (294.82761632902583-5)))

#Calcul de q_c
print(-CP.PropsSI('C', 'T', 294.82761632902583, 'P', 101300, 'Air')*5)
# -



