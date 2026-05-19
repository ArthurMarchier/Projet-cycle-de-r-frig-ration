import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt


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


# +
#Q1.2
def COP():
    Tcond,P3=iteration_T_cond(1e-3,1000)
    P5sv=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    P1=P5sv #On néglige la perte de charge dans les échangeurs
    h1=CP.PropsSI('H','P',P1,'T',-3+273.15,'CO2')
    h4=CP.PropsSI('H','T',Tcond-2,'Q',0,'CO2')
    h5=h4 #détente isenthalpique
    qf=h1-h5

    s1=CP.PropsSI('S','P',P1,'T',-3+273.15,'CO2')
    P2=P3 #on néglige la perte de charges dans l'échangeur
    s2is=s1
    h2is=CP.PropsSI('H','S',s2is,'P',P2,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    return qf/w

print(COP())


# +
def trace1():
    Tcond,P3=iteration_T_cond(1e-3,1000)
    T3sv=Tcond
    T3sl=Tcond
    T4=Tcond-2
    T5=-8+273.15
    T5sv=T5
    T1=T5+5

    P5sv=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    P1=P5sv #On néglige la perte de charge dans les échangeurs
    s1=CP.PropsSI('S','P',P1,'T',T1,'CO2')
    s3sv=CP.PropsSI('S','Q',1,'T',T3sv,'CO2')
    s3sl=CP.PropsSI('S','Q',0,'T',T3sl,'CO2')
    s4=CP.PropsSI('S','Q',0,'T',T4,'CO2')
    s5sv=CP.PropsSI('S','Q',1,'T',T5,'CO2')
    P5=P5sv #Pas de pertes de charge
    h4=CP.PropsSI('H','T',Tcond-2,'Q',0,'CO2')
    h5=h4 #détente isenthalpique
    s5=CP.PropsSI('S','H',h5,'P',P5,'CO2')

    #On traite le cas du point 2
    P2=P3 #on néglige la perte de charges dans l'échangeur
    s2is=s1
    h2is=CP.PropsSI('H','S',s2is,'P',P2,'CO2')
    h1=CP.PropsSI('H','P',P1,'T',T1,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    h2=h1+w
    s2=CP.PropsSI('S','H',h2,'P',P2,'CO2')
    T2=CP.PropsSI('T','H',h2,'P',P2,'CO2')

    s=[s1,s2,s3sv,s3sl,s4,s5,s5sv,s1]
    T=[T1,T2,T3sv,T3sl,T4,T5,T5sv,T1]

    #Pour l'air
    TA1,TA2=15+273.15,20+273.15
    p_atm=101300
    sA1=CP.PropsSI('S','T',TA1,'P',p_atm,'Air')
    sA2=CP.PropsSI('S','T',TA2,'P',p_atm,'Air')
    Tair=[TA1,TA2]
    sair=[sA1,sA2]

    plt.plot(s, T, marker='o', color='green', linewidth=2, markersize=8, label='CO₂')
    plt.plot(sair, Tair, marker='o', color='red', linewidth=2, markersize=8, label='Air')

    for i, (s_val, T_val) in enumerate(zip(s, T)):
        plt.scatter(s_val, T_val, color='green')
        if i in [0,1]:
            plt.text(s_val, T_val, f' {i+1}', fontsize=12, ha='right', va='bottom', color='green')
        elif i==2:
            plt.text(s_val, T_val, f'3sv', fontsize=12, ha='right', va='bottom', color='green')
        elif i==3:
            plt.text(s_val, T_val, f'3sl', fontsize=12, ha='right', va='bottom', color='green')
        elif i in [4,5]:
            plt.text(s_val, T_val, f' {i}', fontsize=12, ha='right', va='bottom', color='green')
        elif i==6:
            plt.text(s_val, T_val, f'5sv', fontsize=12, ha='right', va='bottom', color='green')
            
    plt.scatter(sA1, TA1, color='red')
    plt.text(sA1, TA1, f'A1', fontsize=12, ha='right', va='bottom', color='red')
    plt.scatter(sA2, TA2, color='red')
    plt.text(sA2, TA2, f'A2', fontsize=12, ha='right', va='bottom', color='red')
    
    plt.xlabel('s (J/kg·K)', fontsize=12)
    plt.ylabel('T (K)', fontsize=12)
    plt.title('Diagramme T-s du cycle simple sous-critique', fontsize=14)
    plt.legend(fontsize=10)  # Affiche la légende "Cycle CO₂"
    plt.grid(True, alpha=0.3)  # Grille légère
    
trace1()


# +
def trace1log():
    Tcond,P3=iteration_T_cond(1e-3,1000)
    T3sv=Tcond
    T3sl=Tcond
    T4=Tcond-2
    T5=-8+273.15
    T5sv=T5
    T1=T5+5

    P5sv=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    P1=P5sv #On néglige la perte de charge dans les échangeurs
    s1=CP.PropsSI('S','P',P1,'T',T1,'CO2')
    s3sv=CP.PropsSI('S','Q',1,'T',T3sv,'CO2')
    s3sl=CP.PropsSI('S','Q',0,'T',T3sl,'CO2')
    s4=CP.PropsSI('S','Q',0,'T',T4,'CO2')
    s5sv=CP.PropsSI('S','Q',1,'T',T5,'CO2')
    P5=P5sv #Pas de pertes de charge
    h4=CP.PropsSI('H','T',Tcond-2,'Q',0,'CO2')
    h5=h4 #détente isenthalpique
    s5=CP.PropsSI('S','H',h5,'P',P5,'CO2')

    #On traite le cas du point 2
    P2=P3 #on néglige la perte de charges dans l'échangeur
    s2is=s1
    h2is=CP.PropsSI('H','S',s2is,'P',P2,'CO2')
    h1=CP.PropsSI('H','P',P1,'T',T1,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    h2=h1+w
    s2=CP.PropsSI('S','H',h2,'P',P2,'CO2')
    T2=CP.PropsSI('T','H',h2,'P',P2,'CO2')

    s=[s1,s2,s3sv,s3sl,s4,s5,s5sv,s1]
    T=[T1,T2,T3sv,T3sl,T4,T5,T5sv,T1]

    #Pour l'air
    TA1,TA2=15+273.15,20+273.15
    p_atm=101300
    sA1=CP.PropsSI('S','T',TA1,'P',p_atm,'Air')
    sA2=CP.PropsSI('S','T',TA2,'P',p_atm,'Air')
    Tair=[TA1,TA2]
    sair=[sA1,sA2]

    #Courbe de saturation
    Tsat=np.linspace(0,31+273.15,1000) #La température critique du CO2 est de 31°C
    sliq=CP.PropsSI('S','Q',0,'T',Tsat,'CO2')
    sgaz=CP.PropsSI('S','Q',1,'T',Tsat,'CO2')

    plt.xlim(1e3,4e3)
    plt.ylim(250,350)
    plt.xscale('log')
    plt.plot(s, T, marker='o', color='green', linewidth=2, markersize=8, label='CO₂')
    plt.plot(sair, Tair, marker='o', color='red', linewidth=2, markersize=8, label='Air')
    plt.plot(sliq, Tsat, color='black', linewidth=0.5)
    plt.plot(sgaz, Tsat, color='black', linewidth=0.5)

    for i, (s_val, T_val) in enumerate(zip(s, T)):
        plt.scatter(s_val, T_val, color='green')
        if i in [0,1]:
            plt.text(s_val, T_val, f' {i+1}', fontsize=12, ha='right', va='bottom', color='green')
        elif i==2:
            plt.text(s_val, T_val, f'3sv', fontsize=12, ha='right', va='bottom', color='green')
        elif i==3:
            plt.text(s_val, T_val, f'3sl', fontsize=12, ha='right', va='bottom', color='green')
        elif i in [4,5]:
            plt.text(s_val, T_val, f' {i}', fontsize=12, ha='right', va='bottom', color='green')
        elif i==6:
            plt.text(s_val, T_val, f'5sv', fontsize=12, ha='right', va='bottom', color='green')
            
    plt.scatter(sA1, TA1, color='red')
    plt.text(sA1, TA1, f'A1', fontsize=12, ha='right', va='bottom', color='red')
    plt.scatter(sA2, TA2, color='red')
    plt.text(sA2, TA2, f'A2', fontsize=12, ha='right', va='bottom', color='red')
    
    plt.xlabel('s (J/kg·K)', fontsize=12)
    plt.ylabel('T (K)', fontsize=12)
    plt.title('Diagramme T-s du cycle simple sous-critique avec abscisses en échelle log', fontsize=14)
    plt.legend(fontsize=10)  # Affiche la légende "Cycle CO₂"
    plt.grid(True, alpha=0.3)  # Grille légère
    
trace1log()


# +
#Q2.1
def COP2(P3):
    T3=35+273.15
    h1=CP.PropsSI('H','T',-8+273.15,'Q',1,'CO2')
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    h5=h3 #détente isenthalpique
    qf=h1-h5

    P2=P3 #Pas de perte de charges dans le refroidisseur
    P1=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    s1=CP.PropsSI('S','T',-8+273.15,'Q',1,'CO2')
    s2is=s1
    h2is=CP.PropsSI('H','P',P2,'S',s2is,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    return qf/w
    
print(COP2(120e5))


# +
def maximiseur_COP():
    dP = 10e5
    P_opt = 120e5
    COP_opt = COP2(P_opt)

    while dP > 1e2:
        P_test = P_opt + dP
        COP_test = COP2(P_test)
        
        if COP_test > COP_opt:
            P_opt = P_test
            COP_opt = COP_test
            
        else:
            P_test = P_opt - dP
            COP_test = COP2(P_test)
            
            if COP_test > COP_opt:
                P_opt = P_test
                COP_opt = COP_test
                
            else: #Cas où aucun des deux côtés n'améliore le COP
                dP /= 2

    return P_opt

maximiseur_COP()


# +
def trace2():
    P3=maximiseur_COP()
    T3=35+273.15
    T1=-8+273.15
    T5=T1
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    h4=h3 #détente isenthalpique
    #T4=CP.PropsSI('T','H',h4,'Q',0,'CO2') --> nécessaire de connaître P4 pour calculer

    s1=CP.PropsSI('S','T',T1,'Q',1,'CO2')
    P1=CP.PropsSI('P','T',T1,'Q',1,'CO2')
    s3=CP.PropsSI('S','T',T3,'P',P3,'CO2')
    #s4=CP.PropsSI('S','H',h4,'Q',0,'CO2')
    P5=P1 #Pas de perte de charges dans l'échangeur
    h5=h3 #détente isenthalpique
    s5=CP.PropsSI('S','P',P5,'H',h5,'CO2')

    #Cas du point 2
    P2=P3 #Pas de perte de charges dans le refroidisseur
    s2is=s1
    h2is=CP.PropsSI('H','P',P2,'S',s2is,'CO2')
    h1=CP.PropsSI('H','P',P1,'S',s1,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    h2=h1+w
    s2=CP.PropsSI('S','P',P2,'H',h2,'CO2')
    T2=CP.PropsSI('T','P',P2,'H',h2,'CO2')

    s=[s1,s2,s3,s5,s1]
    T=[T1,T2,T3,T5,T1]
    
    plt.plot(s, T, marker='o', color='green', linewidth=2, markersize=8, label='CO₂')
    for i, (s_val, T_val) in enumerate(zip(s, T)):
        if i==4:
            break
        plt.scatter(s_val, T_val, color='green')
        plt.text(s_val, T_val, f' {i+1}', fontsize=12, ha='right', va='bottom', color='green')

    plt.xlabel('s (J/kg·K)', fontsize=12)
    plt.ylabel('T (K)', fontsize=12)
    plt.title('Diagramme T-s du cycle transcritique simple', fontsize=14)
    plt.legend(fontsize=10)  # Affiche la légende "Cycle CO₂"
    plt.grid(True, alpha=0.3)  # Grille légère
        
trace2()


# +
def trace2():
    P3=maximiseur_COP()
    T3=35+273.15
    T1=-8+273.15
    T5=T1
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    h4=h3 #détente isenthalpique
    #T4=CP.PropsSI('T','H',h4,'Q',0,'CO2') --> nécessaire de connaître P4 pour calculer

    s1=CP.PropsSI('S','T',T1,'Q',1,'CO2')
    P1=CP.PropsSI('P','T',T1,'Q',1,'CO2')
    s3=CP.PropsSI('S','T',T3,'P',P3,'CO2')
    #s4=CP.PropsSI('S','H',h4,'Q',0,'CO2')
    P5=P1 #Pas de perte de charges dans l'échangeur
    h5=h3 #détente isenthalpique
    s5=CP.PropsSI('S','P',P5,'H',h5,'CO2')

    #Cas du point 2
    P2=P3 #Pas de perte de charges dans le refroidisseur
    s2is=s1
    h2is=CP.PropsSI('H','P',P2,'S',s2is,'CO2')
    h1=CP.PropsSI('H','P',P1,'S',s1,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    h2=h1+w
    s2=CP.PropsSI('S','P',P2,'H',h2,'CO2')
    T2=CP.PropsSI('T','P',P2,'H',h2,'CO2')

    s=[s1,s2,s3,s5,s1]
    T=[T1,T2,T3,T5,T1]

    #Courbe de saturation
    Tsat=np.linspace(0,31+273.15,1000) #La température critique du CO2 est de 31°C
    sliq=CP.PropsSI('S','Q',0,'T',Tsat,'CO2')
    sgaz=CP.PropsSI('S','Q',1,'T',Tsat,'CO2')

    plt.xlim(1e3,2.1e3)
    plt.ylim(250,400)
    plt.xscale('log')
    plt.plot(s, T, marker='o', color='green', linewidth=2, markersize=8, label='CO₂')
    plt.plot(sliq, Tsat, color='black', linewidth=0.5)
    plt.plot(sgaz, Tsat, color='black', linewidth=0.5)
    for i, (s_val, T_val) in enumerate(zip(s, T)):
        if i==4:
            break
        plt.scatter(s_val, T_val, color='green')
        plt.text(s_val, T_val, f' {i+1}', fontsize=12, ha='right', va='bottom', color='green')

    plt.xlabel('s (J/kg·K)', fontsize=12)
    plt.ylabel('T (K)', fontsize=12)
    plt.title('Diagramme T-s du cycle transcritique simple avec abscisses en échelle log', fontsize=14)
    plt.legend(fontsize=10)  # Affiche la légende "Cycle CO₂"
    plt.grid(True, alpha=0.3)  # Grille légère
        
trace2()
# -


