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
    while(k < n_iterations and abs(Tk1 - Tk) > ecart_min):
        Tk=Tk1
        cp_air = CP.PropsSI('C', 'T', Tk-5, 'P', p_atm, 'Air')
        cp_co2_liq = CP.PropsSI('CP0MASS', 'T', Tk, 'Q', 0, 'CO2')
        L_v_co2 = CP.PropsSI('L', 'T', Tk, 'Q', 0, 'CO2')
        Tk1= 15 + 273.15 + 5 + ((L_v_co2 + 2*cp_co2_liq)/cp_air)
        k+=1
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
    T1 = -8 + 273.15 + 5
    P1 = CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    h1 = CP.PropsSI('H','T',T1,'P',P1,'CO2')
    s1 = CP.PropsSI('S','T',T1,'P',P1,'CO2')
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    h5=h3 #détente isenthalpique
    qf=h1-h5
    
    P2=P3 #Pas de perte de charges dans le refroidisseur
    s2is=s1
    h2is=CP.PropsSI('H','P',P2,'S',s2is,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    return qf/w
    
print(COP2(120e5))


# +
def maximiseur_COP():
    dP = 1e5
    P_opt = 90e5
    COP_opt = COP2(P_opt)

    while dP > 1e2:
        P_test = P_opt + dP
        COP_test = COP2(P_test)
        
        if COP_test > COP_opt:
            P_opt = P_test
            COP_opt = COP_test
            
        else:
            P_test = P_opt - dP
            if P_test < 70e5:
                P_test = 70e5
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
    T1=-8+273.15+5
    T5=-8+273.15
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    h4=h3 #détente isenthalpique
    #T4=CP.PropsSI('T','H',h4,'Q',0,'CO2') --> nécessaire de connaître P4 pour calculer

    P1=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    s1=CP.PropsSI('S','T',T1,'P',P1,'CO2')
    s3=CP.PropsSI('S','T',T3,'P',P3,'CO2')
    #s4=CP.PropsSI('S','H',h4,'Q',0,'CO2')
    P5=P1 #Pas de perte de charges dans l'échangeur
    h5=h3 #détente isenthalpique
    s5=CP.PropsSI('S','P',P5,'H',h5,'CO2')

    #Cas du point 2
    P2=P3 #Pas de perte de charges dans le refroidisseur
    s2is=s1
    h2is=CP.PropsSI('H','P',P2,'S',s2is,'CO2')
    h1=CP.PropsSI('H','T',T1,'P',P1,'CO2')
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
def trace2log():
    P3=maximiseur_COP()
    T3=35+273.15
    T1=-8+273.15+5
    T5=-8+273.15
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    h4=h3 #détente isenthalpique
    #T4=CP.PropsSI('T','H',h4,'Q',0,'CO2') --> nécessaire de connaître P4 pour calculer

    P1=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    s1=CP.PropsSI('S','T',T1,'P',P1,'CO2')
    s3=CP.PropsSI('S','T',T3,'P',P3,'CO2')
    #s4=CP.PropsSI('S','H',h4,'Q',0,'CO2')
    P5=P1 #Pas de perte de charges dans l'échangeur
    h5=h3 #détente isenthalpique
    s5=CP.PropsSI('S','P',P5,'H',h5,'CO2')

    #Cas du point 2
    P2=P3 #Pas de perte de charges dans le refroidisseur
    s2is=s1
    h2is=CP.PropsSI('H','P',P2,'S',s2is,'CO2')
    h1=CP.PropsSI('H','T',T1,'P',P1,'CO2')
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
        
trace2log()
# +
#Q2.2
def COP1bis(P3): #Cette fois P3 est une donnée du problème
    Tcond=CP.PropsSI('T','P',P3,'Q',1,'CO2')
    P5sv=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    P1=P5sv #On néglige la perte de charge dans les échangeurs
    T1 = -8+273.15+5
    h1=CP.PropsSI('H','T',T1,'P',P1,'CO2')
    s1=CP.PropsSI('S','T',T1,'P',P1,'CO2')

    h4=CP.PropsSI('H','T',Tcond-2,'Q',0,'CO2')
    h5=h4 #détente isenthalpique
    qf=h1-h5
    P2=P3
    s2is=s1
    h2is=CP.PropsSI('H','S',s2is,'P',P2,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    return qf/w


def trace_COP():
    P3_1 = np.linspace(4.5e6,7.3e6,100)   # domaine sous-critique
    P3_2 = np.linspace(7e6,15e6,100)  # domaine transcritique

    COP_1 = [COP1bis(p) for p in P3_1]
    COP_2 = [COP2(p) for p in P3_2]

    plt.plot(P3_1,COP_1,label='COP de la question 1.2')
    plt.plot(P3_2,COP_2,label='COP de la question 2.1')

    plt.xlabel('P (Pa)', fontsize=12)
    plt.ylabel('COP', fontsize=12)
    plt.title('Evolution du COP en fonction de la pression')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.show()

trace_COP()


# +
#Q3.1
def iteration_T_cond_bis(Tamb, ecart_min,n_iterations): #cette fois Tamb est une variable de la fonction
    Tk = Tamb+10
    p_atm= 101300
    cp_air = CP.PropsSI('C', 'T', Tk-5, 'P', p_atm, 'Air')
    cp_co2_liq = CP.PropsSI('CP0MASS', 'T', Tk, 'Q', 0, 'CO2')
    L_v_co2 = CP.PropsSI('L', 'T', Tk, 'Q', 0, 'CO2')
    Tk1= Tamb + 5 + ((L_v_co2 + 2*cp_co2_liq)/cp_air)
    k=0
    while(k<n_iterations and abs(Tk1-Tk)>ecart_min):
        Tk=Tk1
        cp_air = CP.PropsSI('C', 'T', Tk-5, 'P', p_atm, 'Air')
        cp_co2_liq = CP.PropsSI('CP0MASS', 'T', Tk, 'Q', 0, 'CO2')
        L_v_co2 = CP.PropsSI('L', 'T', Tk, 'Q', 0, 'CO2')
        Tk1= Tamb + 5 + ((L_v_co2 + 2*cp_co2_liq)/cp_air)
        k += 1
    p_sat = CP.PropsSI('P', 'T', Tk1, 'Q', 1, 'CO2')
    return Tk1, p_sat

 
def COP1ter(Tamb):
    Tcond,P3=iteration_T_cond_bis(Tamb,1e-3,1000)
    P5sv=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    P1=P5sv #On néglige la perte de charge dans les échangeurs
    h1=CP.PropsSI('H','P',P1,'T',-8+273.15+5,'CO2')
    h4=CP.PropsSI('H','T',Tcond-2,'Q',0,'CO2')
    h5=h4 #détente isenthalpique
    qf=h1-h5

    s1=CP.PropsSI('S','P',P1,'T',-8+273.15+5,'CO2')
    P2=P3 #on néglige la perte de charges dans l'échangeur
    s2is=s1
    h2is=CP.PropsSI('H','S',s2is,'P',P2,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    return qf/w


def COP2bis(P3,Tamb): #Tamb devient une variable
    T3=Tamb+5
    P1=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    h1=CP.PropsSI('H','T',-8+273.15+5,'P',P1,'CO2')
    s1=CP.PropsSI('S','T',-8+273.15+5,'P',P1,'CO2')
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    h5=h3 #détente isenthalpique
    qf=h1-h5
    P2=P3
    s2is=s1
    h2is=CP.PropsSI('H','P',P2,'S',s2is,'CO2')
    wis=h2is-h1
    w=wis/(1-0.121*(P2/P1))
    return qf/w


def maximiseur_COPbis(Tamb):
    dP = 10e5
    P_opt = 120e5
    COP_opt = COP2bis(P_opt,Tamb)

    while dP > 1e2:
        P_test = P_opt - dP
        if P_test < 70e5:   #limite basse
            P_test = 70e5
        COP_test = COP2bis(P_test,Tamb)
        
        if COP_test > COP_opt:
            P_opt = P_test
            COP_opt = COP_test
            
        else:
            P_test = P_opt + dP
            COP_test = COP2bis(P_test,Tamb)
            
            if COP_test > COP_opt:
                P_opt = P_test
                COP_opt = COP_test
                
            else: #Cas où aucun des deux côtés n'améliore le COP
                dP /= 2
    return P_opt

def trace_max_COP():
    Tamb=[10,20,30,40]
    P30=maximiseur_COPbis(30+273.15)
    P40=maximiseur_COPbis(40+273.15)
    T10,P10=iteration_T_cond_bis(10+273.15,1e-3,1000)
    T20,P20=iteration_T_cond_bis(20+273.15,1e-3,1000)
    COP=[COP1ter(10+273.15),COP1ter(20+273.15),COP2bis(P30,30+273.15),COP2bis(P40,40+273.15)]

    plt.plot(Tamb,COP,marker='o', label='COP max')
    plt.xlabel('Tamb (°C)', fontsize=12)
    plt.ylabel('COP maximale', fontsize=12)
    plt.title('COP maximale en fonction de la température ambiante', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.show()

    P_bar=np.array([P10,P20,P30,P40])*1e-5
    plt.plot(Tamb,P_bar,marker='o', label='P3 optimale')
    for i, (T, P) in enumerate(zip(Tamb, P_bar)):
        plt.text(T, P, f'{P:.1f} bar', fontsize=10, ha='left', va='bottom', color='red')
    plt.xlabel('Tamb (°C)', fontsize=12)
    plt.ylabel('P3 (bar)', fontsize=12)
    plt.title('P3 en fonction de la température ambiante', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.show()
    return P10,P20,P30,P40

P10,P20,P30,P40=trace_max_COP()
print("P10 =",P10*1e-5,'bar')
print("P20 =",P20*1e-5,'bar')
print("P30 =",P30*1e-5,'bar')
print("P40 =",P40*1e-5,'bar')


# +
#Q4.1

def COP_propane(Tevap, Tcond):
    COP_carnot = Tevap / (Tcond - Tevap)
    return 0.2 * COP_carnot

def COP4_souscritique(P3,SC):
    Tcond=CP.PropsSI('T','P',P3,'Q',1,'CO2')
    P5sv=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    P1=P5sv #On néglige la perte de charge dans les échangeurs
    T1 = -8+273.15+5
    h1 = CP.PropsSI('H','T',T1,'P',P1,'CO2')
    h4=CP.PropsSI('H','T',Tcond-2,'Q',0,'CO2')
    P4=P3 #On néglige la perte de charge dans les échangeurs
    P4sc=P4 #On néglige la perte de charge dans les échangeurs
    T4sc=Tcond-2-SC
    h4sc=CP.PropsSI('H','T',T4sc,'P',P4sc,'CO2')
    h5=h4sc #détente isenthalpique
    qf=h1-h5

    s1 = CP.PropsSI('S','T',T1,'P',P1,'CO2')
    P2=P3 #on néglige la perte de charges dans l'échangeur
    s2is=s1
    h2is=CP.PropsSI('H','S',s2is,'P',P2,'CO2')
    wis=h2is-h1
    tau=P2/P1
    eta_is=0.3774+0.14*tau-0.02*tau**2+0.001*tau**3
    w=wis/eta_is

    # DMSS
    h3=CP.PropsSI('H','T',Tcond,'Q',0,'CO2')
    Q_DMSS=h3-h4sc
    Tevap_prop=T4sc-5
    Tcond_prop=Tcond
    COP_prop=COP_propane(Tevap_prop,Tcond_prop)
    w_DMSS=Q_DMSS/COP_prop
    
    return qf/(w+w_DMSS)

def COP4_transcritique(P3,Tamb,SC):
    T3=Tamb+5
    P4=P3
    T4=T3-SC
    h4=CP.PropsSI('H','T',T4,'P',P4,'CO2')
    h5=h4
    P1=CP.PropsSI('P','T',-8+273.15,'Q',1,'CO2')
    T1=-8+273.15+5
    h1=CP.PropsSI('H','T',T1,'P',P1,'CO2')
    qf=h1-h5

    s1=CP.PropsSI('S','T',T1,'P',P1,'CO2')
    P2=P3
    s2is=s1
    h2is=CP.PropsSI('H','S',s2is,'P',P2,'CO2')
    wis=h2is-h1
    tau=P2/P1
    eta_is=0.3774+0.14*tau-0.02*tau**2+0.001*tau**3
    w=wis/eta_is

    # DMSS
    h3=CP.PropsSI('H','T',T3,'P',P3,'CO2')
    Q_DMSS=h3-h4
    Tevap_prop=T4-5
    Tcond_prop=Tamb+5    
    COP_prop=COP_propane(Tevap_prop,Tcond_prop)
    w_DMSS=Q_DMSS/COP_prop
    
    return qf/(w+w_DMSS)

print(COP4_transcritique(120e5,30+273.15,8))
print(COP2bis(120e5,30+273.15))
# +
#Q4.2
def trace4_transcritique(Tamb):
    P3 = np.linspace(70e5, 140e5, 200)
    
    COP_avecSC_transcritique = [COP4_transcritique(p, Tamb, 5) for p in P3]
    COP_sansSC_transcritique = [COP2bis(p, Tamb) for p in P3]

    plt.plot(P3*1e-5, COP_avecSC_transcritique, label='COP avec DMSS')
    plt.plot(P3*1e-5, COP_sansSC_transcritique, label='COP sans DMSS')

    plt.xlabel('P (bar)', fontsize=12)
    plt.ylabel('COP', fontsize=12)
    plt.title('Influence de la pression de refoulement sur le COP (régime transcritique)', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.show()

def trace4_souscritique(Tamb):
    P3 = np.linspace(45e5, 70e5, 200)

    COP_avecSC_souscritique = [COP4_souscritique(p, 5) for p in P3]
    COP_sansSC_souscritique = [COP1bis(p) for p in P3]

    plt.plot(P3*1e-5, COP_avecSC_souscritique, label='COP avec DMSS')
    plt.plot(P3*1e-5, COP_sansSC_souscritique, label='COP sans DMSS')

    plt.xlabel('P (bar)', fontsize=12)
    plt.ylabel('COP', fontsize=12)
    plt.title('Influence de la pression de refoulement sur le COP (régime sous-critique)', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.show()

trace4_transcritique(30+273.15)
trace4_souscritique(20+273.15)


# -

#Q4.3
import time

# +
def find_optimal_couple_souscritique_DMSS(Tamb, SC_range=(0, 20, 50)):
    start_time = time.time()
    T_cond,best_P3=iteration_T_cond_bis(Tamb,1e-3,1000)
    SC_min, SC_max, n_SC = SC_range

    best_COP = -np.inf
    best_SC = None

    for i in range(5):
        SC_values = np.linspace(SC_min, SC_max, n_SC)
        dSC=(SC_max-SC_min)/n_SC
        for SC in SC_values:
            current_COP = COP4_souscritique(best_P3, SC)
            if current_COP > best_COP:
                best_COP = current_COP
                best_SC = SC
                
        SC_min, SC_max = best_SC - dSC, best_SC + dSC

    # Temps écoulé
    time_elapsed = time.time() - start_time

    return best_P3, best_SC, time_elapsed


def find_optimal_couple_transcritique_DMSS(Tamb, P3_range=(4.5e6,12e6, 50), SC_range=(0, 20, 50)):
    start_time = time.time()

    P3_min, P3_max, n_P3 = P3_range
    SC_min, SC_max, n_SC = SC_range

    best_COP = -np.inf
    best_P3 = 0
    best_SC = 0

    for i in range(5):
        P3_values = np.linspace(P3_min, P3_max, n_P3)
        SC_values = np.linspace(SC_min, SC_max, n_SC)
        dP3=(P3_max-P3_min)/n_P3
        dSC=(SC_max-SC_min)/n_SC
        for P3 in P3_values:
            for SC in SC_values:
                current_COP = COP4_transcritique(P3,Tamb,SC)
                if current_COP > best_COP:
                    best_COP = current_COP
                    best_P3 = P3
                    best_SC = SC
        P3_min, P3_max = best_P3-dP3,best_P3+dP3
        SC_min,SC_max=best_SC-dSC,best_SC+dSC
        
    # Temps écoulé
    time_elapsed = time.time() - start_time

    return best_P3, best_SC, time_elapsed


# -

print(find_optimal_couple_souscritique_DMSS(283.15, SC_range=(0, 15, 50)))
print(find_optimal_couple_souscritique_DMSS(293.15, SC_range=(0, 15, 50)))
print(find_optimal_couple_transcritique_DMSS(303.15, P3_range=(4.5e6,12e6, 50), SC_range=(0, 15, 50)))
print(find_optimal_couple_transcritique_DMSS(313.15, P3_range=(4.5e6,12e6, 50), SC_range=(0, 15, 50)))


# +
#Q5

def find_optimal_souscritique_sans_DMSS(Tamb):
    Tcond, P3 = iteration_T_cond_bis(Tamb, 1e-3, 1000)
    return P3

def find_optimal_transcritique_sans_DMSS(Tamb, P3_range=(4.5e6,12e6,50)):

    P3_min, P3_max, n_P3 = P3_range
    P3_values = np.linspace(P3_min, P3_max, n_P3)

    best_COP = -np.inf
    best_P3 = 0

    for P3 in P3_values:
        current_COP = COP2bis(P3, Tamb)
        if current_COP > best_COP:
            best_COP = current_COP
            best_P3 = P3

    return best_P3


# +
Tamb_values = [283.15, 293.15, 303.15, 313.15]

for Tamb in Tamb_values:

    print("===== Tamb =", Tamb, "K =====")

    # ----- AVEC DMSS -----
    if Tamb < 303.15:
        P3_DMSS, SC_DMSS, _ = find_optimal_couple_souscritique_DMSS(Tamb)
        COP_DMSS = COP4_souscritique(P3_DMSS, SC_DMSS)

        P3_sans = find_optimal_souscritique_sans_DMSS(Tamb)
        COP_sans = COP1bis(P3_sans)

    else:
        P3_DMSS, SC_DMSS, _ = find_optimal_couple_transcritique_DMSS(Tamb)
        COP_DMSS = COP4_transcritique(P3_DMSS, Tamb, SC_DMSS)

        P3_sans = find_optimal_transcritique_sans_DMSS(Tamb)
        COP_sans = COP2bis(P3_sans, Tamb)

    print("Avec DMSS :", round(P3_DMSS*1e-5,2), "bar | COP =", round(COP_DMSS,3))
    print("Sans DMSS :", round(P3_sans*1e-5,2), "bar | COP =", round(COP_sans,3))
    print()
# -


