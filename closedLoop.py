# Process Control Final project
# By: Aidin HosseiniNasab - MohammadSina GhanbariPakdehi
import scipy.optimize as op
import numpy as np
from math import exp
import scipy.integrate
import matplotlib.pyplot as plt



# proportional controller
def proportional(PV,SP,Kp=22654):
    return Kp*(SP-PV)

def PID(Kp, Ki, Kd, MV_bar=0):
   
    t_prev = 0
    I = 0
    
    # initial control
    MV = MV_bar
    
    while True:
        # yield MV, wait for new t, PV, SP
        t, PV, SP = yield MV
        
        # PID calculations
        e = SP - PV
        
        P = Kp*e
        I = I + Ki*e*(t - t_prev)
        D = 0
       
        MV = MV_bar + P + I + D
        
        # update stored data for next iteration
        e_prev = e
        t_prev = t


# model ode's are in this function
def model (t,z,K,PT):

    CA1, T1, CA2,T2 = z
    u= [100,1,350,98.99]
    q ,CAf,Tf,qc= u
    
    Tcf=350
    V1=100
    V2=100
    hA1= 1.67e5
    hA2 = 1.67e5
    Cp=0.239
    Cpc=0.239
    deltaH = -4.78e4
    rho = 1000
    rhoc = 1000
    EdivR = 1e4
    k0 = 7.2e10

    CA2setPoint = 0.005
    qc_ss = 98.99

    if t>=PT:
        Tf*=K

######################################## ATTENTION ###################################################    
#If you want to activate PID mode, clear "#"" symbol in the beginning of line 69 & make lin 70 as comment!

    qc = proportional(CA2,CA2setPoint)+qc_ss
    #qc=controller.send([t,CA2,CA2setPoint])

    if (qc>500): print("ERR")
    
    dCA1 = (q/V1)*(CAf-CA1)-k0*CA1*exp(-EdivR/T1)

    dT1=(q/V1)*(Tf-T1)-(k0*deltaH*CA1)/(rho*Cp)*exp(-EdivR/T1)+(rhoc*Cpc)/(rho*Cp*V1)*qc*(1-exp(-hA1/(rhoc*Cpc*qc)))*(Tcf-T1)

    dCA2=(q/V1)*(CA1-CA2)- k0*CA2*exp(-EdivR/T2)

    dT2 = (q/V2)*(T1-T2)-(k0*deltaH*CA2)/(rho*Cp)*exp(-EdivR/T2)+(rhoc*Cpc)/(rho*Cp*V2)*qc*(1-exp(-hA2/(rhoc*Cpc*qc)))*(T1-T2+exp(-hA1/(rhoc*Cpc*qc))*(Tcf-T1))

    return [dCA1,dT1,dCA2,dT2]

x0 = np.array([0.08506,442,0.005,449.9116]) #initial values
K=1.1  # disturbance
PT = 5 # u(t-PT)

#solves the 4 ode
#Kcu = 22654 , Pu=0.7337-> KI=
k=[7079.375,4385.85]## Tyreus-Luyben
#k=[22654,0] #to find Kcu
kp,ki=k
controller = PID(kp,ki,0,98.99) #define controller
controller.send(None)

ans = scipy.integrate.solve_ivp(model,[0,20],x0,args=(K,PT),dense_output=True,method='LSODA',rtol=1e-7,atol=1e-7)

t = np.linspace(0, 20, 3000)

z= ans.sol(t)

title_String = '+10% Tf closed loop'

#plotting the graphs
plt.plot(t,z[2])
plt.xlabel('Time (min)')
plt.ylabel('CA2 (mol/lit)')
plt.legend(['CA2'])
plt.title(title_String)
plt.grid('True')

plt.figure(2)
plt.plot(t,z[0])
plt.xlabel('Time (min)')
plt.ylabel('CA1 (mol/lit)')
plt.legend(['CA1'])
plt.title(title_String)
plt.grid('True')

plt.figure(3)
plt.plot(t,z[1])
plt.xlabel('Time (min)')
plt.ylabel('T1 (K)')
plt.legend(['T1'])
plt.title(title_String)
plt.grid('True')

plt.figure(4)
plt.plot(t,z[3])
plt.xlabel('Time (min)')
plt.ylabel('T2 (K)')
plt.legend(['T2'])
plt.title(title_String)
plt.grid('True')

plt.show()
