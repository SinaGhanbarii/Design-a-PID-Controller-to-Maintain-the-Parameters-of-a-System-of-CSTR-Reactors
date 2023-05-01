# Process Control Final project
# By: Aidin HosseiniNasab - MohammadSina GhanbariPakdehi
import numpy as np
import scipy.optimize as op
import math
import matplotlib.pyplot as plt
q = 100 #lit/min
CAf = 1 #mol / lit
Tf = 350 #K
Tcf=350 #K
V1=100 #lit
V2=100 #lit
hA1= 1.67e5 #j/min*K
hA2 = 1.67e5 #j/min*K
Cp=0.239 #j/g*K
Cpc=0.239 #j/g*K
deltaH = -4.78e4 # j/mol
rho = 1000 #g/lit
rhoC = 1000 #g/lit
EdivR = 1e4 #K
k0 = 7.2e10 #1/min
Qcmax = 500 # lit/min
Ca2_ss=0.005 #mol/lit

num=150  
Tpoint = np.arange(350,350+num,1) 
T2point = np.array([])
CaPoint = np.array([])
QcPoint = np.array([])
nTPoint = np.array([])

#iterating T1 [350,350+num);
for T1 in Tpoint:
    def eq1(x):
        return (q/V1)*(CAf-x)-k0*math.exp((-1*EdivR)/T1)*x
    
    Ca1= op.root(eq1,0.1, method='lm').x ##CA1 would be calculated.
    CaPoint = np.append(CaPoint,Ca1)

    def eq3(T2):
        return (q/V1)*(Ca1-Ca2_ss)-k0*math.exp((-1*EdivR)/T2)*Ca2_ss

    T2=op.root(eq3,400.0, method='lm').x #T2 would be calculated!
    T2point = np.append(T2point,T2)
    

    def eq2(qc):
        return (q/V1)*(Tf-T1)-((k0*deltaH*Ca1)/(rho*Cp))*math.exp((-1*EdivR)/T1)+((rhoC*Cpc*qc)/(rho*Cp*V1))*(1-math.exp((-1*hA1)/(rhoC*Cpc*qc)))*(Tcf-T1)

    qc=op.root(eq2,100.0, method='lm').x #qc would be calculated!
    QcPoint = np.append(QcPoint,qc)

    def eq4(nT1):
        return (q/V2)*(nT1-T2)-((k0*deltaH*Ca2_ss)/(rho*Cp))*math.exp((-1*EdivR)/T2)+((rhoC*Cpc*qc)/(rho*Cp*V2))*(1-math.exp((-1*hA2)/(rhoC*Cpc*qc)))*(nT1-T2+math.exp((-1*hA1)/(rhoC*Cpc*qc))*(Tcf-nT1))

    nT1 = op.root(eq4,400.0, method='lm').x # new T1 would be calculated!
    nTPoint = np.append(nTPoint,nT1)
    if T1 == 442:
        print(T1,' ',Ca1,' ',T2,' ',qc,' ',nT1)
    if T1 == 378:
        print(T1,' ',Ca1,' ',T2,' ',qc,' ',nT1)
    if T1 == 358:
        print(T1,' ',Ca1,' ',T2,' ',qc,' ',nT1)

    


#plotting the results
plt.plot(Tpoint,nTPoint,color = 'r',label='new T1 vs T1')
plt.plot(Tpoint,Tpoint,color = 'g' , label='x=y')
plt.xlabel('T1(K)')
plt.ylabel('new T1(K)')
plt.grid(True)
plt.xticks(np.arange(min(Tpoint), max(Tpoint)+1, 3))
plt.legend()
plt.show()



