import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Set parameters
N = 1000
I0 = 1
gamma = -np.log(1 - 0.3)
mu = -np.log(1 - 0.1)
n = 7
etaA = 2/3
etaB = 2/3
tauA = 1/2
tauB = 5/6

maxT = 80
vs = np.linspace(0, 999, 1000)

#Calculate beta parameters
betaAA = etaA * tauA * n
betaAB = etaA * tauB * n
betaBA = etaB * tauA * n
betaBB = etaB * tauB * n

#Create list for storing final fatality numbers for each v
finalDs = []

for v in vs:

    #Set up equations (where X = [SA, SB, IA, IB, R, D])
    def dX_dt(X, t):
        return [-betaAA*X[2]*X[0]/N - betaBA*X[3]*X[0]/N, -betaAB*X[2]*X[1]/N - betaBB*X[3]*X[1]/N, betaAA*X[2]*X[0]/N + betaBA*X[3]*X[0]/N - gamma*X[2] - mu*X[2], betaAB*X[2]*X[1]/N + betaBB*X[3]*X[1]/N - gamma*X[3] - mu*X[3], gamma*X[2] + gamma*X[3], mu*X[2] + mu*X[3]]

    #Calculate number of vaccines to be given
    vaccinesToB = np.min([v, (N - I0)/2])
    vaccinesToA = v - vaccinesToB

    #Calculate appropraite initial conditions
    SA0 = (N - I0)/2 - vaccinesToA
    SB0 = (N - I0)/2 - vaccinesToB
    R0 = vaccinesToA + vaccinesToB
    
    #Solve equations over specified time scale and with specified initial conditions
    plotT = range(0, maxT + 1)
    X0 = [SA0, SB0, I0/2, I0/2, R0, 0]
    sol = odeint(dX_dt, X0, plotT)
    S = sol[:,0] + sol[:,1]
    I = sol[:,2] + sol[:,3]
    R = sol[:,4]
    D = sol[:,5]

    #Check that final number of deceased appears to be converging. If it is not, we should simulate for longer.
    if (abs(D[maxT] - D[maxT - 1]) > 0.1):
        print("WARNING: Possible non-convergence of fatality numbers when v = " + str(v) + ". Consider increasing maxT.")

    #Add final fatlity numbers to list
    finalDs.append(D[maxT])
    print("Completed simulation for v = " + str(v))

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()
#plt.title('Change in Final Fatality Numbers with v in SIRD ODE model')
#Plot final fatality numbers against v
plt.plot(vs, finalDs, 'k-', linewidth = 2.5)
plt.xlabel('v')
plt.ylabel('Final number of deceased individuals')
plt.xlim([0,1000])
plt.ylim([0, 300])

plt.show()
