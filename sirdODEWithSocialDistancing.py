import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Set parameters
N = 1000
I0 = 1
gamma = -np.log(1 - 0.3)
mu = -np.log(1 - 0.1)
etaA = 2/3
etaB = 2/3
tauA = 2/3
tauB = 2/3

maxT = 80
ns = np.linspace(0, 10, 21)

#Create list for storing final fatality numbers for each n
finalDs = []

for n in ns:

    #Calculate beta parameters
    betaAA = etaA * tauA * n
    betaAB = etaA * tauB * n
    betaBA = etaB * tauA * n
    betaBB = etaB * tauB * n

    #Set up equations (where X = [SA, SB, IA, IB, R, D])
    def dX_dt(X, t):
        return [-betaAA*X[2]*X[0]/N - betaBA*X[3]*X[0]/N, -betaAB*X[2]*X[1]/N - betaBB*X[3]*X[1]/N, betaAA*X[2]*X[0]/N + betaBA*X[3]*X[0]/N - gamma*X[2] - mu*X[2], betaAB*X[2]*X[1]/N + betaBB*X[3]*X[1]/N - gamma*X[3] - mu*X[3], gamma*X[2] + gamma*X[3], mu*X[2] + mu*X[3]]

    #Solve equations over specified time scale and with specified initial conditions
    plotT = range(0, maxT + 1)
    X0 = [(N - I0)/2, (N - I0)/2, I0/2, I0/2, 0, 0]
    sol = odeint(dX_dt, X0, plotT)
    S = sol[:,0] + sol[:,1]
    I = sol[:,2] + sol[:,3]
    R = sol[:,4]
    D = sol[:,5]

    #Check that final number of deceased appears to be converging. If it is not, we should simulate for longer.
    if (abs(D[maxT] - D[maxT - 1]) > 0.1):
        print("WARNING: Possible non-convergence of fatality numbers when n = " + str(n) + ". Consider increasing maxT.")

    #Add final fatlity numbers to list
    finalDs.append(D[maxT])

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()
plt.title('Change in Final Fatality Numbers with n in SIRD ODE model')
#Plot final fatality numbers against n
plt.plot(ns, finalDs, 'k-', linewidth = 2.5)
plt.xlabel('n')
plt.ylabel('Final number of deceased individuals')
plt.ylim([0, 300])

plt.show()

