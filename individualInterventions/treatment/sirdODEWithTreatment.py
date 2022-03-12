import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Set parameters
N = 1000
I0 = 1
gamma = -np.log(1 - 0.3)
mu = -np.log(1 - 0.1)
gammaHat = -np.log(1 - 0.6)
muHat = -np.log(1 - 0.05)
n = 7
etaA = 2/3
etaB = 2/3
tauA = 2/3
tauB = 2/3

maxT = 80
treatmentsAvailableVals = np.linspace(0, 1000, 101)

#Define functions determining how many treatments should be given to each group based on the total number of currently available treatments
def fA(currentlyAvailable):
    return currentlyAvailable/2

def fB(currentlyAvailable):
    return currentlyAvailable/2

#Calculate beta parameters
betaAA = etaA * tauA * n
betaAB = etaA * tauB * n
betaBA = etaB * tauA * n
betaBB = etaB * tauB * n

#Create list for storing final fatality numbers for each number of available treatments
finalDs = []

for treatmentsAvailable in treatmentsAvailableVals:
    
    #Set up equations (where X = [SA, SB, IA, IB, TA, TB, R, D])
    def dX_dt(X, t):
        return [-betaAA*X[2]*X[0]/N - betaBA*X[3]*X[0]/N,
                -betaAB*X[2]*X[1]/N - betaBB*X[3]*X[1]/N,
                betaAA*X[2]*X[0]/N + betaBA*X[3]*X[0]/N - np.min([fA(treatmentsAvailable - (X[4] + X[5])), X[2]]) - np.min([X[2] - np.min([fA(treatmentsAvailable - (X[4] + X[5])), X[2]]), fB(treatmentsAvailable - (X[4] + X[5])) - np.min([fB(treatmentsAvailable - (X[4] + X[5])), X[3]])]) - gamma*X[2] - mu*X[2],
                betaAB*X[2]*X[1]/N + betaBB*X[3]*X[1]/N - np.min([fB(treatmentsAvailable - (X[4] + X[5])), X[3]]) - np.min([X[3] - np.min([fB(treatmentsAvailable - (X[4] + X[5])), X[3]]), fA(treatmentsAvailable - (X[4] + X[5])) - np.min([fA(treatmentsAvailable - (X[4] + X[5])), X[2]])]) - gamma*X[3] - mu*X[3],
                np.min([fA(treatmentsAvailable - (X[4] + X[5])), X[2]]) + np.min([X[2] - np.min([fA(treatmentsAvailable - (X[4] + X[5])), X[2]]), fB(treatmentsAvailable - (X[4] + X[5])) - np.min([fB(treatmentsAvailable - (X[4] + X[5])), X[3]])]) - gammaHat*X[4] - muHat*X[4],
                np.min([fB(treatmentsAvailable - (X[4] + X[5])), X[3]]) - np.min([X[3] - np.min([fB(treatmentsAvailable - (X[4] + X[5])), X[3]]), fA(treatmentsAvailable - (X[4] + X[5])) - np.min([fA(treatmentsAvailable - (X[4] + X[5])), X[2]])]) - gammaHat*X[5] - muHat*X[5],
                gamma*X[2] + gamma*X[3] + gammaHat*X[4] + gammaHat*X[5],
                mu*X[2] + mu*X[3] + muHat*X[4] + muHat*X[5]]
    
    #Solve equations over specified time scale and with specified initial conditions
    plotT = range(0, maxT + 1)
    X0 = [(N - I0)/2, (N - I0)/2, I0/2, I0/2, 0, 0, 0, 0]
    sol = odeint(dX_dt, X0, plotT)
    S = sol[:,0] + sol[:,1]
    I = sol[:,2] + sol[:,3]
    T = sol[:,4] + sol[:,5]
    R = sol[:,6]
    D = sol[:,7]

    #Check that final number of deceased appears to be converging. If it is not, we should simulate for longer.
    if (abs(D[maxT] - D[maxT - 1]) > 0.1):
        print("WARNING: Possible non-convergence of fatality numbers when T-Hat = " + str(treatmentsAvailable) + ". Consider increasing maxT.")

    #Add final fatlity numbers to list
    finalDs.append(D[maxT])
    print("Completed simulation for T-Hat = " + str(treatmentsAvailable))

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()
#plt.title('Change in Final Fatality Numbers with v in SIRD ODE model')
#Plot final fatality numbers against v
plt.plot(treatmentsAvailableVals, finalDs, 'k-', linewidth = 2.5)
plt.xlabel(r'$\hat{T}$')
plt.ylabel('Final number of deceased individuals')
plt.xlim([0,1000])
plt.ylim([0, 300])

plt.show()
