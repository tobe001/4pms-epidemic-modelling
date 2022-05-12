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
tauA = 2/3
tauB = 2/3

maxT = 80
vs = np.linspace(0, 500, 51)

#Define functions to calculate number of vaccines to be intended for groups A and B as a function of the total number available and the number of susceptibles in each compartment
def fA(currentlyAvailable, SA, SB):
    return (SA / (SA + SB)) * currentlyAvailable

def fB(currentlyAvailable, SA, SB):
    return (SB / (SA + SB)) * currentlyAvailable

#Define functions to calculate number of vaccines provided to each group from the number intended for each one as a function of the total number available and the number of susceptibles in each group
def numVaccinesAA(currentlyAvailable, SA, SB):
    return np.min([fA(currentlyAvailable, SA, SB), SA])

def numVaccinesBB(currentlyAvailable, SA, SB):
    return np.min([fB(currentlyAvailable, SA, SB), SB])

def numVaccinesBA(currentlyAvailable, SA, SB):
    return np.min([SA - numVaccinesAA(currentlyAvailable, SA, SB), fB(currentlyAvailable, SA, SB) - numVaccinesBB(currentlyAvailable, SA, SB)])

def numVaccinesAB(currentlyAvailable, SA, SB):
    return np.min([SB - numVaccinesBB(currentlyAvailable, SA, SB), fA(currentlyAvailable, SA, SB) - numVaccinesAA(currentlyAvailable, SA, SB)])

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
        return [-betaAA*X[2]*X[0]/N - betaBA*X[3]*X[0]/N - numVaccinesAA(v, X[0], X[1]) - numVaccinesBA(v, X[0], X[1]),
                -betaAB*X[2]*X[1]/N - betaBB*X[3]*X[1]/N - numVaccinesBB(v, X[0], X[1]) - numVaccinesAB(v, X[0], X[1]),
                betaAA*X[2]*X[0]/N + betaBA*X[3]*X[0]/N - gamma*X[2] - mu*X[2],
                betaAB*X[2]*X[1]/N + betaBB*X[3]*X[1]/N - gamma*X[3] - mu*X[3],
                gamma*X[2] + gamma*X[3] + numVaccinesAA(v, X[0], X[1]) + numVaccinesBA(v, X[0], X[1]) + numVaccinesBB(v, X[0], X[1]) + numVaccinesAB(v, X[0], X[1]),
                mu*X[2] + mu*X[3]]

    #Specify initial conditions
    SA0 = (N - I0)/2
    SB0 = (N - I0)/2
    IA0 = I0/2
    IB0 = I0/2
    R0 = 0
    D0 = 0
    X0 = [SA0, SB0, IA0, IB0, R0, D0]
    
    #Solve equations over specified time scale and with specified initial conditions
    plotT = range(0, maxT + 1)
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
#Plot final fatality numbers against v
plt.plot(vs, finalDs, 'k-', linewidth = 2.5)
plt.xlabel(r'$\tilde{v}$')
plt.ylabel('Final number of deceased individuals')
plt.xlim([0,500])
plt.ylim([0, 300])

plt.show()
