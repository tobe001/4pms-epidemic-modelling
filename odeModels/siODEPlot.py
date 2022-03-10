import numpy as np
import matplotlib.pyplot as plt

#Set parameters for two different plots
beta1 = 0.5
beta2 = 2
#Population size
N = 1000
#Number of initial infectious cases
I0 = 1

#Compute I(t) and S(t) for t from 0 to 40 for both plots using analytically solved equation
plotT = np.linspace(0, 40, 82)
#Plot 1
I1 = N*I0/((N - I0)*np.exp(-beta1*plotT) + I0)
S1 = N - I1
#Plot 2
I2 = N*I0/((N - I0)*np.exp(-beta2*plotT) + I0)
S2 = N - I2

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()
plt.suptitle('Number of infectious and susceptible individuals over time in SI ODE model')
#Plot 1
plt.subplot(1, 2, 1)
plt.plot(plotT, I1, 'r-', plotT, S1, 'b-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible'])
plt.text(30, 700, r'$\beta = $' + str(beta1))
#Plot 2
plt.subplot(1, 2, 2)
plt.plot(plotT, I2, 'r-', plotT, S2, 'b-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible'])
plt.text(30, 700, r'$\beta = $' + str(beta2))

plt.show()
