import numpy as np
import matplotlib.pyplot as plt

#Set parameters
N = 1000
I0 = 1
beta = 2
#Time scale over which to plot
plotT = np.linspace(0, 40, 82)

#Compute I(t) and S(t) for given t using analytically solved equation
I = N*I0/((N - I0)*np.exp(-beta*plotT) + I0)
S = N - I

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()

plt.plot(plotT, I, 'r-', plotT, S, 'b-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible'])
plt.text(30, 700, r'$\beta = $' + str(beta))

plt.show()
