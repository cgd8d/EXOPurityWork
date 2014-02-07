'''
Likelihood method
Correct the data; then fit it against a model function, ie:

P(E, drift_time) = (2*pi*sigma^2)^(-1/2) *
                   exp(-[E exp( drift_time /L) - E_0]^2/[2*sigma^2])

data is read in as an array of (charge_energy, drift_time) tuples.

To run this script, you should already have a csv file containing events.
Do: python <csv filename>
'''

import math
import scipy.optimize
import csv
import sys
import matplotlib.pyplot

def Prob(L, E, drift_time, sigma, E_0):
    purity_corr = math.exp(drift_time/L)
    ret_val = math.exp(-pow(E*purity_corr - E_0, 2)/(2*sigma*sigma))
    ret_val /= math.sqrt(2*math.pi)*sigma
    return ret_val

def NegLogLikelihood(x, EData, DriftTimeData):
    L = x[0]
    sigma = x[1]
    E_0 = x[2]
    ret_val = 0.
    if len(EData) != len(DriftTimeData):
        print "Energy and drift time data are mismatched."
        sys.exit(1)
    for i in range(len(EData)):
        ret_val -= math.log(Prob(L, EData[i], DriftTimeData[i], sigma, E_0))
    return ret_val

# Import data from csv file.  Store it in (charge_energy, drift_time) tuples.
ChargeEnergyData = []
DriftTimeData = []
with open(sys.argv[1], 'rb') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if float(row[1]) < 2200: continue # Skip for now, need to figure out why it crept in.
        ChargeEnergyData.append(float(row[1]))
        DriftTimeData.append(float(row[0]))

xopt, NLL_vals = scipy.optimize.fmin(NegLogLikelihood,
                                     (1000., 100., 2600.),
                                     args = (ChargeEnergyData, DriftTimeData),
                                     retall = True)
print "Result of optimization:"
print "Lifetime =", xopt[0], "us."
print "Charge sigma =", xopt[1], "keV."
print "Charge peak =", xopt[2], "keV."

# Graphics to show the result.
def PeakPositionFunc(drift_time):
    return xopt[2] * math.exp(-drift_time/xopt[0])
funcxpts = [0.1*x for x in range(1200)]
funcypts = [PeakPositionFunc(x) for x in funcxpts]
matplotlib.pyplot.plot(DriftTimeData, ChargeEnergyData, 'b.')
matplotlib.pyplot.plot(funcxpts, funcypts, 'r-')
matplotlib.pyplot.show()

# Find the error bars on the lifetime.
# We look for the lifetimes where the NLL is 0.5 worse than optimal;
# those correspond to 1-sigma error bars.
#optimal_NLL = NLL_vals[-1]
#def NegLogLikelihood_floatparameters(L):
#    xopt, NLL_vals = 
# lower_bound = 
# Use 
