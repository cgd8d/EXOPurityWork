'''
Likelihood method
Correct the data; then fit it against a model function, ie:

P(E, drift_time) = (2*pi*sigma^2)^(-1/2) *
                   exp(-[E exp( drift_time /L) - E_0]^2/[2*sigma^2])

data is read in as an array of (charge_energy, drift_time) tuples.

To run this script, you should already have a csv file containing events.
Do: python Method2.py <csv filename>
'''

import math
import scipy.optimize
import numpy
import csv
import sys
import matplotlib.pyplot

def Prob(L, E, drift_time, sigma, E_0):
    purity_corr = math.exp(drift_time/L)
    ret_val = math.exp(-pow(E*purity_corr - E_0, 2)/(2*sigma*sigma))
    ret_val /= math.sqrt(2*math.pi)*sigma
    return ret_val

def NegLogLikelihood(L, sigma, E_0, EData, DriftTimeData):
    if len(EData) != len(DriftTimeData):
        print "Energy and drift time data are mismatched."
        sys.exit(1)
    ret_val = 0.
    for i in range(len(EData)):
        ret_val -= math.log(Prob(L, EData[i], DriftTimeData[i], sigma, E_0))
    return ret_val

def NLL_L_sigma_E0(x, EData, DriftTimeData): # For floating all three parameters.
    return NegLogLikelihood(x[0], x[1], x[2], EData, DriftTimeData)

def NLL_sigma_E0(x, L, EData, DriftTimeData): # For floating only sigma and E_0.
    return NegLogLikelihood(L, x[0], x[1], EData, DriftTimeData)

# Import data from csv file.
ChargeEnergyData = []
DriftTimeData = []
with open(sys.argv[1], 'rb') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if float(row[1]) < 2200: continue # Skip for now, need to figure out why it crept in.
        ChargeEnergyData.append(float(row[1]))
        DriftTimeData.append(float(row[0]))

xopt = scipy.optimize.fmin(NLL_L_sigma_E0,
                                     (1000., 100., 2600.),
                                     args = (ChargeEnergyData, DriftTimeData))
print "Result of optimization:"
print "Lifetime =", xopt[0], "us."
print "Charge sigma =", xopt[1], "keV."
print "Charge peak =", xopt[2], "keV."

# Graphics to show the result.
def PeakPositionFunc(drift_time):
    return xopt[2] * math.exp(-drift_time/xopt[0])
funcxpts = numpy.arange(0, 120, 0.1)
funcypts = [PeakPositionFunc(x) for x in funcxpts]
matplotlib.pyplot.plot(DriftTimeData, ChargeEnergyData, 'b.')
matplotlib.pyplot.plot(funcxpts, funcypts, 'r-')
matplotlib.pyplot.show()

# Plot NLL vs (fixed) lifetime.
NLL_xpts = numpy.arange(5000, 15000, 1000)
NLL_ypts = [scipy.optimize.fmin(NLL_sigma_E0,
                                (100., 2600.),
                                args = (x, ChargeEnergyData, DriftTimeData),
                                retall = True,
                                full_output = True)[1] for x in NLL_xpts]
matplotlib.pyplot.plot(NLL_xpts, NLL_ypts, 'ro')
matplotlib.pyplot.show()
