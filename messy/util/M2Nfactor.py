###############################################################################
# This simple script provides the mass-to-number conversion factor for a
# lognormal distribution, given the median diameter D [nm], sigma [-] and
# the specific density rho (kg/m3), according to the equation
#
#     M2N = 6 / pi / exp(4.5 ln^2(sigma)) / D^3 / rho
#
# #############################################################################

import sys, math

if (len(sys.argv) < 4 or sys.argv[1] == '-h'):
    print 'Usage:'
    print '  M2Nfactor.py <diameter> <sigma> <rho>'
    quit()

args = map(float, sys.argv[1:])
D = float(args[0]) * 1.e-9
sigma = float(args[1])
rho = float(args[2])

fac = 6. / math.pi / math.exp(4.5 * math.log(sigma)**2) / D**3 / rho

print("%10.6e"% fac)
