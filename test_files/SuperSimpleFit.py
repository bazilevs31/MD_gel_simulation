from pylab import *
ion()
import fit
from numpy import random, exp
random.seed(0)

# Create some data to fit
x = arange(0, 10, .1)
# A gaussian of height 10, width 2, centered at zero. With noise.
# y = 10*exp(-x**2/8) + (random.rand(100) - 0.5)
y = (1. - exp(-20.*x**3)) + (random.rand(len(x)) - 0.5)*0.1

# # No need to provide first guess at parameters for fit.gaus
# (xf, yf), params, err, chi = fit.fit(fit.gaus, x,y)

# print "N:    %.2f +/- %.3f" % (params[0], err[0])
# print "N:    %.2f +/- %.3f" % (params[1], err[1])
# print "N:    %.2f +/- %.3f" % (params[2], err[2])


def example_function(params, x):
    kappa, n = params
    return (1. - exp(-kappa*x**n))

# It will still try to guess parameters, but they are dumb!
(xf,yf), p, e, chi = fit.fit(example_function, x,y)
plot(x,y, 'bo', label='Data')
# plot(xf,yf, 'r-', label='Fit')
# errorbar(xf,yf,yerr=e,'r-', label='Fit')
errorbar(xf,yf,yerr=chi)

legend()

# results = fit.fit(example_function, x, y, default_pars = [1, 12, 10, 1, 1, 1])
# plot(results[0][0], results[0][1], 'r--')

# Fit a sub-range:

# clf()
# results = fit.fit(fit.gaus, x, y, data_range=[0, 23])
# plot(results[0][0], results[0][1], 'r-.')
