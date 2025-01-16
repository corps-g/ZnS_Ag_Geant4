from nice_plot import *

""" Simulate the Hornyak button's and the layer detector's responses to Cf252
    source.
"""

# histogram files
hb = '../half_inch_hb/results/cf252/_h1_0.csv'
half_inch_layer = '../layer/run/cf252/results/pmma0.02_zsa12.00_units59/_h1_0.csv'
one_inch_layer = '../layer/run/cf252/results/pmma0.02_zsa12.00_units119/_h1_0.csv'
four_cm_layer = '../layer/run/cf252/results/pmma0.02_zsa12.00_units188/_h1_0.csv'

def read_histogram(afile):
    # read the spectrum. Skip the last beyond max bin, start from the [0, 0.9)
    # bin
    spec = np.loadtxt(afile, delimiter = ',', skiprows = 8)[0 : -1, 0]
    return spec

history = 1e9
hb = read_histogram(hb) / history
hil = read_histogram(half_inch_layer) / history
oil = read_histogram(one_inch_layer) / history
fcl = read_histogram(four_cm_layer) / history

bin = np.arange(401) * 5 + 0.1

plt.figure()
plt.step(bin, hb, 'k-', label = 'Hornyak button')
plt.step(bin, hil, '-', color = 'yellowgreen', label = 'MLFD 12mm')
plt.step(bin, oil, 'r-', label = 'MLFD 25.4mm')
plt.step(bin, fcl, '-', color = 'deepskyblue', label = 'MLFD 40mm')
plt.yscale('log')
plt.legend(loc = 0)
plt.xlim(xmax = 2e3)
plt.xlabel('Number of optical photons')
plt.ylabel('Norm. number of events')
plt.grid()
plt.savefig('cf_spec.png')
