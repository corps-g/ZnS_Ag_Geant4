from nice_plot import *
import os

def read_histogram(file, history):
    hist = np.loadtxt(file, delimiter = ',', skiprows = 8)[:, 0]
    """In Geant4 analysis file, there is a last bin computing number of events 
    with optical photons beyond the set edge. Simple add it into the last
    set bin to maintain same length as the set bin array.  
    """
    hist[-2] += hist[-1]
    hist = hist[0 : -1]
    
    s = 0 
    eff = [] 
    for h in hist:
        s += h 
        eff.append(history - s)
    return (np.array(eff) / history, hist / history)

def cal_sig(eff, history):
    sig = ((eff - eff ** 2) / (history - 1)) ** .5
    return sig

def cal_divide_sig(var1, sig1, var2, sig2):
    # calculate the division sigma: divide = var1 / var2
    return ((1. / var2) ** 2 * sig1**2 + \
              (var1 / var2**2)**2 * sig2 ** 2) ** .5
              
def find_iend(sn):
    for i in range(len(sn)):
        if np.isinf(sn[i]):
            return i
    return len(sn)

def main():
    history = 1e7
    op = [1]
    for i in range(0, 400):
        op.append((i+1) * 5. + .1)
    lengths = [1, 2, 2.54] + range(3, 41)
    # best for sn 100 and 500
    pmma = 0.018016
    zsa = 12.915497
    snlimit = 100.

    lib = []
    for length in lengths:
        tmp = [length]
        units = int(length / (pmma + zsa * 1e-4))
        nfile = '../results/neutron/pmma%.2f_zsa%.2f_units%i/_h1_0.csv' % \
        (pmma, zsa, units)
        neff = read_histogram(nfile, history)[0]
        
        scinche = '../results/gamma/pmma%.2f_zsa%.2f_units%i_scin_che_10g/_h1_0.csv' % \
        (pmma, zsa, units)
        geff = read_histogram(scinche, history)[0]
        
        sn = neff / geff
        
        tmp.append(neff)
        tmp.append(geff)
        tmp.append(sn)
        
        lib.append(tmp)
    
    neff_sn = [] # parametric study
    sn_value = []
    for data in lib:
        neff, geff, sn = data[1 : ]
        
        for i in range(len(sn)):
            if sn[i] > snlimit:
                break
        
        sni1 = sn[i-1]
        sni = sn[i]
        neffi1 = neff[i-1]
        neffi = neff[i]
        
        # linear interpolation between the two bins
        tmp = neffi1 - (neffi1 - neffi) / (sni - sni1) * (snlimit - sni1)
        neff_sn.append(tmp)
        sn_value.append(sn[i])
    neff_sn = np.array(neff_sn)
    sig = cal_sig(neff_sn, history)
    print 'neff at length'
    print neff_sn
    print 'sn values at lengths'
    print sn_value
    
    
    plt.errorbar(lengths, neff_sn * 100, yerr = sig * 100, \
                 fmt = '-k', ecolor = 'k')
    plt.grid()
    plt.xlabel('Detector length (cm)')
    plt.ylabel('Neutron-detection efficiency (%)')
    plt.savefig('length_test_layer.pdf')
    
if __name__ == '__main__':
    main()
    
