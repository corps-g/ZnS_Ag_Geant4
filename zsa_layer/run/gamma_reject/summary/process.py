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
    # layer plot
    history = 1e6
    op = [1]
    for i in range(0, 100):
        op.append((i+1) * 20. + .1)
    length = 5.
    
    # layered detector
    print '\nLayered detector'
    pmma = 0.018016
    zsa = 12.915497
    snlimit = 100.
    gpes = [1]+range(10, 101, 10)
    n = int(length / (pmma + zsa * 1e-4))
    nfile = '../results/layer/neutron/pmma%.2f_zsa%.2f_units%i/_h1_0.csv' % \
    (pmma, zsa, n)
    neff = read_histogram(nfile, history)[0]
    lib = []
    for gpe in gpes: 
        gfile = '../results/layer/gamma/pmma%.2f_zsa%.2f_units%i_scin_che_%ig/_h1_0.csv' % \
                      (pmma, zsa, n, gpe)
        geff = read_histogram(gfile, history*10)[0]
        sn = neff / geff
        for i in range(len(sn)):
            if sn[i] > 100:
                break
#         sni = sn[i]
#         sni1 = sn[i-1]
#         neffi = neff[i]
#         neffi1 = neff[i-1]
#         factor = (neffi1 - neffi) / (sni - sni1)
#         ineff = neff[i-1] - factor * (100 - sni1)
        ineff = neff[i]
        sig_neff = cal_sig(ineff, history)
        igeff = geff[i]
        sig_geff = cal_sig(igeff, history * 10)
        isn = sn[i]
        sig_sn = cal_divide_sig(ineff, sig_neff, igeff, sig_geff)
        lib.append([ineff, sig_neff, isn, sig_sn])
        print 'gpe = %i, neff = %.4e, sig_neff = %.4e, sn = %.2f, sig_sn = %.2f' % \
        (gpe, ineff, sig_neff, isn, sig_sn)
        
    lib = np.array(lib)
    plt.errorbar(gpes, lib[:, 0]*100, yerr = lib[:, 1]*100, fmt = 'k-s', 
                 ecolor = 'k', label='layered', ms=8)
    plt.xlim(xmax=102)
    
    # immersion plot
    print '\nImmersion detector'
    history = 1e6
    op = [1]
    for i in range(0, 100):
        op.append((i+1) * 20. + .1)
    length = 5.
    
    gpes = [1] + range(10, 101, 10)
    wr = 0.12
    nfile = '../results/immersion/neutron/lenz%.2f_wr%.2f/_h1_0.csv' % \
    (length, wr)
    neff = read_histogram(nfile, history)[0]
    lib = []
    for gpe in gpes: 
        gfile = '../results/immersion/gamma/lenz%.2f_wr%.2f_scin_che_%ig/_h1_0.csv' % \
                      (length, wr, gpe)
        geff = read_histogram(gfile, history*10)[0]
        sn = neff / geff
        
        for i in range(len(sn)):
            if sn[i] > 100:
                break
        
        ineff = neff[i]
        sig_neff = cal_sig(ineff, history)
        igeff = geff[i]
        sig_geff = cal_sig(igeff, history * 10)
        isn = sn[i]
        sig_sn = cal_divide_sig(ineff, sig_neff, igeff, sig_geff)
        lib.append([ineff, sig_neff, isn, sig_sn])
        print 'gpe = %i, neff = %.4e, sig_neff = %.4e, sn = %.2f, sig_sn = %.2f' % \
        (gpe, ineff, sig_neff, isn, sig_sn)
    lib = np.array(lib)
    plt.errorbar(gpes, lib[:, 0]*100, yerr = lib[:, 1]*100, fmt = 'k-*', 
                 ecolor = 'k', label='homogenized', ms=8)
    plt.xlim(xmax = 102)
    plt.xlabel('Number of gamma rays per event')
    plt.ylabel('Neutron-detection efficiency (%)')
    plt.grid()
    plt.legend(loc = 0)
    plt.yscale('log')
    plt.savefig('gamma_rejection_compare.pdf')
    
    
if __name__ == '__main__':
    main()
    
