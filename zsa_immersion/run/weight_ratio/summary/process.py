from nice_plot import *
import os
from scipy.signal import savgol_filter

def read_histogram(file, history):
    hist = np.loadtxt(file, delimiter = ',', skiprows = 8)[:, 0]
    """In Geant4 analysis file, there is a last bin computing number of events 
    with optical photons beyond the set edge. Simple add it into the last
    set bin to maintain same length as the set bin array.  
    """
    hist[-2] += hist[-1]
    hist = hist[0 : -1]
    eff = []
    s = 0
    for h in hist:
        s += h 
        eff.append(history - s)
    return np.array(eff) / history

def cal_sig(eff, history):
    sig = ((eff - eff ** 2) / (history - 1)) ** .5
    return sig

def main():
    history = 1e7
    op = [1]
    for i in range(0, 400):
        op.append((i+1) * 5. + .1)
    lenz = 5.
    a = np.linspace(.01, .3, 30)
#     b = np.arange(.4, .7, .1)
    # the running of weight ratio = 0.5 and 0.6 is not completed 
    b = np.arange(.4, .7, .1)
    wrs = np.append(a, b)
    gpe = 10
    lib = []
    for wr in wrs:
        # get data 
        tmp = [wr]
        nfile = '../results/neutron/lenz%.2f_wr%.2f/_h1_0.csv' % (lenz, wr)
        neff = read_histogram(nfile, history)
        gfile = '../results/gamma/lenz%.2f_wr%.2f_scin_che_10g/_h1_0.csv' \
                % (lenz, wr)
        geff = read_histogram(gfile, history)
        sn = neff / geff
        tmp.append(neff)
        tmp.append(geff)
        tmp.append(sn)
        # wr, neff, geff, sn
        lib.append(tmp)
    assert len(lib) == len(wrs)
    
    snlimit = 100.
    sn_value = []
    neff_sn = []
    for l in lib:
        neff, geff, sn = l[1 : ]
        for i in range(len(sn)):
            if sn[i] > snlimit:
                break
        if i != len(sn) - 1:
            sni1 = sn[i-1]
            sni = sn[i]
            neffi1 = neff[i-1]
            neffi = neff[i]
            # linear interplate between i and i-1 to get neutron-detection
            # efficiency at S/N ratio = 100
            tmp = neffi1 - (neffi1 - neffi) / (sni - sni1) * \
            (snlimit - sni1)
            neff_sn.append(tmp)
            # store the S/N ratio at i for inspectation
            sn_value.append(sn[i])
    neff_sn = np.array(neff_sn)
    print 'sn values'
    print sn_value
    loc = np.argmax(neff_sn)
    best_case = neff_sn[loc]
    sig = cal_sig(best_case, history)
    print 'best case, snlimit = %i, wr = %.4f, neff = %.4e, sig = %.4e' % \
    (snlimit, wrs[loc], best_case, sig)
        
    # parametric comparison
    plt.figure()
    sig = cal_sig(neff_sn, history)
    plt.figure()
    plt.errorbar(wrs * 100, neff_sn * 100, yerr = sig*100, 
                 fmt = 'none', ecolor='k')
    fit = np.polyval(np.polyfit(wrs  * 100, neff_sn  * 100, 14), 
                     wrs *100)
    plt.plot(wrs  * 100, fit, 'k-')
    plt.grid()
    plt.xlim(xmax = wrs[-1] * 100 + 2)
    plt.xlabel('ZnS(Ag) mass ratio (%)')
    plt.ylabel('Neutron-detection efficiency (%)')
    plt.savefig('immersion_wr.png')
    
if __name__ == '__main__':
    main()
    
