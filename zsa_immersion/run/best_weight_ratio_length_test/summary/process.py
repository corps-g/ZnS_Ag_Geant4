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
    return np.array(eff) / history

def cal_sig(eff, history):
    sig = ((eff - eff ** 2) / (history - 1)) ** .5
    return sig

def main():
    history = 1e7
    op = [1]
    for i in range(0, 400):
        op.append((i+1) * 5. + .1)
    wr = 0.12
    a = [1., 2., 2.54]
    b = np.arange(3., 41.)
    lenzs = np.append(a, b)
    gpe = 10
    lib = []
    for lenz in lenzs:
        # get data 
        tmp = [lenz]
        nfile = '../results/neutron/lenz%.2f_wr%.2f/_h1_0.csv' % (lenz, wr)
        neff = read_histogram(nfile, history)
        gfile = '../results/gamma/lenz%.2f_wr%.2f_scin_che_10g/_h1_0.csv' \
                % (lenz, wr)
        geff = read_histogram(gfile, history)
        sn = neff / geff
        tmp.append(neff)
        tmp.append(geff)
        tmp.append(sn)
        
        # lenz, neff, geff, sn
        lib.append(tmp)
        
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
            tmp = neffi1 - (neffi1 - neffi) / (sni - sni1) * \
            (snlimit - sni1)
            neff_sn.append(tmp)
            sn_value.append(sn[i])
    neff_sn = np.array(neff_sn)
    
    print 'neff applying LLD'
    print neff_sn
            
    sig = cal_sig(neff_sn, history)
    plt.figure()
    plt.errorbar(lenzs, neff_sn * 100, yerr = sig*100, 
                 fmt = 'none', ecolor='k')
    fit = np.polyval(np.polyfit(lenzs, neff_sn * 100, 6), lenzs)
    plt.plot(lenzs, fit, 'k-')
    plt.grid()
    plt.xlabel('Detector length (cm)')
    plt.ylabel('Neutron-detection efficiency (%)')
    plt.savefig('immersion_length.png')
        
    
if __name__ == '__main__':
    main()
    
