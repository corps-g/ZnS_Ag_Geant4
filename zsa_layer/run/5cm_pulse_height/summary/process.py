from nice_plot import *
import os

def read_histogram(file, history):
    s = 0 
    hist = []
    with open(file, 'r') as f:
        lines = f.readlines()[8 : -1]
        for line in lines:
            id = line.index(',')
            hist.append(float(line[0 : id]))
    eff = [] 
    for h in hist:
        s += h 
        eff.append(history - s)
    return (np.array(eff) / history, np.array(hist) / history)

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
    nhist = 1e7
    ghist = 5e7
    op = [1]
    for i in range(0, 400):
        op.append((i+1) * 5. + .1)
    length = 5.
    # best for sn 100 and 500
    pmma = 0.018016
    zsa = 12.915497
    snlimit = 100.

    units = int(length / (pmma + zsa * 1e-4))
    nfile = '../results/neutron/pmma%.2f_zsa%.2f_units%i/_h1_0.csv' % \
    (pmma, zsa, units)
    gfileSC = '../results/gamma/pmma%.2f_zsa%.2f_units%i_scin_che_10g/_h1_0.csv' % \
    (pmma, zsa, units)
    gfileS = '../results/gamma/pmma%.2f_zsa%.2f_units%i_scin_10g/_h1_0.csv' % \
    (pmma, zsa, units)
    gfileC = '../results/gamma/pmma%.2f_zsa%.2f_units%i_che_10g/_h1_0.csv' % \
    (pmma, zsa, units)
    
    # pulse height
    nph = read_histogram(nfile, nhist)[1]
    scph = read_histogram(gfileSC, ghist)[1]
    sph = read_histogram(gfileS, ghist)[1]
    cph = read_histogram(gfileC, ghist)[1]
    
    # pulse height plot
    plt.figure()
    plt.step(op, nph, 'k-', label = 'neutron')
    plt.step(op, cph, 'g-', label = 'gamma Che.')
    plt.step(op, sph, 'b-', label = 'gamma scin.')
    plt.step(op, scph, 'r-', label = 'gamma Che.+scin.')
    plt.grid()
    plt.xlabel('Number of OPs')
    plt.ylabel('Normalized number of events')
    plt.yscale('log')
    plt.xlim(xmax=2e3)
    plt.legend(loc=0)
    plt.savefig('layer_5cm_ph.pdf')

    
    # neutron detection efficiency as a function of LLD
    neff = read_histogram(nfile, nhist)[0]
    geff = read_histogram(gfileSC, ghist)[0]
    sn = neff / geff
    for i in range(len(sn)):
        if sn[i] > 100.:
            break
    print 'S/N = %.2f, sni1=%.2f, LLD = %.2f, neff = %.4e, neffi1=%.4e' \
    % (sn[i],sn[i-1], op[i], neff[i], neff[i-1])
    sig = cal_sig(neff, nhist)
    gsig = cal_sig(geff, ghist)
    sn_sig = cal_divide_sig(neff, sig, geff, gsig)
    relerr = sig / neff
    print max(relerr[0 : i+1])
    plt.figure()
    # neff plot
    plt.errorbar(op, neff*100, yerr=sig*100, fmt='k-', ecolor = 'k')
    # marker sn 100
    plt.plot(op[i], neff[i]*100, 'sk', ms=8)
    plt.arrow(op[60], neff[60]*100, -30, 0, 
              head_width=.1, head_length=10, color='k')
    plt.grid()
    plt.xlabel('LLD (number of OPs)')
    plt.ylabel('Neutron-detection efficiency (%)')
    plt.xlim(xmax=500)
    # sn plot
    plt.twinx()
    plt.plot(op, sn, 'k:')
    plt.fill_between(op, sn-sn_sig, sn+sn_sig, facecolor='k', alpha=.4)
    plt.plot(op[i], sn[i], 'sk', ms=8)
    plt.plot((op[i], op[i]), (sn[i], 7), 'k--')
    plt.arrow(op[60], sn[60], 30, 0, 
              head_width=60, head_length=10, color='k')
    plt.xlim(xmax=500)
    plt.yscale('log')
    plt.ylim(ymax=1e4)
    plt.ylabel('S/N ratio')
    plt.savefig('layer_5cm_neff.pdf')
    
    
if __name__ == '__main__':
    main()
    
