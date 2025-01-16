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
    return np.array(eff) / history, np.array(hist) / history

def cal_sig(eff, history):
    sig = ((eff - eff ** 2) / (history - 1)) ** .5
    return sig

def cal_divide_sig(var1, sig1, var2, sig2):
    # calculate the division sigma: divide = var1 / var2
    return ((1. / var2) ** 2 * sig1**2 + \
              (var1 / var2**2)**2 * sig2 ** 2) ** .5

def main():
    history = 1e7
    op = [1]
    for i in range(0, 400):
        op.append((i+1) * 5. + .1)
    wr = 0.12
    lenz = 5.
    gpe = 10
    lib = []
    # get data 
    nfile = '../results/neutron/lenz%.2f_wr%.2f/_h1_0.csv' % (lenz, wr)
    neff, nspec = read_histogram(nfile, history)
    nsig = cal_sig(neff, history)
    
    scfile = '../results/gamma/lenz%.2f_wr%.2f_scin_che_10g/_h1_0.csv' \
            % (lenz, wr)
    sceff, scspec = read_histogram(scfile, history)
    sc_sig = cal_sig(sceff, history)
    
    scinfile = '../results/gamma/lenz%.2f_wr%.2f_scin_10g/_h1_0.csv' \
            % (lenz, wr)
    scineff, scinspec = read_histogram(scinfile, history)
    
    chefile = '../results/gamma/lenz%.2f_wr%.2f_che_10g/_h1_0.csv' \
            % (lenz, wr)
    cheeff, chespec = read_histogram(chefile, history)
    
    sn = neff / sceff
    sn_sig = cal_divide_sig(neff, nsig, sceff, sc_sig)
    
    for i in range(len(sn)):
        if sn[i] > 100:
            break
    
    print '\nConsidering the scintillation and Cherenkov noises'
    print 'At i, LLD = %.i, sn =%.2f +- %.2f, neff = %.4e +- %.4e' %\
    (op[i], sn[i], sn_sig[i], neff[i], nsig[i])
    print 'At i-1, LLD = %.i, sn =%.2f +- %.2f, neff = %.4e +- %.4e' %\
    (op[i-1], sn[i-1], sn_sig[i-1], neff[i-1], nsig[i-1])
    
    # linear interpolation between sn[i] and sn[i-1] to compute efficiency at
    # S/N ratio = 100
    factor = abs(neff[i] - neff[i-1]) / abs(sn[i] - sn[i-1]) 
    diff = factor * abs(100. - sn[i-1])
    neff100 = neff[i-1] - diff
    print 'Linear interpolate neutron detection efficiency at S/N 100 = %.4e' \
    % neff100
    
    # pulse height plot
    plt.figure()
    plt.step(op, nspec, 'k-', label = 'neutron')
    plt.step(op, chespec, 'g-', label = 'gamma Che.')
    plt.step(op, scinspec, 'b-', label = 'gamma scin.')
    plt.step(op, scspec, 'r-', label = 'gamma Che.+scin.')
    plt.grid()
    plt.legend(loc = 0)
    plt.ylabel('Normalized number of events')
    plt.yscale('log')
    plt.xlabel('Number of OPs')
    plt.xlim(xmax=2e3)
    plt.savefig('immersion_pulse_height.png')
    
    # neff plot
    plt.figure()
    plt.errorbar(op, neff*100, yerr = nsig * 100, fmt='k-', ecolor='k')
    # sn 100 marker, use data at (i-1) since sn[i-1] = 99.45 is more close to 
    # 100
    plt.plot(op[i-1], neff[i-1]*100, 'ks', ms=8)
    plt.arrow(op[60], neff[60]*100, -40, 0, 
              head_width=.05, head_length=10, color='k')
    plt.grid()
    plt.xlabel('LLD (number of OPs)')
    plt.ylabel('Neutron-detection efficiency (%)')
    plt.xlim(xmax=500)
    # sn plot
    plt.twinx()
    plt.plot(op, sn, ':k')
    plt.fill_between(op, sn-sn_sig, sn+sn_sig, color='k', alpha=.4)
    
    # sn 100 marker
    plt.plot(op[i-1], sn[i-1], 'ks', ms=8)
    # connect the two sn 100 markers
    plt.plot((op[i-1], op[i-1]), (sn[i-1], 6.2), 'k--')
    
    plt.arrow(op[60], sn[60], 40, 0, 
              head_width=80, head_length=10, color='k')
    plt.xlim(xmax=500)
    plt.yscale('log')
    plt.ylim(ymax=1e4)
    plt.ylabel('S/N ratio')
    plt.savefig('immersion_5cm_neff.png')
        
    
if __name__ == '__main__':
    main()
    
