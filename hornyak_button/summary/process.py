from nice_plot import *
import os

def read_histogram(file, history):
    s = 0 
#     hist = []
#     with open(file, 'r') as f:
#         lines = f.readlines()[8 : -1]
#         for line in lines:
#             id = line.index(',')
#             hist.append(float(line[0 : id]))
    hist = np.loadtxt(file, delimiter = ',', skiprows = 8)[:, 0]
    """In Geant4 analysis file, there is a last bin computing number of events 
    with optical photons beyond the set edge. Simple add it into the last
    set bin to maintain same length as the set bin array.  
    """
    hist[-2] += hist[-1]
    hist = hist[0 : -1]
    eff = [] 
    for h in hist:
        s += h 
        eff.append(history - s)
    return hist / history, np.array(eff) / history

def cal_sig(eff, history):
    sig = ((eff - eff ** 2) / (history - 1)) ** .5
    return sig

def cal_divide_sig(var1, sig1, var2, sig2):
    # calculate the division sigma: divide = var1 / var2
    return ((1. / var2) ** 2 * sig1**2 + \
              (var1 / var2**2)**2 * sig2 ** 2) ** .5
    

def main():
    ghist = 5e7
    nhist = 1e7
    op = [1]
    for i in range(0, 400):
        op.append((i+1) * 5. + .1)
    gpe = 69
        
    neutronFile = '../results/neutron/abs13.25_sa1.0/neutron_h1_0.csv'
    gcheFile = '../results/gamma' + \
    '/abs13.25_sa1.0_che_%ig/gamma_h1_0.csv'%gpe
    gscinFile = '../results/gamma' + \
    '/abs13.25_sa1.0_scin_%ig/gamma_h1_0.csv'%10
    gscincheFile = '../results/gamma' + \
    '/abs13.25_sa1.0_scin_che_%ig/gamma_h1_0.csv'%gpe
    
    nspec, neff = read_histogram(neutronFile, nhist)
    chespec, gcheeff = read_histogram(gcheFile, ghist)
    scinspec, gscineff = read_histogram(gscinFile, ghist)
    scinchespec, gscincheeff = read_histogram(gscincheFile, ghist)
    
    # pulse height plot
    plt.step(op, nspec, 'k-', label = 'neutron')
    
    plt.step(op, chespec, 'g-', label = 'gamma Che.')
    
    plt.step(op, scinspec, 'b-', label = 'gamma scin.')
    
    plt.step(op, scinchespec, 'r-', label = 'gamma Che.+scin.')
    plt.yscale('log')
    plt.grid()
    plt.xlabel('Number of OPs')
    plt.ylabel('Normalized number of events')
    plt.xlim(xmax = 2e3)
    plt.legend(loc = 0)
    plt.savefig('pulse_height_hb.pdf')
    plt.savefig('pulse_height_hb.png')

    
    nsig = cal_sig(neff, nhist)
    scin_sig = cal_sig(gscineff, ghist)
    che_sig = cal_sig(gcheeff, ghist)
    scinche_sig = cal_sig(gscincheeff, ghist)
    
    snscin = neff / gscineff
    snscinche = neff / gscincheeff
    
    sn_sig = cal_divide_sig(neff, nsig, gscineff, scin_sig)
    sn_sig_scinche = cal_divide_sig(neff, nsig, gscincheeff, scinche_sig)
    
    # setting LLD
    snlimit = 100.
    for i in range(len(snscin)):
        if snscin[i] > snlimit:
            break
    iscin=i
    print '\nreject scintillation noise'
    print 'At i, sn=%.2f +- %.2f, LLD=%i, neff=%.4e +- %.4e'%\
    (snscin[i], sn_sig[i],  op[i], neff[i], nsig[i])
    print 'At i-1, sn=%.2f +- %.2f, LLD=%i, neff=%.4e +- %.4e'%\
    (snscin[i-1], sn_sig[i-1], op[i-1], neff[i-1], nsig[i-1])
#     # linear interpolation
#     factor = abs((neff[i-1] - neff[i]) / (snscin[i] - snscin[i-1]))
#     neff100 = neff[i-1] - factor * (100. - snscin[i-1])
#     op100 = op[i-1] + 5. * (100. - snscin[i-1]) * factor
#     print 'After linear interpolation, neff = %.4e, op = %.4f' % \
#     (neff100, op100)
    
    print 'reject scin.+Che. noise'
    for i in range(len(snscinche)):
        if snscinche[i] > snlimit:
            break
    
    iscinche=i
    print 'At i, sn=%.2f +- %.2f, LLD=%i, neff=%.4e +- %.4e'%\
    (snscinche[i], sn_sig_scinche[i],  op[i], neff[i], nsig[i])
    print 'At i-1, sn=%.2f +- %.2f, LLD=%i, neff=%.4e +- %.4e\n'%\
    (snscinche[i-1], sn_sig_scinche[i-1], op[i-1], neff[i-1], nsig[i-1])
    
    # neutron detection efficiency plot
    plt.figure()
    xmax=800
    plt.errorbar(op, neff * 100, 
                 yerr = nsig*100, 
                 fmt = 'k-', ecolor = 'k')
    
    # label the neutron detection efficiency at S/N 100
    # scintillation noise
    plt.plot(op[iscin], neff[iscin]*100, 'k*')
    # scintillation + Cherenkov noise
    plt.plot(op[iscinche], neff[iscinche]*100, 'ks', ms=8)
    
    # add point-left arrow labeling the neutron-detection efficiency curve
    plt.arrow(op[60], neff[60]*100, -60, 0, 
              head_width=.01, head_length=20, color='k')
    
    plt.xlim(xmax = xmax)
    plt.ylim(ymin = 3e-2)
    plt.grid()
    plt.xlabel('LLD (number of OPs)')
    plt.ylabel('Neutron-detection efficiency (%)')
    
    # S/N ratio plot
    plt.twinx()
    # S/N ratio counts scintillation noise
    plt.plot(op, snscin, '--k')
    plt.fill_between(op, snscin-sn_sig, snscin+sn_sig, color='k',
                     facecolor='k', alpha=.4)
    
    # label S/N = 100 point
    plt.plot(op[iscin], snscin[iscin], 'k*' )
    # line connecting S/N 100 and corresponding neutron efficiency
    plt.plot((op[iscin], op[iscin]), (snscin[iscin], .8), 'k--')
    
    plt.annotate('S/N ratio\nscin. only', (op[60], snscin[60]),
                 (op[60]-120, snscin[60]+1.2e3),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1,
                                 headwidth=3),
                 ha='center', va='bottom')
    
    # S/N ratio considering scintillation + Cherenkov noise
    plt.plot(op, snscinche, ':k', ms=8)
    plt.fill_between(op, snscinche-sn_sig_scinche, snscinche+sn_sig_scinche,
                     facecolor='k', alpha=.4)
    plt.annotate('S/N ratio\nChe.+scin.', (op[70], snscinche[70]),
                 (op[70]+80, .5),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1,
                                 headwidth=3),
                 ha='left', va='bottom')
    
    # label S/N = 100 point and connect it with neutron-detection efficiency
    plt.plot(op[iscinche], snscinche[iscinche], 'sk', ms=8)
    plt.plot((op[iscinche], op[iscinche]), (snscinche[iscinche], 2e-2), 
             'k--')
    
    plt.ylabel('S/N ratio')
    plt.xlim(xmax = xmax)
    plt.ylim((1e-2, 1e4))
    plt.yscale('log')
    plt.savefig('neff_hb.pdf')
    plt.savefig('neff_hb.png')
    
   
    
if __name__ == '__main__':
    main()
    
