from nice_plot import *
import os

def read_histogram(file, history = int(1e6)):
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
    return np.array(eff) / history

def main():
    history = 1e6
    op = [1]
    for i in range(0, 100):
        op.append((i+1) * 20. + .1)
    lenPMMA = np.logspace(-2, np.log(2)/np.log(10), 10)  # 0.01 ~ 2 cm
    lenZSA = np.logspace(0, 2, 10)
    length = 5. 
    gpe = 10
    lib = []
    for pmma in lenPMMA:
        for zsa in lenZSA:
            # eliminate fail cases
            units = int(length / (pmma + zsa * 1e-4))
            motherDir = '../results/neutron'
            ndir = '/pmma%.2f_zsa%.2f_units%i' % \
            (pmma, zsa, units)
            files = os.listdir(motherDir + ndir)
            if not '_h1_0.csv' in files:
                print 'neutron fail, pmma=%.2f, zsa=%.2f'%(pmma, zsa)
                continue
            # get data 
            nfile = motherDir + ndir + '/_h1_0.csv'
            neff = read_histogram(nfile, history)
            gdir = '../results/gamma/pmma%.2f_zsa%.2f_units%i_scin_che_10g'\
                    %(pmma, zsa, units)
            if not '_h1_0.csv' in os.listdir(gdir):
                print 'gamma fail, pmma=%.2f, zsa=%.2f'%(pmma, zsa)
                continue
            gfile = gdir + '/_h1_0.csv'
            geff = read_histogram(gfile, history)
            sn = neff / geff
            tmp = [pmma, zsa]
            tmp.append(neff)
            tmp.append(geff)
            tmp.append(sn)
            # pmma, zsa, neff, geff, sn
            lib.append(tmp)
    print 'length of lib = %i'%len(lib)
    neff_sn = [] # parametric study
    snlimit = [100]
    data = {}
    for limit in snlimit:
        neff_sn = []
        for l in lib:
            neff = l[2]
            geff = l[3]
            sn = l[4]
            iend = 0
            for i in range(len(sn)):
                if np.isinf(sn[i]):
                    iend = i 
            if iend == 0:
                iend = len(neff)
#             for i in range(len(sn)):
#                 if sn[i] > limit:
#                     break
#             if i == len(sn) - 1:
#                 print 'sn not meet, pmma=%.2f, zsa=%.2f'%(l[0], l[1])
#             else:
#                 sni1 = sn[i-1]
#                 sni = sn[i]
#                 neffi1 = neff[i-1]
#                 neffi = neff[i]
#                 tmp = neffi1 - (neffi1 - neffi) / (sni - sni1) * (limit - sni1)
#                 neff_sn.append([l[0], l[1], tmp])
            obj = neff[np.ix_(sn[0 : iend] > limit)]
            if len(obj) != 0:
                maxloc = np.argmax(obj)
                maxeff = obj[maxloc]
                neff_sn.append([l[0], l[1], maxeff])
            else:
                print 'pmma=%.2f, zsa=%.2f'%(l[0], l[1])
        data[limit] = np.array(neff_sn)
    
    for sn in data.keys():
        neff_sn = data[sn] 
        print 'length of neff_sn = %i'%len(neff_sn)
        loc = np.argmax(neff_sn[:, 2])
        best_case = neff_sn[loc]
        print 'best case, snlimit = %i, pmma=%.4f, zsa=%.4f, maxeff = %.4e' % \
        (sn, best_case[0], best_case[1], best_case[2])
        
    # print out the cases with threshhold efficiency
    neff_sn = data[100]
    neff_2 = neff_sn[np.ix_(neff_sn[:, -1]>.02)]
    print 'number of cases satisfying > 2%% = %i'%len(neff_2)
    for item in neff_2:
        print 'pmma=%.2f, zsa=%i, neff = %.4e'%(item[0]*10, item[1], item[2])
    exit()
    
    
    # parametric comparison
    for sn in data.keys():
        neff_sn = data[sn]
        # a 7 x 9 matrix
        ans = np.zeros((7, 9))
        for p in range(0, 7):
            pmma = np.logspace(-2, np.log(2)/np.log(10), 10)[p]
            for z in range(1, 10):
                zsa = np.logspace(0, 2, 10)[z]
                for item in neff_sn:
                    if abs(item[0] - pmma) < 1e-7 and abs(item[1]-zsa)<1e-7:
                        
                        ans[p, z-1] = item[-1]
        print ans
        ans = np.transpose(ans)
        plt.figure()
        plt.imshow(ans, origin='lower', cmap=plt.cm.Greys,interpolation='nearest')
        plt.colorbar(orientation='vertical')
        pmmaticks = []
        for thick in lenPMMA[0 : 7]*10:
            pmmaticks.append('%.1f'%thick)
        zsaticks = []
        for thick in lenZSA[1 : ]:
            zsaticks.append('%.1f'%thick)
        plt.xticks(range(7), pmmaticks)
        plt.yticks(range(9), zsaticks)
        plt.xlabel('Thickness of PMMA (mm)')
        plt.ylabel('Thickness of ZnS(Ag) ($\mu$m)')
        plt.savefig('thickness_opt_layer.pdf')
                
    # parametric comparison
    for sn in data.keys():
        neff_sn = data[sn]
        tpmma = neff_sn[:, 0]
        tzsa = neff_sn[:, 1]
        maxeff = neff_sn[:, 2]
        plt.figure()
        plt.tripcolor(tpmma, tzsa, maxeff, cmap = plt.cm.Greys)
        plt.xscale('log')
        plt.yscale('log')
        plt.colorbar()
        plt.xlabel('Thickness of PMMA (cm)')
        plt.ylabel('Thickness of ZnS(Ag) ($\mu$m)')
        plt.xlim(xmax = .18)
        plt.ylim(ymin = 1.66810054)
        plt.savefig('sn%i_layer.pdf'%sn)
        
    # pmma = 0.06cm case
    print 'thickness of pmma = %.2f cm' % lenPMMA[3]
    
    sn100 = data[100]
    pmma6 = sn100[np.ix_(sn100[:, 0] == lenPMMA[3])]
    loc = np.argmax(pmma6[:, -1])
    print 'max eff. at pmma %.2f cm, zsa %.2f um, eff = %.4e' % \
    (pmma6[loc, 0], pmma6[loc, 1], pmma6[loc, 2])
    plt.figure()
    plt.plot(pmma6[:, 1], pmma6[:, -1] * 100, '*-k')
    plt.xlabel('Thickness of ZnS(Ag) [um]')
    plt.ylabel('Detection efficiency (sn > 100) [%]')
    plt.grid()
    plt.savefig('pmma6.pdf')
    
    pmma, zsa = pmma6[loc, 0], pmma6[loc, 1]
    for entry in lib:
        if entry[0] == pmma and entry[1] == zsa:
            data = entry
    pmma, zsa, neff, geff, sn = data[:]
    for i in range(len(sn)):
        value = sn[i]
        if np.isinf(value):
            iend = i
            break
    fig, ax1 = plt.subplots()
    plt.plot(op[0 : iend], neff[0 : iend]*100, 'k-', label = 'neutron eff.')
    plt.plot(op[0 : iend], geff[0 : iend]*100, 'k:', label = 'gamma eff.')
    plt.ylabel('Detection efficiency [%]')
    plt.xlabel('Optical photon discriminator')
    plt.grid()
    ax2 = ax1.twinx()
    plt.plot(op[0 : iend], sn[0 : iend], 'k--')
    plt.ylabel('sn ratio')
    plt.yscale('log')
    plt.arrow(200, 100, 30, 0, width=5, head_width=15, fc = 'k')
    plt.savefig('neff_sn_pmma6.pdf')
    
    
if __name__ == '__main__':
    main()
    
