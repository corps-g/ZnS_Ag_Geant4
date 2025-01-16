import numpy as np
import sys
import os
import scipy.integrate as integrate

def gamma_spec():
    # 0.086 to 8 MeV 
    erg = np.logspace(-1.0655015487564323, 0.90308998699194343, 1000)
    prob = []
    for e in erg:
        if e < .3:
            prob.append(38.13 * (e - .085) * np.e ** (1.648 * e))
        elif e > .3 and e < 1.:
            prob.append(26.8 * np.e ** (-2.3 * e))
        elif e > 1. and e < 8.:
            prob.append(8. * np.e ** (-1.1 * e))
    return erg, prob 

def Spectrum(erg):
    return np.exp(-erg / 0.988) * np.sinh((2.249 * erg) ** 0.5)

def neutron_spec():
    # Construct neutron energy spectrum
    erg = np.linspace(0., 20, 1000)
    hist = [0.]
    for i in range(1, len(erg)):
        ans = integrate.quad(Spectrum, erg[i - 1], erg[i])
        if ans[1] > 1e-12:
            print "Bad result, i = ", i
            break
        hist.append(ans[0])
    return erg, hist

def writer(history, absLen, sa, particle, scin, cherenkov, gperEvent):
    history = int(history)
    assert particle in ['neutron', 'gamma']
    
    
    lenz = 0.5 * 2.54
    halfx, halfy = 5. / 8. * 2.54 / 2., 7. / 64. * 2.54 / 2.
    
    if particle == 'neutron':
        dirName = './results/neutron/abs%.2f_sa%.1f' % \
        (absLen, sa)
    elif particle == 'gamma':
        if cherenkov and not scin:
            dirName = './results/gamma/abs%.2f_sa%.1f_che_%ig' % \
            (absLen, sa, gperEvent)
        elif scin and not cherenkov:
            dirName = './results/gamma/abs%.2f_sa%.1f_scin_%ig' % \
            (absLen, sa, gperEvent)
        elif scin and cherenkov:
            dirName = './results/gamma/abs%.2f_sa%.1f_scin_che_%ig' % \
            (absLen, sa, gperEvent)
    os.system('rm -r ' + dirName)
    os.system('mkdir ' + dirName)
    with open(dirName + '/run.mac', 'w') as f:
        if particle == "gamma":
            if scin and not cherenkov:
                f.write('# Disable cherenkov\n')
                f.write('/process/optical/processActivation Cerenkov false\n\n')
            elif cherenkov and not scin:
                f.write('# Disable scintillation\n')
                f.write('/process/optical/processActivation Scintillation false\n\n')
                
        # Geometry definition
        f.write('# Initialize kernel\n/run/initialize\n')
        
        f.write('/control/verbose 0\n/run/verbose 0\n/event/verbose 0\n' +
                '/tracking/verbose 0\n')
        
        # GPS definition
        f.write('#----------------------------------------------------------\n')
        f.write('# GPS definition\n')
        f.write('#----------------------------------------------------------\n')
        if particle == 'neutron':
            f.write('\n# neutron energy spectrum\n')
            f.write('/gps/particle neutron\n')
            f.write('/gps/ene/type User\n/gps/hist/type energy\n')
            erg, hist = neutron_spec()
            for i in range(len(erg)):
                f.write('/gps/hist/point %.9f %.9f\n' % (erg[i], hist[i]))
            f.write('/gps/direction 0 0 1\n')
            f.write('\n# Position\n')
            f.write('/gps/pos/type Plane\n/gps/pos/shape Rectangle\n')
            srcZ = -lenz * 0.5 - 1e-6
            f.write('/gps/pos/centre 0 0 %.9f cm\n' % srcZ)
            f.write('/gps/pos/halfx ' + str(halfx) + ' cm\n')
            f.write('/gps/pos/halfy ' + str(halfy) + ' cm\n')
            f.write('/analysis/setFileName neutron\n')
        else:
            f.write('/gps/particle gamma\n')
            # energy sampling
            s = '\n# energy sampling\n'
            s += '/gps/ene/type Arb\n'
            s += '/gps/hist/type arb\n'
            erg, prob = gamma_spec()
            for i in range(len(erg)):
                e = erg[i]
                p = prob[i]
                s += '/gps/hist/point %.8f %.8f\n' % (e, p)
            s += '/gps/hist/inter Lin\n'
            s += """# gamma position sampling
/gps/pos/type Plane
/gps/pos/shape Rectangle
/gps/pos/centre 0 0 SRCZ cm
/gps/pos/halfx HALFX cm
/gps/pos/halfy HALFY cm
# half isotropic gamma
/gps/ang/type iso
/gps/ang/maxtheta 90 deg
# multiple gammas per event
/gps/number GPEREVENT
/analysis/setFileName gamma
"""
            srcz = -.5 * lenz - 1e-6
            s = s.replace('SRCZ', '%.9f' % srcz)
            
            if cherenkov:
                s += '/gps/pos/confine srcPV'
                diameter = 3. / 4. * 2.54;
                s = s.replace('HALFX', '%.9f' % diameter)
                s = s.replace('HALFY', '%.9f' % diameter)
            else:
                s = s.replace('HALFX', '%.9f' % halfx)
                s = s.replace('HALFY', '%.9f' % halfy)
            s = s.replace('GPEREVENT', '%i'%gperEvent)
            f.write(s)
            
        f.write('\n/run/beamOn ' + str(history))
        
    return dirName 

def submit(dir, absLen, sa):
    os.system('cp ~/submit/submit.qsub ' + dir + '/submit.qsub')
    with open(dir + '/submit.qsub', 'a') as f:
        f.write(os.path.abspath('.') + '/build/exampleB1 run.mac ' + 
                    str(absLen) + ' ' + str(sa) + ' > out.txt')
    os.system('cd ' + dir + " && qsub submit.qsub")
    return 
    

if __name__ == '__main__':
    rootDir = os.listdir('.')
    if 'results' not in rootDir:
        os.system('mkdir results')
        os.system('mkdir results/neutron')
        os.system('mkdir results/gamma')
         
    ghist = 5e7
    nhist = 1e7
    frac = 0.01
    absLen = -61. / np.log(frac)
    sa = 1.
    
    gperEvent = 69
    ndir = writer(nhist, absLen, sa, 'neutron', True, True, gperEvent)
    submit(ndir, absLen, sa)
    
    scin_che = writer(ghist, absLen, sa, 'gamma', True, True, gperEvent)
    submit(scin_che, absLen, sa)
    che = writer(ghist, absLen, sa, 'gamma', False, True, gperEvent)
    submit(che, absLen, sa)
    gperEvent = 10
    scin = writer(ghist, absLen, sa, 'gamma', True, False, gperEvent)
    submit(scin, absLen, sa)
    
