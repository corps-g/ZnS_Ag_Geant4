import numpy as np
import scipy.integrate as integrate
import sys
import os

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

format = '{:.9f}'.format

def writer(lenPMMA, lenZSA, nUnit, particle, scin, che, history, gperEvent):
    """UI mac file writer
    
    Parameters:
        lenPMMA: double in cm
            Thick of PMMA layer
        lenZSA: double in um
            Thickness of ZSA layer
        nUnit: int
            # of PMMA+ZSA unit
        particle: string
            source particle
                'neutron': fission spectrum
                'gamma':   8 MeV
        history: int
            # of MC histories
    """
    
    lenPMMA = float(lenPMMA)
    lenZSA = float(lenZSA)
    nUnit = int(nUnit)
    assert nUnit * (lenPMMA + lenZSA * 1e-4) <= 20. 
    assert particle in ['neutron', 'gamma']
    history = int(history)
    
    if particle == 'neutron':
        dirName = './results/neutron/pmma%.2f_zsa%.2f_units%i' % \
                  (lenPMMA, lenZSA, nUnit)
    else:
        if scin and che:
            dirName = './results/gamma/pmma%.2f_zsa%.2f_units%i_scin_che_%ig' % \
                      (lenPMMA, lenZSA, nUnit, gperEvent)
        elif scin and not che:
            dirName = './results/gamma/pmma%.2f_zsa%.2f_units%i_scin_%ig' % \
                      (lenPMMA, lenZSA, nUnit, gperEvent)
        elif che and not scin:
            dirName = './results/gamma/pmma%.2f_zsa%.2f_units%i_che_%ig' % \
                      (lenPMMA, lenZSA, nUnit, gperEvent)
    os.system('rm -r ' + dirName)
    os.system('mkdir ' + dirName)
    with open(dirName + '/run.mac', 'w') as f:
        if particle == "gamma" and scin and not che:
            f.write('# Only scin\n')
            f.write('/process/optical/processActivation Cerenkov false\n\n')
        elif particle == 'gamma' and che and not scin:
            f.write('# Only cherenkov\n')
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
            f.write('/gps/particle neutron\n')
            f.write('/gps/direction 0 0 1\n')
        else:
            f.write('/gps/particle gamma\n')
            # isotropic gamma background shoot to ditector position
            s = """/gps/ang/type iso
/gps/ang/maxtheta 90 deg
"""
            f.write(s)
        
        f.write('\n# Position\n')
        f.write('/gps/pos/type Plane\n/gps/pos/shape Rectangle\n')
        srcZ = - 0.00001
        f.write('/gps/pos/centre 0 0 ' + format(srcZ) + ' cm\n')
        halfx, halfy = .1255, .4445 
        f.write('/gps/pos/halfx ' + str(halfx) + ' cm\n')
        f.write('/gps/pos/halfy ' + str(halfy) + ' cm\n')
        
        # Energy spectrum
        # Histogram can only be input by point, not file allowed
        # Point spectrum can do both.
        if particle == 'neutron':
            f.write('\n# neutron energy spectrum\n')
            f.write('/gps/ene/type User\n/gps/hist/type energy\n')
            erg, hist = neutron_spec()
            for i in range(len(erg)):
                f.write('/gps/hist/point ' + '{:.9f}'.format(erg[i]) + ' ' + 
                        format(hist[i]) + '\n')
        else:
            # gamma energy sampling
            s = '\n# gamma energy sampling\n'
            s += '/gps/ene/type Arb\n'
            s += '/gps/hist/type arb\n'
            erg, prob = gamma_spec()
            for i in range(len(erg)):
                e = erg[i]
                p = prob[i]
                s += '/gps/hist/point %.8f %.8f\n' % (e, p)
            s += '/gps/hist/inter Lin\n'
            f.write(s)
            f.write('/gps/number %i' % gperEvent)
        f.write('\n/run/beamOn %i'%history)
    return dirName

def submit(dirName, pmma, zsa, n):
    # Submit file
    os.system('cp ~/submit/submit.qsub ' + dirName + '/')
    with open(dirName + '/submit.qsub', 'a') as f:
#         f.write('\n#PBS -N r%.2f_%.2f_%i' % (pmma, zsa, n))
        f.write('\n%s/build/exampleB1 run.mac %.6f %.6f %i > out.txt' % \
                (os.path.abspath('.'), pmma, zsa, n))
    # Submit job
    os.system('cd ' + dirName + ' && qsub submit.qsub')

def submit_local_mathine(dirName, l1, l2, n):
    l1 = float(l1)
    l2 = float(l2)
    n = int(n)
    os.system('cd ' + dirName + ' && ' + os.path.abspath('.') + \
              '/build/exampleB1 run.mac ' + str(l1) + ' ' + \
              str(l2) + ' ' + str(n) + ' > out.txt')
                
if __name__ == '__main__':
    rootDir = os.listdir('.')
    if 'results' not in rootDir:
        os.system('mkdir results')
        os.system('mkdir ./results/gamma')
        
    history = 1e6
    length = 5. 
    lenpmma = np.logspace(-2, np.log(2)/np.log(10), 10)
    lenzsa = np.logspace(0, 2, 10)
    gperEvent = [10]
    for pmma in lenpmma:
        for zsa in lenzsa:
            for gamma in gperEvent:
                n = int(length / (pmma + zsa * 1e-4))
                dirg = writer(pmma, zsa, n, 'gamma', True, True, history, gamma)
                submit(dirg, pmma, zsa, n)
        
