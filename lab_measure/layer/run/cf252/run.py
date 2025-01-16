import numpy as np
import os

def writer(lenPMMA, lenZSA, nUnit, history):
    """UI mac file writer
    
    Parameters:
        lenPMMA: double in cm
            Thickness of PMMA layer
        lenZSA: double in um
            Thickness of ZSA layer
        nUnit: int
            # of PMMA+ZSA unit
        history: int
            # of MC histories
    """
    
    lenPMMA = float(lenPMMA)
    lenZSA = float(lenZSA)
    nUnit = int(nUnit)
    
    # make the results dir
    rootDir = os.listdir('.')
    if 'results' not in rootDir:
        os.system('mkdir results')
    
    # make the sub result dir
    dirName = './results/pmma%.2f_zsa%.2f_units%i' % \
              (lenPMMA, lenZSA, nUnit)
    os.system('rm -r ' + dirName)
    os.system('mkdir ' + dirName)
    
    with open(dirName + '/run.mac', 'w') as f:
        # initialize
        s = """# Initialize kernel
/run/initialize
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

/run/beamOn %i
""" % history
        f.write(s)
    return dirName

def submit(dirName, pmma, zsa, n):
    # copy submit file to specific result dir
    os.system('cp ~/submit/submit.qsub ' + dirName)
    
    # write the run command
    with open(dirName + '/submit.qsub', 'a') as f:
        exe = os.path.abspath('../..') + '/build/exampleB1'
        f.write('\n%s run.mac %.6f %.6f %i > out.txt' % \
                (exe, pmma, zsa, n))
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
    history = 1e9
    # half inch and one inch layer detectors, in cm
    lengths = [0.5 * 2.54, 2.54, 4.0]
    # PMMA layer thickness in cm
    pmma = 0.02
    # ZnS(Ag) layer thickness in um
    zsa = 12.0
    
    for length in lengths:
        # run the append 4-cm length
        if length > 3.0:
            n = int(length / (pmma + zsa * 1e-4))
            adir = writer(pmma, zsa, n, history)
            submit(adir, pmma, zsa, n)

