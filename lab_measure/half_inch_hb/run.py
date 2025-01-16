import numpy as np
import os
import socket

def writer(history):    
    # result dir
    dirName = './results/cf252'
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
/run/beamOn %i""" % history
        f.write(s)
    # eliminate the preceeding '.'
    return dirName

def submit(dir, absLen, sa):
    os.system('cp ~/submit/submit.qsub ' + dir + '/submit.qsub')
    with open(dir + '/submit.qsub', 'a') as f:
        f.write(os.path.abspath('.') + '/build/exampleB1 run.mac ' + 
                    str(absLen) + ' ' + str(sa) + ' > out.txt')
    os.system('cd ' + dir + " && qsub submit.qsub")
    return

def submit_local(dir, absLen, sa):
    command = 'cd DIR && EXE run.mac ABSLEN SA > out.txt'
    command = command.replace('DIR', dir)
    exe = os.path.abspath('.') + '/build/exampleB1'
    command = command.replace('EXE', exe)
    command = command.replace('ABSLEN', '%.6f' % absLen)
    command = command.replace('SA', '%6f' % sa)
    os.system(command)
    
    

if __name__ == '__main__':
    rootDir = os.listdir('.')
    if 'results' not in rootDir:
        os.system('mkdir results')
    
    # mean free path of optical photon in ZnS(Ag)
    frac = 0.01
    absLen = -61. / np.log(frac)
    sa = 1.
    
    # auto determine running on office pc (for test run)
    # or super computer (real run)
    if socket.gethostname() == 'corps-g-1':
        history = 100
        aDir = writer(history)
        submit_local(aDir, absLen, sa)
        
    else:
        history = 1e9
        aDir = writer(history)
        submit(aDir, absLen, sa)

    
    
