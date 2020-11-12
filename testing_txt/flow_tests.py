import sys
import os
import re
import numpy as np

evtth = .1          # 100 eV for WFI Geant4 simulations
splitth = evtth    # same as evtth for WFI Geant4 simulations
xydep_min = -513
xydep_max = 513


# hash of codes indexed by particle type indexed
ptypes = {
     'proton': 0, 'gamma': 1, 'electron': 2, 'neutron': 3, 
     'pi+': 4, 'e+': 5, 'pi-': 6, 'nu_mu': 7, 
     'anti_nu_mu': 8, 'nu_e': 9, 'kaon+': 10, 'mu+': 11, 
     'deuteron': 12, 'kaon0L': 13, 'lambda': 14, 'kaon-': 15, 
     'mu-': 16, 'kaon0S': 17, 'alpha': 18, 'anti_proton': 19,
     'triton': 20, 'anti_neutron': 21, 'sigma-': 22,
     'sigma+': 23, 'He3': 24, 'anti_lambda': 25, 'anti_nu_e': 26,
     'anti_sigma-': 27, 'xi0': 28, 'anti_sigma+': 29, 'xi-': 30,
     'anti_xi0': 31, 'C12': 32, 'anti_xi-': 33, 'Li6': 34,
     'Al27': 35, 'O16': 36, 'Ne19': 37, 'Mg24': 38,
     'Li7': 39, 'He6': 40, 'Be8': 41, 'Be10': 42
}

def match(regex: str, string: str):
    return re.compile(regex).search(string)

def flow_test1(rf_name = 'flow_test1.txt'):
    with open(rf_name, 'r') as FH:
        f = FH.readlines()
    
    g = []
    with open(infile, 'r') as dec:
        for line in dec:
            if match('^\s*#', line): #skip comments 
                continue
            if match('^\s*$', line): #skip blank lines:
                continue
            if not match(',', line):
                continue
            g.append(line)
            fields = line.rstrip().split(',')
            if match('[a-zA-Z]', fields[0]):
                g.append("if\n")
            else:
                g.append("else\n")
                
    if len(f) != len(g):
        raise ValueError(f"The outputs had different number of lines {len(f)} vs {len(g)}")
        
    for i in range(len(g)):
        if f[i] != g[i]:
            raise ValueError(f"Different values found at index {i}: {f[i]}, {g[i]}")
            
    print("All good python")
    

def flow_test2(rf_name = 'flow_test2.txt'):
    with open(rf_name, 'r') as FH:
        f = FH.readlines()
    
    g = []
    with open(infile, 'r') as dec:
        for line in dec:
            if match('^\s*#', line): #skip comments
                continue
            if match('^\s*$', line): #skip blank lines:
                continue
            if not match(',', line):
                continue
            g.append(line)
            fields = line.rstrip().split(',')
            if match('[a-zA-Z]', fields[0]):
                g.append("if\n")
            else:
                if float(fields[2]) <= splitth:
                    continue 
                
                tmp_x, tmp_y = int(fields[0]), int(fields[1])
                if tmp_x<xydep_min or tmp_y<xydep_min or tmp_x>xydep_max or tmp_y>xydep_max:
                    continue # skip it if it's outside the 512x512 region of a quad
                g.append("else\n")
                
    if len(f) != len(g):
        raise TypeError
        
    for i in range(len(g)):
        if f[i] != g[i]:
            raise ValueError(f"Different values found at index {i}: {f[i]}, {g[i]}")
            
    print("All good python")
    
def flow_test3(rf_name = 'flow_test3.txt'):
    state = 0
    f = []
    perl_ptype = {}
    perl_cproc = {}
    with open(rf_name, 'r') as FH:
        for line in FH:
            #print(line, end = '')
            if line == "next\n":
                state += 1
                continue
            if state == 0:
                f.append(line)
            if state == 1:
                a = line.split()
                a[0] = a[0][1:]
                a[-1] = a[-1][:-1]
                perl_eid = np.array(a, dtype = int)
            if state == 2:
                key, value = line.split(',')
                perl_ptype[int(key)] = value.rstrip()
            if state == 3:
                key, value = line.split(',')
                perl_cproc[int(key)] = value.rstrip()
            if state == 4:
                deps = line.split(',')[:-1]
                perl_xdep = np.array(deps, dtype = int)
            if state == 5:
                deps = line.split(',')[:-1]
                perl_ydep = np.array(deps, dtype = int)
            if state == 6:
                deps = line.split(',')[:-1]
                perl_endep = np.array(deps, dtype = float)
            if state == 7:
                deps = line.split(',')[:-1]
                perl_rundep = np.array(deps, dtype = int)        
            if state == 8:
                deps = line.split(',')[:-1]
                perl_detectordep = np.array(deps, dtype = int)
            if state == 9:
                deps = line.split(',')[:-1]
                perl_eiddep = np.array(deps, dtype = int)
            if state == 10:
                deps = line.split(',')[:-1]
                perl_framedep = np.array(deps, dtype = int)
            if state == 11:
                deps = line.split(',')[:-1]
                perl_piddep = np.array(deps, dtype = int)
            if state == 12:
                deps = line.split(',')[:-1]
                perl_ptypedep = np.array(deps, dtype = int)
            if state == 13:
                ids = line.split(',')[:-1]
                perl_blobid = np.array(ids, dtype = int)
    
    eid = np.zeros(0, dtype=int)        # primary ID
    xdep = np.zeros(0, dtype=int)
    ydep = np.zeros(0, dtype=int)
    endep = np.zeros(0, dtype=float)
    rundep = np.zeros(0, dtype=int)
    detectordep = np.zeros(0, dtype=int)
    eiddep = np.zeros(0, dtype=int)
    framedep = np.zeros(0, dtype=int)
    piddep = np.zeros(0, dtype=int)
    ptypedep = np.zeros(0, dtype=int)
    blobid = np.zeros(0, dtype=int)

    ptype = {}
    cproc = {}
    rc = match('([0-9]+)_detector([0-9]+)', infile) #extracts the detector name
    this_run = int(rc.group(1))
    this_detector = int(rc.group(2))  
     
    g = []
    with open(infile, 'r') as IN:
        for line in IN: #switched to a for loop because of built-in __iter__ method
            if match('^\s*#', line): #skip comments #added ability to have arbritrary whitespace before '#'
                continue
            if match('^\s*$', line): #skip blank lines:
                continue
            if not match(',', line): #could be if ',' not in line
                continue

            g.append(line)
            fields = line.rstrip().split(',')

            if match('[a-zA-Z]', fields[0]):
                this_eid = int(fields[1])
                eid = np.append(eid, int(this_eid))
                ptype[int(fields[2])] = fields[0]
                cproc[int(fields[2])] = fields[4]
                g.append("if\n")
            else:
                if float(fields[2]) <= splitth:
                    continue

                tmp_x, tmp_y = int(fields[0]), int(fields[1])
                if tmp_x<xydep_min or tmp_y<xydep_min or tmp_x>xydep_max or tmp_y>xydep_max:
                    continue # skip it if it's outside the 512x512 region of a quad

                xdep = np.append(xdep, tmp_x)
                ydep = np.append(ydep, tmp_y)
                endep = np.append(endep, float(fields[2]))
                rundep = np.append(rundep, this_run)
                detectordep = np.append(detectordep, this_detector)
                eiddep = np.append(eiddep, this_eid)
                framedep = np.append(framedep, 0)
                piddep = np.append(piddep, int(fields[3]))
                ptypedep = np.append(ptypedep, ptypes[ptype[int(fields[3])]])
                blobid = np.append(blobid, 0)#np.array(0, dtype = int));
                g.append("else\n")
                
    to_test = {'control': [f, g], 'eid': [perl_eid, eid], 
               'ptype': [perl_ptype, ptype], 'cproc': [perl_cproc, cproc],
               'xdep': [perl_xdep, xdep], 'ydep': [perl_ydep, ydep],
               'endep': [perl_endep, endep], 'rundep': [perl_rundep, rundep],
               'detectordep': [perl_detectordep, detectordep], 'eiddep': [perl_eiddep, eiddep],
               'framedep': [perl_framedep, framedep], 'piddep': [perl_piddep, piddep],
               'ptypedep': [perl_ptypedep, ptypedep], 'blobid': [perl_blobid, blobid]}
                
    for test, items in to_test.items():        
        if len(items[0]) != len(items[1]):
            raise ValueError(f"Different lenghts of {test}: {len(items[0])} vs {len(items[1])}")

        if isinstance(items[0], dict):
            for key in items[0]:
                if items[0][key] != items[1][key]:
                    raise ValueError(f"Different values of {test} found at key {key}: {items[0][key]}, {items[1][key]}")
        else:
            for i in range(len(items[0])):
                if items[0][i] != items[1][i]:
                    raise ValueError(f"Different values of {test} found at index {i}: {items[0][i]}, {items[1][i]}")

    
    print("All good python")
    
do_tests = sys.argv[1]
inp_file = sys.argv[2]
cwd = os.getcwd()
input_dir = cwd[:cwd.rindex('/') + 1]
infile =  input_dir + 'input/' + inp_file; 

eval(f"flow_test{do_tests}()")
