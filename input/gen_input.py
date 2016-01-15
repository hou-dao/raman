import sys
sys.path.append('../scripts')

from math import sqrt
import json
import numpy as np
import armadillo as arma
import syst
import bath
from unitconverter import cm2au, fs2au, kt2au

if __name__ == '__main__':

    with open('default.json') as f:
        ini = json.load(f)

    temp = 298.0*kt2au
    
    # syst
    mat1 = np.random.rand(4,4)
    mat2 = np.random.rand(4,4)

    hams = (mat1+mat1.T)*(1.0+0.0J)*temp

    qmat = (mat2+mat2.T)*0.5
    qmds = np.array([np.dot(qmat,qmat)])*(1.0+0.0J)

    syst.init (ini['syst'],hams,qmds)
   
    # bath
    ini['bath']['temp'] = temp
   #ini['bath']['jsdr'] = [{'lams':1.0,'omgs':5.0,'gams':1.0}]
    ini['bath']['jdru'] = [{'lamd':50.0*cm2au,'gamd':100.0*cm2au}]
    ini['bath']['nmod'] = qmds.shape[0]
    ini['bath']['pade'] = 2
    ini['bath']['npsd'] = 1
    bath.init (ini['bath'])
    
    # hidx
    ini['hidx']['lmax'] = 15
    ini['hidx']['nmax'] = 30000
    ini['hidx']['ferr'] = 2.0e-6

    #proprho
    jsonInit = {"deom":ini,
                "nonlinear":{
                    "w1max": 800*cm2au,  # cm-1 convert to au
                    "w2max": 800*cm2au,  # cm-1 convert to au
                    "dt":  1.0*fs2au,    # fs   convert to au
                    "nt1": 512,
                    "nt2": 512,
                    "nk": 16,
                    "staticErr": 2.0e-9,
                    "sch_hei": "hei",
                    "dip0": "inp_dip0.mat",
                    "dip1": "inp_dip1.mat",
                    "dip2": "inp_dip2.mat",
                }
            }
    dipo0 = qmat*(1.0+0.0J)
    dipo1 = np.dot(qmat,qmat)*(1.0+0.0J)
    dipo2 = qmat*(1.0+0.0J)
    arma.save(dipo0,"inp_dip0.mat")
    arma.save(dipo1,"inp_dip1.mat")
    arma.save(dipo2,"inp_dip2.mat")

    
    with open('input.json','w') as f:
        json.dump(jsonInit,f,indent=4) 
