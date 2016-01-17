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
    hams = np.diag(np.array([100.0,300.0,500.0]))*(1.0+0.0J)*cm2au

   #qmd1 = np.loadtxt('./qmod_0')
    qmd1 = np.zeros((3,3),dtype=float)
    qmd1[0,1] = qmd1[1,0] = 1.0/sqrt(2.0)
    qmd1[1,2] = qmd1[2,1] = 1.0
    qmds = np.array([np.dot(qmd1,qmd1)])*(1.0+0.0J)

    syst.init (ini['syst'],hams,qmds)
   
    # bath
    ini['bath']['temp'] = temp
    ini['bath']['mode'] = [
            {
                'jdru': [{'lamd':50.0*cm2au,'gamd':50.0*cm2au}],
                'jsdr': []
            }
        ]
    ini['bath']['nmod'] = len(ini['bath']['mode'])
    ini['bath']['pade'] = 2
    ini['bath']['npsd'] = 1
    bath.init (ini['bath'])
    
    # hidx
    ini['hidx']['lmax'] = 30
    ini['hidx']['nmax'] = 30000
    ini['hidx']['ferr'] = 2.0e-6

    #proprho
    jsonInit = {"deom":ini,
                "nonlinear":{
                    "w1max": 1200*cm2au,  # cm-1 convert to au
                    "w2max": 1200*cm2au,  # cm-1 convert to au
                    "dt":  1.0*fs2au,    # fs   convert to au
                    "nt1": 512,
                    "nt2": 512,
                    "nk": 16,
                    "staticErr": 2.0e-8,
                    "sch_hei": "hei",
                    "dip0": "inp_dip0.mat",
                    "dip1": "inp_dip1.mat",
                    "dip2": "inp_dip2.mat",
                }
            }
   #jsonInit = {"deom":ini,
   #            "linear":{
   #                "wmax": 1200*cm2au,  # cm-1 convert to au
   #                "dt":  1.0*fs2au,    # fs   convert to au
   #                "nt": 512,
   #                "nk": 16,
   #                "staticErr": 2.0e-9,
   #                "dip0": "inp_dip0.mat",
   #                "dip1": "inp_dip1.mat",
   #            }
   #        }
    dipo0 = (qmd1)*(1.0+0.0J)
    dipo1 = (qmd1)*(1.0+0.0J)
    dipo2 = np.dot(qmd1,qmd1)*(1.0+0.0J)
    arma.save(dipo0,"inp_dip0.mat")
    arma.save(dipo1,"inp_dip1.mat")
    arma.save(dipo2,"inp_dip2.mat")

    
    with open('input.json','w') as f:
        json.dump(jsonInit,f,indent=4) 
