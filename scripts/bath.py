#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt
from BoseFermiExpansion import PSD
import armadillo as arma 

def jwdru (omg,jdru):
    lamd, gamd = jdru['lamd'], jdru['gamd']
    return 2.0*lamd*gamd*omg/(omg**2+gamd**2)

def jwsdr (omg,jsdr):
    lams, omgs, gams = jsdr['lams'], jsdr['omgs'], jsdr['gams']
    return 2.0*lams*omgs**2*gams*omg/((omg**2-omgs**2)**2+(omg*gams)**2)

def fBose (x,pole,resi,rn,tn):
    return 1/x+0.5+rn*x+tn*x**3+sum(2.0*resi[i]*x/(x**2+pole[i]**2) for i in xrange(len(pole)))
    
def init (inidic):
        
    nmod = inidic['nmod']
    npsd = inidic['npsd']
    pade = inidic['pade']
    temp = inidic['temp']

    ndru = max([len(m['jdru']) for m in inidic['mode']])
    nsdr = max([len(m['jsdr']) for m in inidic['mode']])
    nper = ndru+2*nsdr+npsd

    delr = np.zeros(nmod,dtype=float)
    expn = np.zeros((nmod,nper),dtype=complex)
    etal = np.zeros((nmod,nper),dtype=complex)
    etar = np.zeros((nmod,nper),dtype=complex)
    etaa = np.zeros((nmod,nper),dtype=float)

    pole, resi, rn, tn = PSD (npsd,BoseFermi=1,pade=pade)
    
    for m in xrange(nmod):

        jdru = inidic['mode'][m]['jdru']
    
        for idru in xrange(len(jdru)):
            lamd, gamd = jdru[idru]['lamd'], jdru[idru]['gamd']
            iper = idru
            expn[m,iper] = gamd
            etal[m,iper] = -2.J*lamd*gamd*fBose(-1.J*gamd/temp,pole,resi,rn,tn)
            etar[m,iper] = etal[m,iper].conj()
            etaa[m,iper] = abs(etal[m,iper])
            delr[m] += 2.*lamd*gamd/temp*rn
        
        jsdr = inidic['mode'][m]['jsdr']

        for isdr in xrange(len(jsdr)):
            lams, omgs, gams = jsdr[isdr]['lams'], jsdr[isdr]['omgs'], jsdr[isdr]['gams']
            iper, jper = ndru+isdr*2, ndru+isdr*2+1
            etaBO = 2.*lams*omgs*omgs*gams
            Delta = omgs*omgs-gams*gams/4.0
            delr[m] +=  etaBO*tn/(temp**3)
            if Delta > 0:
                OmgB = sqrt(Delta)
                expn[m,iper] = 0.5*gams+1.J*OmgB
                expn[m,jper] = 0.5*gams-1.J*OmgB
            elif Delta < 0:
                OmgB = sqrt(-Delta)
                expn[m,iper] = 0.5*gams+OmgB
                expn[m,jper] = 0.5*gams-OmgB           
            else:
                raise ValueError("Not prepared for Delta=0")
            z1, z2 = -1.J*expn[m,iper], -1.J*expn[m,jper]
            etal[m,iper] = -2.J*etaBO*z1/(2.*z1*(z1+z2)*(z1-z2))*fBose(z1/temp,pole,resi,rn,tn)
            etal[m,jper] = -2.J*etaBO*z2/(2.*z2*(z2+z1)*(z2-z1))*fBose(z2/temp,pole,resi,rn,tn)
            if Delta > 0:
                etar[m,iper] = etal[m,jper].conj()
                etar[m,jper] = etal[m,iper].conj()
                etaa[m,iper] = etaa[m,jper] = sqrt(abs(etal[m,iper])*abs(etal[m,jper]))
            elif Delta < 0:
                etar[m,iper] = etal[m,iper].conj()
                etar[m,jper] = etal[m,jper].conj()
                etaa[m,iper] = abs(etal[m,iper])
                etaa[m,jper] = abs(etal[m,jper])
            else:
                raise ValueError("Not prepared for Delta=0")

        for ipsd in xrange(npsd):
        
            iper = ndru+nsdr*2+ipsd
            zomg = -1.J*pole[ipsd]*temp
            jsum = sum(jwdru(zomg,x) for x in jdru)+sum(jwsdr(zomg,x) for x in jsdr)
            expn[m,iper] = pole[ipsd]*temp            
            etal[m,iper] = -2.J*resi[ipsd]*temp*jsum   
            etar[m,iper] = etal[m,iper].conj()
            etaa[m,iper] = abs(etal[m,iper])  
                            
    arma.save (etal,inidic['etalFile'])
    arma.save (etar,inidic['etarFile'])
    arma.save (etaa,inidic['etaaFile'])
    arma.save (expn,inidic['expnFile'])
    arma.save (delr,inidic['delrFile'])
            
            
if __name__ == '__main__':

    inidic = {
        "nmod": 1,
        "npsd": 0,
        "pade": 2,
        "temp": 1.0,
        "jdru": [{"lamd":1.0,"gamd":1.0}],
        "jsdr": [],
	"etalFile": "inp_etal.mat",
	"etarFile": "inp_etar.mat",
	"etaaFile": "inp_etaa.mat",
	"expnFile": "inp_expn.mat",
	"delrFile": "inp_delr.mat"
    }
    init(inidic)
    
