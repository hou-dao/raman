from math import sqrt
import numpy as np
from unitconvert import fs2au, cm2au, kt2au
import syst
import bath

class hidx:

    def __init__ (self, nind, nmax):

        self.keys = np.array(['!'*nind if x>0 else '0'*nind for x in xrange(nmax)])
        self.locs = {'0'*nind:0,}
        self.gams = np.zeros(nmax,dtype=float)
        self.levs = np.array([-9527 if x>0 else 0 for x in xrange(nmax)])
        self.nrho = 1
        self.tier = 0


class deom:

    def __init__ (self, InitJson):

        self.type = InitJson['proj']['type']
        self.logFile = InitJson['proj']['logFile']
        self.outFile = InitJson['proj']['outFile']
        self.nsys =   int(InitJson['syst']['nsys'])
        self.nmod =   int(InitJson['syst']['nmod'])
        self.npsd =   int(InitJson['bath']['npsd'])
        self.iniSite = int(InitJson['prop']['iniSite'])
        self.nmax =   int(InitJson['prop']['nmax'])
        self.lmax =   int(InitJson['prop']['lmax'])
        self.tolr = float(InitJson['prop']['toler'])
        self.nstep=   int(InitJson['prop']['nstep'])
        self.nchk =   int(InitJson['prop']['nchk'])
        self.dt   = float(InitJson['prop']['dt'])
        self.tem  = float(InitJson['bath']['tem'])
        self.nper = self.npsd+1
        if self.temUnit == 'K':
            self.tem *= kt2au
        elif self.temUnit == 'cm-1':
            self.tem *= cm2au
        self.nind = self.nmod*self.nper
        self.hams, self.qmds = syst.init()
        self.etal, self.etar, self.etaa, self.expn, self.wres = bath.init()

    def propagation (self):

        indx = hidx(self.nind,self.nmax)
        rhot = np.zeros((self.nmax,self.nsys,self.nsys),dtype=complex)
        rhot[0,self.iniSite-1,self.iniSite-1] = 1.0+0.0J

        dt = self.dt 
        self.dt *= fs2au

        flog = open(self.logFile,'a')
        fout = open(self.outFile,'w')
        for k in xrange(self.nstep):
            if k%self.nchk == 0:
                print >> flog, 'time step:', k, 'nrho =', indx.nrho, 'tier =', indx.tier
                print >> fout, '%16.6e'%(k*dt), ''.join('%16.6e'%x for x in rhot[0].diagonal().real)
            self.rk4(rhot,indx,k,self.dt)
        fout.close()
        flog.close()


    def filter (self,rhos,indx):
        nreal, ntier = 0, 0
        for iado in xrange(1,indx.nrho):
            rho, key, gam, lev = rhos[iado], indx.keys[iado], indx.gams[iado], indx.levs[iado]
            if np.any(np.abs(rho)>self.tolr):
                nreal += 1
                if nreal != iado:
                    indx.locs[key] = nreal
                    rhos[nreal], indx.keys[nreal], indx.gams[nreal], indx.levs[nreal] = rho, key, gam, lev
                ntier = max(lev,ntier)
            else:
                del indx.locs[key]
        indx.nrho, indx.tier = nreal+1, ntier


    def rk4 (self,rhot,indx,k,dt):

        n0 = indx.nrho
        dt2,dt6 = dt*0.5,dt/6.0;

        rho1 = np.empty(rhot.shape,dtype=complex)
        rho2 = np.empty(rhot.shape,dtype=complex)
        rho3 = np.empty(rhot.shape,dtype=complex)

        n1 = self.rem(rho1,rhot,indx,k*dt)
        rho3[:n0] = rhot[:n0]+rho1[:n0]*dt2
        if n1>n0:
            rho3[n0:n1] = rho1[n0:n1]*dt2

        n2 = self.rem(rho2,rho3,indx,(k+0.5)*dt)
        rho1[:n1] += rho2[:n1]*2
        if n2>n1:
            rho1[n1:n2] = rho2[n1:n2]*2
        rho3[:n0] = rhot[:n0]+rho2[:n0]*dt2
        if n2>n0:
            rho3[n0:n2] = rho2[n0:n2]*dt2

        n3 = self.rem(rho2,rho3,indx,(k+0.5)*dt)
        rho1[:n2] += rho2[:n2]*2
        if n3>n2:
            rho1[n2:n3] = rho2[n2:n3]*2
        rho3[:n0] = rhot[:n0]+rho2[:n0]*dt
        if n3>n0:
            rho3[n0:n3] = rho2[n0:n3]*dt

        n4 = self.rem(rho2,rho3,indx,(k+1.0)*dt)
        rho1[:n3] += rho2[:n3]
        if n4>n3:
            rho1[n3:n4] = rho2[n3:n4]
        rhot[:n0] = rhot[:n0]+rho1[:n0]*dt6
        if n4>n0:
            rhot[n0:n4] = rho2[n0:n4]*dt6

        self.filter(rhot,indx)


    def rem (self,drho,rhot,indx,t):

        nsys,nmod,nper,nind,nmax,lmax,nado = self.nsys,self.nmod,self.nper,self.nind,self.nmax,self.lmax,indx.nrho
        hams,qmds = self.hams,self.qmds
        etal,etar,etaa,expn,wres = self.etal,self.etar,self.etaa,self.expn.flat,self.wres

        drho[:indx.nrho] = np.array([-1.J*(np.dot(hams,rhot[i])-np.dot(rhot[i],hams))\
                -indx.gams[i]*rhot[i] for i in xrange(indx.nrho)])
       
        for irho in xrange(indx.nrho):

            if np.any(np.abs(rhot[irho])>self.tolr):

                rhoi = rhot[irho]
                ikey, damp, tier = indx.keys[irho], indx.gams[irho], indx.levs[irho]

                qact = np.array([[np.dot(q,rhoi),np.dot(rhoi,q)] for q in qmds])
                drho[irho] += -np.sum(wres[m]*(np.dot(qmds[m],qact[m,0]-qact[m,1]) \
                              -np.dot(qact[m,0]-qact[m,1],qmds[m])) for m in xrange(nmod))

                for mp in xrange(nind):
                    m,p,n = mp/nper,mp%nper,ord(ikey[mp])-48
                    el,er,ea = etal[m,p],etar[m,p],etaa[m,p]
                    if tier<lmax:
                        sca = -1.0J*sqrt((n+1.0)/ea)
                        cjl,cjr = sca*el,sca*er
                        key = ikey[:mp]+chr(n+49)+ikey[mp+1:]
                        try:
                            loc = indx.locs[key]
                            drho[loc] += cjl*qact[m,0]-cjr*qact[m,1]
                        except:
                            if nado+1>nmax:
                                raise ValueError('Hierarchy space overflow!')
                                sys.exit(1)
                            indx.locs[key] = nado
                            indx.keys[nado], indx.gams[nado], indx.levs[nado] = key, damp+expn[mp], tier+1
                            drho[nado] = cjl*qact[m,0]-cjr*qact[m,1]
                            nado += 1
                    if n>0: 
                        sca = -1.0J*sqrt(n*ea)
                        key = ikey[:mp]+chr(n+47)+ikey[mp+1:]
                        try:
                            loc = indx.locs[key]
                            drho[loc] += sca*(qact[m,0]-qact[m,1])
                        except:
                            if nado+1>nmax:
                                raise ValueError('Hierarchy space overflow!')
                                sys.exit(1)
                            indx.locs[key] = nado
                            indx.keys[nado], indx.gams[nado], indx.levs[nado] = key, damp+expn[mp], tier-1
                            drho[nado] = sca*(qact[m,0]-qact[m,1])
                            nado += 1
        indx.nrho = nado
        return nado
