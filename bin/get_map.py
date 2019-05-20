#!/usr/bin/env python
 
## this will get a fucking map the hard way
 
import astropy.io.fits as pf
import numpy as np
from numpy import sin,cos,exp
import healpy as hp
from scipy.sparse import dok_matrix
from scipy.sparse.linalg import lsqr, lsmr
import argparse ## SY 27/2/19
import configargparse ## SY 27/2/19

class Quasars:
    def __init__(self,Nside):
        self.Nside=Nside
        print ("Loading quasars")
        da=pf.open(zcat)[1] ## SY 19/2/19
        self.ra=da.data['RA']
        self.dec=da.data['DEC']
        self.phi=self.ra/180*np.pi
        self.theta=np.pi/2-self.dec/180*np.pi
        self.tid=da.data['THING_ID']
        self.ndx={}
        for i,t in enumerate(self.tid):
            self.ndx[t]=i
        print ("Healpix indices...")
        self.pixels=hp.ang2pix(Nside,self.theta,self.phi)
 
             
    def getCover(self,Nside):
        ### find covering pixels
        uniqpixels=np.array(list(set(self.pixels)))
        return uniqpixels
 
class Data:
    def __init__ (self,q,Nside):
        self.q=q
        self.Nside=Nside
 
        print ("Loading data...")
        da=pf.open(mapin)[1]       ## SY 18/2/19
        self.tid1=da.data['THID1']
        self.tid2=da.data['THID2']
        self.signal=da.data['SKAPPA']
        self.weight=da.data['WKAPPA']
        print ("Finding indices...")
        self.i1=np.array([q.ndx[tid] for tid in self.tid1])
        self.i2=np.array([q.ndx[tid] for tid in self.tid2])
        print ("Finding healpix indices")
        ## we will use idiotic midpoint.
        self.hi1=q.pixels[self.i1]
        self.hi2=q.pixels[self.i2]
                                  
         
    def len(self):
        return len(self.tid1)
     
class Analysis:
    def __init__ (self, Nside):
        #Nside=128
        print(srad, rtmax, rpmax)
        self.sradius=srad/180*np.pi ##search radius around map pixel ## SY 8/3/19
        self.pixarea=hp.nside2pixarea(Nside) ## SY 27/2/19
        self.Nside=Nside
        self.q=Quasars(Nside)
        self.d=Data(self.q,Nside)
 
 
        self.pixid=self.q.getCover(Nside)
        self.pixtheta,self.pixphi=hp.pix2ang(Nside,self.pixid)
        self.Np=len(self.pixid) ## number of pixels
        self.Nd=self.d.len()
 
 
        A=self.getAMatrix()
        b=self.d.signal*np.sqrt(self.d.weight) ## (we suppressed by weight)
        print ("running solver")
        mp=lsmr(A,b,atol=1e-06, btol=1e-06, show=True) ## SY 15/4/19
        #mp=lsmr(A,b,show=True)
        #mp=lsqr(A,b,show=True) ## SY 13/4/19
        print (mp[0])
        print (mp[1])
        print (mp[2])
        print (mp[3])
        print (mp[4])
        print (mp[5])
        print (mp[6])
        print (mp[7])
        nside  = np.array(Nside)
        pixels = self.pixid
        kappas = np.array(mp[0])
        if realn == 0:
            #np.savez('kappa_opt_srad/kappa_{}_rt{}_rp{}_nside{}_srad{}_{}'.format \
            np.savez('kappa_opt_srad/kappa_{}_rt{}_rp{}_nside{}_srad{}_{}'.format \
              (maptype, rtmax, rpmax, nside, srad, damping), nside, pixels, kappas)
        else:
            #np.savez('kappa_realisations/kappa_{}_rt{}_rp{}_nside{}_srad{}_{}'.format \
            #     (maptype, rtmax, rpmax, nside, srad, realn), nside, pixels, kappas)
            np.savez('est_opt_{}/kappa{}'.format(maptype, realn), nside, pixels, kappas)
           
         
    def getAMatrix(self):
        '''A is a sparse matrix consisting of response functions. It has
        dimensions #Healpixels by #QSO pairs.'''

        A = dok_matrix((self.Np, self.Nd), dtype=np.float32)
        sigma = damping * np.sqrt(4*np.pi/(12*self.Nside**2))
        print ("sigma=",sigma/np.pi*180,'deg')
         
        ## we loop over pixels
        for i,hpix in enumerate(self.pixid):
            ## first find nearby healpixels
            theta,phi=self.pixtheta[i],self.pixphi[i]
            mvec=(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
            neipixels=hp.query_disc(self.Nside,mvec ,self.sradius)
            assert(hpix in neipixels)
            B=np.zeros(self.Nd,type(False))
            for neipix in neipixels:
                B= B | (neipix==self.d.hi1) | (neipix==self.d.hi2)
            s=np.where(B)[0]
            ### so at this point we have for map pixel i, the list of data 
            ### pixels that are close enough to matter.
            ### we need to loop over them and get the relevant matrix elements
            qtheta1=self.q.theta[self.d.i1[s]]
            qphi1=self.q.phi[self.d.i1[s]]
            qtheta2=self.q.theta[self.d.i2[s]]
            qphi2=self.q.phi[self.d.i2[s]]
            ## we are going to employ 3 vectors as a foolproof method
            dx1=sin(qtheta1)*cos(qphi1)-mvec[0]
            dy1=sin(qtheta1)*sin(qphi1)-mvec[1]
            dz1=cos(qtheta1)-mvec[2]
            dr1=np.sqrt(dx1**2+dy1**2+dz1**2)
            ## we have response1 in direction dr1
            response1=1/dr1*(1-exp(-dr1**2/(2*sigma**2)))
            #response1=1/dr1
            ## ditto for q2
            dx2=sin(qtheta2)*cos(qphi2)-mvec[0]
            dy2=sin(qtheta2)*sin(qphi2)-mvec[1]
            dz2=cos(qtheta2)-mvec[2]
            dr2=np.sqrt(dx2**2+dy2**2+dz2**2)
            ## we have response2 in direction dr2
            response2=1/dr2*(1-exp(-dr2**2/(2*sigma**2)))
            #response2=1/dr2
            ## now we take the difference 
            dxr=dx1*response1-dx2*response2
            dyr=dy1*response1-dy2*response2
            dzr=dz1*response1-dz2*response2
            ## the difference in vector
            dx=dx1-dx2
            dy=dy1-dy2
            dz=dz1-dz2
            ## total response is movement/distance
            totresponse= np.sqrt(dxr*dxr+dyr*dyr+dzr*dzr)/np.sqrt(dx*dx+dy*dy+dz*dz) #SY 17/4/19
            #totresponse= (dxr*dx+dyr*dy+dzr*dz)/(dx*dx+dy*dy+dz*dz)
            totresponse*=np.sqrt(self.d.weight[s])  ## we downweigh response by weight
            A[i,s]=totresponse
            if (i%100==0):
                print (i)
 
        print ("Transposing matrix.")
        A=A.transpose()
        print ("A.tocsr")
        A=A.tocsr()
        return A
     
 
 
## main
parser = configargparse.ArgParser()

parser.add('--mapin', required=True, type=str, \
               help='pathname to kappalist fits file')
parser.add('--zcat', required=True,  type=str,\
               help='quasar catalogue name') 
parser.add('--maptype', required=True, type=str, \
               help='name of dataset (eg noisy)')
parser.add('--rtmax', required=True, type=int, default=40, \
               help='maximum transverse separation')
parser.add('--rpmax', required=True, type=int, default=10, \
               help='maximum radial separation')
parser.add('--nside', required=False, type=int, default=256, \
               help='resolution of map')
parser.add('--sradius', required=False, type=float, default=0.8, \
               help='search radius')
parser.add('--damping', required=False, type=float, default=1., \
               help='damping factor for sigma')
parser.add('--realisation', required=False, type=int, default=0, \
               help='realisation number')
args, unknown = parser.parse_known_args()

mapin   = args.mapin
zcat    = args.zcat
maptype = args.maptype
rtmax   = args.rtmax
rpmax   = args.rpmax
Nside   = args.nside
srad    = args.sradius
realn   = args.realisation
damping = args.damping

Analysis(Nside)
