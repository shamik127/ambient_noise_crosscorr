#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
import obspy.core as oc
import scipy.signal as ss
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import figure, show, rc
from mpl_toolkits.basemap import Basemap


def snr(dt, cclags, ccmat, umin, umax):

    offsets = map(lambda x: x.stats.distance*1e-3, ccmat.traces)
    rsd = np.array(offsets)

    tstart = rsd/umax
    tend = rsd/umin

    posbr_tl = cclags[len(cclags)//2:]
    
    pb_iws = np.searchsorted(posbr_tl,tstart)
    pb_iwe = np.searchsorted(posbr_tl,tend)

    egy_win_pb = np.empty(len(ccmat.traces))
    egy_owin_pb = np.empty(len(ccmat.traces))

    for tn in range(len(ccmat.traces)):
        posbr_cc = ccmat.traces[tn].data
        #print(posbr_cc[0:50])
        cccopy = posbr_cc.copy()
        pbwin = posbr_cc[pb_iws[tn]:pb_iwe[tn]+1]

        egy_win_pb[tn] = np.mean(np.square(pbwin))
        cccopy[pb_iws[tn]:pb_iwe[tn]+1] = 0.0
        egy_owin_pb[tn] = np.sum(np.square(cccopy))/( len(cccopy)- len(pbwin) )

    
    #print(egy_win_pb/egy_owin_pb)
    return np.mean(egy_win_pb/egy_owin_pb)


class do_one_brec:
    def __init__(self,ccmat,rlist,azs,dist,angles):
        a = 2
        self.azbin = azs
        self.distances = dist
        self.angles = angles
        self.rlist = rlist
        self.ccm = ccmat
        self.nrec = len(rlist)

        self.cum_arr = np.cumsum(np.arange( (self.nrec-1), 0, -1))
        self.cum_arr = np.insert(self.cum_arr, 0, 0)
        self.cum_arr = np.insert(self.cum_arr, -1, self.cum_arr[-1])


    def get_rec_coord(self, rec):
        x = crd[crd[:,0] == rec, 2]
        y = crd[crd[:,0] == rec, 3] 
        return x,y

    def getsection_with_azm(self,brec,azm):
        try:
            bri = list(self.rlist).index(brec)
        except ValueError:
            sys.exit('Could not find the base receiver.')

        upperRecs = np.arange(self.cum_arr[bri], self.cum_arr[bri+1])
        
        lowerRecs = np.array([], dtype=np.int32)
        for i in range(0,bri):
            lowerRecs = np.insert(lowerRecs,0, (self.cum_arr[bri-(i+1)]+i))


        

        # Base receiver angles with the upper receivers
        upperAngs = self.angles[upperRecs]

        # Base receiver angles with the lower receivers (the angles data contains angles corresponding to 
        # lower to upper receiver)
        lowerAngs = self.angles[lowerRecs] + 180

        # index positions in the angs array needed for the given angle
        upperAngIdsWithinBin = np.argwhere(np.abs(upperAngs-azm) < self.azbin/2.0).flatten()
        lowerAngIdsWithinBin = np.argwhere(np.abs(lowerAngs-azm) < self.azbin/2.0).flatten()

        upperRecsIdsWithinBin = upperRecs[upperAngIdsWithinBin]
        lowerRecsIdsWithinBin = lowerRecs[lowerAngIdsWithinBin]

        #return np.append( (lowerAngIdsWithinBin + 1), (upperAngIdsWithinBin + (bri+2)) )
        
        ccml = self.ccm.shape[0]
        upperTrs = [oc.trace.Trace(data=self.ccm[ccml//2+1:,x], header={'distance': self.distances[x]}) for x in upperRecsIdsWithinBin]
        lowerTrs = [oc.trace.Trace(data=np.flip(self.ccm[0:ccml//2-1,x]), header={'distance': self.distances[x]}) for x in lowerRecsIdsWithinBin]

        final_pnsep= oc.stream.Stream(traces=upperTrs) + oc.stream.Stream(traces=lowerTrs)
        #print(len(final_pnsep.traces))
        egy = snr(si,cclags,final_pnsep,1.5,6.0)
        
        return egy


def stn_pair_distances(rlist):
    

    nrecs = len(rlist)

    distances = np.array([])
    angles = np.array([])

    for bid in range(0, nrecs-1):
    
        base_x = crd.loc[bid].x
        base_y = crd.loc[bid].y

        recs = np.arange( (bid+1), nrecs ) 

        for rec in recs:
            rec_x = crd.loc[rec].x
            rec_y = crd.loc[rec].y

            dist = np.sqrt( (base_x - rec_x)**2 + (base_y - rec_y)**2 )
            qang = np.arctan2(rec_y - base_y, rec_x - base_x)
            
            if (qang < 0.0):
                qang += 2*np.pi

            rpazim_d = qang*180/np.pi
    
            distances = np.append(distances, dist)
            angles = np.append(angles, rpazim_d)

    return distances, angles


def recpair_angles(rlist):
    angles = {}

    nrecs = len(rlist)

    for bid in range(0, nrecs-1):
        angles[bid+1] = {}
        base_x = crd.loc[bid].x
        base_y = crd.loc[bid].y

        recs = np.arange(0, nrecs ) 
        recs = np.delete(recs, bid)

        for rec in recs:
            rec_x = crd.loc[rec].x
            rec_y = crd.loc[rec].y

            qang = np.arctan2(rec_y - base_y, rec_x - base_x)
            
            if (qang < 0.0):
                qang += 2*np.pi

            rpazim_d = qang*180/np.pi
            angles[bid+1][rec+1] = rpazim_d

    return angles

def recpair_distances(rlist):
    distances = {}

    nrecs = len(rlist)

    for bid in range(0, nrecs-1):
        distances[bid+1] = {}
        base_x = crd.loc[bid].x
        base_y = crd.loc[bid].y

        recs = np.arange(0, nrecs ) 
        recs = np.delete(recs, bid)

        for rec in recs:
            rec_x = crd.loc[rec].x
            rec_y = crd.loc[rec].y

            dist = np.sqrt( (base_x - rec_x)**2 + (base_y - rec_y)**2 )
            distances[bid+1][rec+1] = dist

    return distances


crd = pd.read_csv('coordinates_receivers_recnum.csv', sep='\t', header=None)
crd.columns = ['rec', 'lat', 'x', 'y']

loaded = np.load('stack_manypairs_100_zz.npz')
reclist=loaded['reclist']
fpband=loaded['fpband']
si=loaded['si'][0]
cclags=loaded['cclags']
cookie=loaded['cookie']
cclags=np.array(map(lambda x: float("%.2f" %x), cclags))

del loaded

# angles = recpair_angles(reclist)
# distances = recpair_distances(reclist)

# np.save('pairangles.npy', angles)
# np.save('pairdistances.npy', distances)

distances, angles = stn_pair_distances(reclist)

angs = np.arange(5, 360, 10)
egy = do_one_brec(cookie,reclist,10,distances, angles)

master_snr = np.empty([len(reclist), len(angs)])
for rec in range(len(reclist)):
    for j,ang in enumerate(angs):
        avg_snr = egy.getsection_with_azm(rec+1, ang)
        master_snr[rec][j] = avg_snr

np.savetxt('snr_data_100.txt', master_snr, delimiter=',')
  
# egy.getsection_with_azm(288, 5)


#print(angs[5])

# rec_dir = {}
# for rec in range(len(reclist)):
#     rec_dir[rec+1] = {}
#     for ang in angs:
#         rec_dir[rec+1][ang] = egy.getsection_with_azm(rec+1, ang)


# np.save('rec_rel.npy', rec_dir)
#print(rec_dir)
#egy.test_azm(3, 15)

#print(master_snr)
