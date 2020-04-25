#!/usr/bin/python

import os
import sys
import numpy as np
import obspy.core as oc
import scipy.signal as ss
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import figure, show, rc
from mpl_toolkits.basemap import Basemap


def snr(dt, cclags, ccmat, umin, umax):

    offsets = map(lambda x: x.stats.distance*1e-3, ccmat.traces)
    rsd = np.array(offsets)
    
    hlen = (len(cclags)-1)/2

    tstart = rsd/umax
    tend = rsd/umin

    pb_iws = np.searchsorted(cclags,tstart)
    pb_iwe = np.searchsorted(cclags,tend)

    egy_win_pb = np.empty(len(ccmat.traces))
    egy_owin_pb = np.empty(len(ccmat.traces))


    for tn in range(len(ccmat.traces)):
        cccopy = ccmat.traces[tn].data
        pbwin = ccmat.traces[tn].data[pb_iws[tn]:pb_iwe[tn]+1]

        egy_win_pb[tn] = np.sum(np.square(pbwin))
        cccopy[pb_iws[tn]:pb_iwe[tn]+1] = 0.0
        egy_owin_pb[tn] = np.sum(np.square(cccopy))

    return np.mean(egy_win_pb/egy_owin_pb)


class do_one_brec:
    def __init__(self,ccmat,rlist,azs,dist,angles):
        a = 2
        self.azbin = azs
        self.distances = dist
        self.angles = angles
        self.rlist = rlist
        self.ccm = ccmat


    def get_rec_coord(self, rec):
        x = crd[crd[:,0] == rec, 2]
        y = crd[crd[:,0] == rec, 3] 
        return x,y


    def getsection_with_azm(self,brec,azm):
        try:
            bri = list(self.rlist).index(brec)
        except ValueError:
            sys.exit('Could not find the base receiver.')

        base_x, base_y = self.get_rec_coord(bri)

        nrec = len(self.rlist)

        cum_arr = np.cumsum(np.arange( (nrec-1), 0, -1))
        cum_arr = np.insert(cum_arr, 0, 0)

        rel_rec_idx = np.arange(cum_arr[bri], cum_arr[bri+1])

        angs = self.angles[rel_rec_idx]


        ang_args = np.argwhere(np.abs(angs) < self.azbin/2)

        if azm==0 or azm==360:
            set1=np.concatenate( (np.argwhere(np.abs(angs-azm) < self.azbin/2), np.argwhere(np.abs(360-angs-azm) < self.azbin/2)), axis=0)[:,0]
        else:
            set1=np.argwhere(np.abs(angs-azm) < self.azbin/2)[:,0]

        set2=np.concatenate((np.argwhere(np.abs(angs-180-azm) < self.azbin/2), np.argwhere(np.abs(angs+180-azm) < self.azbin/2)), axis=0)[:,0]

        ang_mod_idx_set1 = rel_rec_idx[set1]
        ang_mod_idx_set2 = rel_rec_idx[set2]

        trs1 = [oc.trace.Trace(data=self.ccm[:,x], header={'distance': self.distances[x]}) for x in ang_mod_idx_set1]
        trs2 = [oc.trace.Trace(data=self.ccm[:,x], header={'distance': self.distances[x]}) for x in ang_mod_idx_set2]

        final_pnsep=oc.stream.Stream(traces=trs1) + oc.stream.Stream(traces=trs2)
        egy = snr(si,cclags,final_pnsep,1.5,6.0)
        
        return egy


def stn_pair_distances(rlist):
    
    nrecs = len(rlist)
    cordfl = np.genfromtxt('coordinates_receivers_recnum.csv')

    distances = np.array([])
    angles = np.array([])

    for bidx in range(1, nrecs):
        base_x = cordfl[cordfl[:,0] == bidx, 2]
        base_y = cordfl[cordfl[:,0] == bidx, 3]

        for rec in range((bidx + 1), (nrecs+1)):
            rec_x = cordfl[cordfl[:,0] == rec, 2]
            rec_y = cordfl[cordfl[:,0] == rec, 3]

            dist = np.sqrt( (base_x - rec_x)**2 + (base_y - rec_y)**2 )
            qang = np.arctan2(rec_x - base_x, rec_y - base_y)

            qang[qang < 0.0] += 2*np.pi
            rpazim_d = qang*180/np.pi

            distances = np.append(distances, dist)
            angles = np.append(angles, rpazim_d)

    return distances, angles

"""
class geo_plot:
    def __init__(self):
        buildingdf = np.genfromtxt('coordinates_rec_latlon_id_gcy.csv', delimiter='\t')
        self.lat  = buildingdf[:,0]
        self.lon  = buildingdf[:,1]

        # determine range to print based on min, max lat and lon of the data
        margin = 2 # buffer to add to the range
        lat_min = min(self.lat) - margin
        lat_max = max(self.lat) + margin
        lon_min = min(self.lon) - margin
        lon_max = max(self.lon) + margin


        print(lat_min, lat_max, lon_min, lon_max)

        # create map using BASEMAP
        self.m = Basemap(llcrnrlon=lon_min,
                    llcrnrlat=lat_min,
                    urcrnrlon=lon_max,
                    urcrnrlat=lat_max,
                    lat_0=(lat_max - lat_min)/2,
                    lon_0=(lon_max-lon_min)/2,
                    projection='merc',
                    resolution = 'h',
                    area_thresh=10000.,
                    )
        self.m.drawcoastlines()
        self.m.drawcountries()
        self.m.drawstates()
        self.m.drawmapboundary(fill_color='#46bcec')
        self.m.fillcontinents(color = 'white',lake_color='#46bcec')


        self.m.drawparallels(np.arange(lat_min,lat_max,0.5),labels=[1,1,0,0])
        self.m.drawmeridians(np.arange(lon_min,lon_max,0.8),labels=[0,0,0,1])
        plt.title('Receivers location')

    def draw(self):
        for lo, la in zip(self.lon, self.lat):
            lons, lats = self.m(lo, la)
            self.m.scatter(lons, lats, marker = 'o', color='r', zorder=5)
            plt.savefig("rec_on_map", dpi=96)
            plt.show()

geo = geo_plot()
geo.draw()
"""

crd = np.genfromtxt('coordinates_receivers_recnum.csv')
loaded = np.load('stack_manypairs_99_zz.npz')
reclist=loaded['reclist']
fpband=loaded['fpband']
si=loaded['si'][0]
cclags=loaded['cclags']
cookie=loaded['cookie']
cclags=np.array(map(lambda x: float("%.2f" %x), cclags))

distances, angles = stn_pair_distances(reclist)
print(angles)

angs = np.arange(0, 370, 30)
egy = do_one_brec(cookie,reclist,5,distances, angles)

master_snr = np.empty([len(reclist), len(angs)])
for i in range(len(reclist)-1):
    for j,ang in enumerate(angs):
        avg_snr = egy.getsection_with_azm(i+1, ang)
        master_snr[i][j] = avg_snr


print(master_snr[288,:])

np.savetxt('snr_data.txt', master_snr, delimiter=',')

"""
fig = figure(figsize=(8,8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

N = len(angs)
theta = np.arange(0.0, 2*np.pi, 2*np.pi/N)
radii = [10]*N
width = [np.pi/10]*N
bars = ax.bar(theta, radii, width=width, bottom=0.0)
for r,bar,snr in zip(radii, bars, master_snr[21]):
    bar.set_alpha(0.5)
    if snr == np.nan or snr == 0.0:
        snr = 0
        bar.set_alpha(0.0)

    bar.set_facecolor( cm.jet(snr))
    

show()
"""