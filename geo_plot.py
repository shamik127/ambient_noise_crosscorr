#!/usr/bin/python

from os import path, makedirs
import sys
import numpy as np
from collections import OrderedDict
import obspy.core as oc
import pandas as pd
#import scipy.signal as ss
#import matplotlib.cm as cm
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import figure, show, rc

def snr(dt, cclags, ccmat, umin, umax, rel_recs):

    offsets = map(lambda x: x.stats.distance*1e-3, ccmat.traces)
    rsd = np.array(offsets)   
    tstart = rsd/umax
    tend = rsd/umin
    posbr_tl = cclags[len(cclags)//2+1:]
    pb_iws = np.searchsorted(posbr_tl,tstart)
    pb_iwe = np.searchsorted(posbr_tl,tend)
    egy_win_pb = np.empty(len(ccmat.traces))
    egy_owin_pb = np.empty(len(ccmat.traces))

    snrall={}

    for tn in range(len(ccmat.traces)):
        posbr_cc = ccmat.traces[tn].data
        cccopy = posbr_cc.copy()
        pbwin = posbr_cc[pb_iws[tn]:pb_iwe[tn]+1]
        egy_win_pb[tn] = np.mean(np.square(pbwin))
        cccopy[pb_iws[tn]:pb_iwe[tn]+1] = 0.0
        egy_owin_pb[tn] = np.sum(np.square(cccopy))/( len(cccopy) - len(pbwin) )
        #snrall[rel_recs[tn]] = egy_win_pb[tn]/egy_owin_pb[tn] * rsd[tn]/30 * kernel[tn]
        snrall[rel_recs[tn]] = egy_win_pb[tn]/egy_owin_pb[tn]

    if (gaussw):
        kernel = gaussian_kernel(rsd,gaussw) 
        for i,key in enumerate(snrall.keys()):
            snrall[key] *= kernel[i]

    if (linw):
        for i,key in enumerate(snrall.keys()):
            snrall[key] *= rsd[i]/30

    return snrall


class do_one_brec:
    def __init__(self,ccmat,rlist,azms,azbin):
        a = 2
        self.azms=azms
        self.azbin = azbin
        self.rlist = rlist
        self.ccm = ccmat
        self.nrec = len(reclist)

        self.cum_arr = np.cumsum(np.arange( (self.nrec-1), 0, -1))
        self.cum_arr = np.insert(self.cum_arr, 0, 0)
        self.cum_arr = np.insert(self.cum_arr, -1, self.cum_arr[-1])

        self.master_snr=np.empty([len(self.rlist), len(self.azms)])
        self.snrall={}
        self.avgsnrdict={}
        self.avgsnr=np.empty([len(self.rlist), len(self.azms)])
        self.get_snr()


    def get_snr(self):
        print('\nCalculating snr...')
        for i,brec in enumerate(self.rlist):
            #print(brec)
            try:
                bri = list(reclist).index(brec)
            except ValueError:
                sys.exit('Could not find the base receiver.')

            self.snrall[brec]={}

            upperRecs = np.arange(self.cum_arr[bri], self.cum_arr[bri+1])
            lowerRecs = np.array([], dtype=np.int32)
            for i in range(0,bri):
                lowerRecs = np.insert(lowerRecs,0, (self.cum_arr[bri-(i+1)]+i))
            upperAngs = angles[upperRecs]            
            lowerAngs = angles[lowerRecs] + 180

            for i in range(len(lowerAngs)):
                if (lowerAngs[i] > 360):
                    lowerAngs[i] = lowerAngs[i] - 360
            #print(azms)
            for j,azm in enumerate(self.azms):
                upperAngIdsWithinBin = np.argwhere(np.abs(upperAngs-azm) < self.azbin/2.0).flatten()
                lowerAngIdsWithinBin = np.argwhere(np.abs(lowerAngs-azm) < self.azbin/2.0).flatten()

                upperRecsIdsWithinBin = upperRecs[upperAngIdsWithinBin]
                lowerRecsIdsWithinBin = lowerRecs[lowerAngIdsWithinBin]

                #print(lowerAngIdsWithinBin + 1)
                #print(upperAngIdsWithinBin + (bri+2))
                rel_recs = np.append( (lowerAngIdsWithinBin + 1), (upperAngIdsWithinBin + (bri+2)) )
                #print(rel_recs)
                # print(rel_recs)
                #print('--------------------')
                
                ccml = self.ccm.shape[0]
                upperTrs = [oc.trace.Trace(data=self.ccm[ccml//2+1:,x], header={'distance': distances[x]}) for x in upperRecsIdsWithinBin]
                lowerTrs = [oc.trace.Trace(data=np.flip(self.ccm)[ccml//2+1:,x], header={'distance': distances[x]}) for x in lowerRecsIdsWithinBin]
                final_pnsep= oc.stream.Stream(traces=lowerTrs) + oc.stream.Stream(traces=upperTrs)
            
                snrall = snr(si,cclags,final_pnsep,1.5,6.0, rel_recs)
                self.snrall[brec][azm]=OrderedDict(sorted(snrall.items()))
            
        self.avg_snr()
        print('\nDone!')


    def get_all_snr(self):
        return self.snrall


    def avg_snr(self):
        for i,rec in enumerate(self.rlist):
            self.avgsnrdict[rec]={}
            for j,azm in enumerate(self.azms):
                avgsnrval=np.mean(self.snrall[rec][azm].values())
                self.avgsnrdict[rec][azm]=avgsnrval
                self.avgsnr[i][j]=avgsnrval

        self.select_percentile()


    def select_percentile(self):
        num=np.nanpercentile(self.avgsnr, 90)
        if not path.exists('percentile.npy'):
            perc_dir={}
        else:
            perc_dir=np.load('percentile.npy', allow_pickle=True).item()

        perc_dir[name]=num
        np.save('percentile.npy', perc_dir)
        

    def get_avgsnrdict(self):
        return self.avgsnrdict


    def get_avgsnr(self):
        return self.avgsnr


    def save_files(self):
        np.savetxt('Results/%s/avgsnr.txt' %(name), self.avgsnr, delimiter=',')
        np.save('Results/%s/avgsnrdict.npy' %(name), self.avgsnrdict)



########################################################################
# Helper functions
########################################################################
def gaussian_kernel(x, sigma):
    kernel = np.exp(-0.5 * (x**2) / np.square(sigma)) / (sigma*np.sqrt(2*np.pi))
    return kernel

# https://stackoverflow.com/a/61151205

def confirm_input(question, default="no"):

    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '{}}'".format(default))

    while True:
        print(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            print("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def stn_pair_distances(rlist):
    crd = pd.read_csv('coordinates_receivers_recnum.csv', sep='\t', header=None)
    crd.columns = ['rec', 'lat', 'x', 'y']
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
            if (angVert):
                qang = np.arctan2(rec_x - base_x, rec_y - base_y)
            else:
                qang = np.arctan2(rec_y - base_y, rec_x - base_x)
            
            if (qang < 0.0):
                qang += 2*np.pi

            rpazim_d = qang*180/np.pi
            distances = np.append(distances, dist)
            angles = np.append(angles, rpazim_d)

    np.save('Results/%s/distances.npy' %(name), distances)
    np.save('Results/%s/angles.npy' %(name), angles)

    #return distances, angles


########################################################################
# Code execution
########################################################################
if __name__ == '__main__':

    rewrite = True
    angVert = False

    nargs=len(sys.argv)
    print("\n-----------------------------------------------------------------------------")
    rec_input=raw_input("\nEnter receiver numbers (leave empty to include all receivers): ")

    # Azimuth bin
    print("\n-----------------------------------------------------------------------------")
    azbin=raw_input("\nDesired angle bin width (default=10): ")
    azbin = int(azbin) if azbin else 10
    # Desired receiver angles

    print("\n-----------------------------------------------------------------------------")
    rec_angs=raw_input("\nEnter the desired angle (leave empty to include all angles): ")
    if len(rec_angs.split())>0:
        azms=np.array(rec_angs.split(), dtype=np.int16)
        if np.sum(np.array([int(azm)%5 for azm in azms])) != 0:
            sys.exit("One or more invalid angle(s). Angles should multiple of %d" %(azbin//2))
    else:
        azms = np.arange(azbin//2, 360, azbin)

    gaussw=None
    if confirm_input("\nAdd gaussian weights? ", 'yes'):
        gaussw=input('Enter the value of sigma: ')
        gaussw=gaussw if gaussw else 10

    linw=False
    if confirm_input("\nAdd linear weights? ", 'yes'):
        linw=True


    for fn in range(1,nargs):
        if not path.isfile(sys.argv[fn]):
            sys.exit("File not found: %s" %(sys.argv[fn]))
        print "\nReading ", sys.argv[fn]

        name=sys.argv[fn][0:-4]
        name=name if not gaussw else name+str(gaussw)
        if not path.exists('Results/%s' %(name)):
            makedirs('Results/%s' %(name))

        loaded = np.load(sys.argv[fn])
        reclist=loaded['reclist']
        si=loaded['si'][0]
        cclags=loaded['cclags']
        cookie=loaded['cookie']
        cclags=np.array(map(lambda x: float("%.2f" %x), cclags))
        print sys.argv[fn], "reading done!"

        if len(rec_input.split())>0:
            rlist=np.array(rec_input.split(), dtype=np.int16)
        else:
            rlist=reclist
        
        del loaded
        if (rewrite):
            stn_pair_distances(reclist)

        print "\nReading angles and distances.."
        angles=np.load('Results/%s/angles.npy' %(name))
        distances=np.load('Results/%s/distances.npy' %(name))
        print 'Done!\n'

        #########################################################################################
        # Result
        #########################################################################################
        print("\n##############################################################################")
        print('Results for file %s.npz' %(name))
        print("##############################################################################\n")

        rso = do_one_brec(cookie,rlist,azms,azbin)
        rso.save_files()

    def dict_to_csv(flname):
        fl=np.load(flname, allow_pickle=True).item()
        df=pd.DataFrame.from_dict(fl, orient='index')
        df.to_csv('Results/percentile.csv')

    dict_to_csv('percentile.npy')