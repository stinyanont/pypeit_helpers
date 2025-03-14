#Convert pypeit output in FITS into a plot and ascii file
import astropy.io.fits as fits
import astropy.io.ascii as asci
import numpy as np 
import os, sys, glob, argparse
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import pdb

if __name__ == '__main__':
    #read files in


    parser = argparse.ArgumentParser()

    parser.add_argument('filename')     # files to be converted
    parser.add_argument('--minwl')      # minimum wavelength
    parser.add_argument('--obj_name')
    parser.add_argument('--smooth_size')

    args = parser.parse_args()

    files = glob.glob(args.filename)
    files.sort()

    if args.minwl is not None:
        minwl = float(args.minwl)
    else:
        minwl = None
    obj_name = args.obj_name
    if args.smooth_size is not None:
        smooth_size = int(args.smooth_size)
    else:
        smooth_size = None

    print(minwl, obj_name, smooth_size)

    # if len(sys.argv) > 2:
    #     obj_name = sys.argv[2]
    # else:
    #     obj_name = None

    # if len(sys.argv) > 3:
    #     smooth_size = int(sys.argv[3])
    # else:
    #     smooth_size = None

    # if len(sys.argv) > 4 :
    #     minwl= float(sys.argv[4])
    # else:
    #     minwl = None

    for fn in files:
        hdr = fits.open(fn)
        if obj_name is None:
            obj = (hdr[0].header['TARGET']).strip()
        else:
            obj = obj_name
        
        tel = hdr[0].header['TELESCOP']
        inst = hdr[0].header['INSTRUME']
        date_obj = Time(hdr[0].header['MJD'], format = 'mjd')
        date = date_obj.isot[0:10]
        timeobs = date_obj.isot[11:]

        spec_data = hdr[1].data
        spec = Table([spec_data['wave'], spec_data['flux'], np.sqrt(spec_data['ivar'])], names = ('wavelength', 'flux', 'fluxerr'))
        sn_coord = SkyCoord(ra = hdr[0].header['RA']*u.deg, dec = hdr[0].header['DEC']*u.deg)
        ###Table to write
        spec.meta['comments'] = ['wavelength flux fluxerr', 'GROUPS UCSC,YSE,KITS', 'SNID %s'%(obj), 'OBS_GROUP UCSC', 
									'OBS_DATE %s %s'%(date, timeobs), 'INSTRUMENT %s'%inst,
									'RA %f'%sn_coord.ra.deg, 'DEC %f'%sn_coord.dec.deg]


        out_fn = obj+'_'+tel+'_'+inst+'_'+date+'.png'
        out_txt = obj+'_'+tel+'_'+inst+'_'+date+'.flm'
        asci.write(spec, out_txt, format = 'no_header', delimiter = ' ', overwrite = True)

        ##Deal with range

        #Clip bogus points with wl = 0
        if np.min(spec['wavelength']) < 2000:
            spec = spec[spec['wavelength'] > 4000]     
        # good_snr = spec[1,:]/spec[2,:] > 5
        if np.max(spec['wavelength'] > 15000):
            YJ = np.logical_and(spec['wavelength'] > 10000, spec['wavelength'] < 13000)
            H  = np.logical_and(spec['wavelength'] > 14200, spec['wavelength'] < 17500)
            K  = np.logical_and(spec['wavelength'] > 20000, spec['wavelength'] < 24000)

            # good_snr = np.logical_and(good, spec[0,:] , 13500)
            good_snr = np.logical_or(YJ, H)
            good_snr = np.logical_or(good_snr, K)
        # elif np.min(spec['wavelength'] < 3000):
        #     good_snr = spec['wavelength'] > 3000
        else:
            good_snr = (spec['flux'] > 5*spec['fluxerr'])

        if minwl is not None:
            good_snr = np.logical_and(good_snr, spec['wavelength'] > minwl)
        # pdb.set_trace()
        max_y = np.nanmax(median_filter(spec['flux'][good_snr], size= 5))
        min_y = np.nanmin(median_filter(spec['flux'][good_snr], size= 5))
        print(fn, min_y, max_y)
        fig, ax = plt.subplots(1,1,figsize = (15,8))

        #smooth with median 
        if (smooth_size is None) | (smooth_size == 0):
            ax.step(spec['wavelength'], spec['flux'], lw=1)
        else:
            ax.step(spec['wavelength'], median_filter(spec['flux'], size = smooth_size))
        ax.tick_params(labelsize = 18)
        ax.set_ylabel(r'$F_\lambda$ ($\rm erg\, s^{-1} \, cm^{-2}\, \AA^{-1}$)', fontsize = 18)
        ax.set_xlabel(r'Wavelength ($\rm \AA$)', fontsize = 18)
        ax.set_title(obj, fontsize = 18)
        # if min_y <= 0:
        # 	ax.set_ylim([1.2*min_y, 1.1*max_y])
        # else:
        # 	ax.set_ylim([0, 1.1*max_y])
        exp = int(np.log10(max_y))
        ax.set_ylim([-0.1*(max_y), 1.1*max_y])
        # ax.set_ylim([-0.1*(max_y), 1.1*max_y])

        # if minwl is not None:
            # ax.set_xlim
        fig.savefig(out_fn, bbox_inches = 'tight')


