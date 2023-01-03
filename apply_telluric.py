import sys
from pypeit.core.telluric import general_spec_reader, save_coadd1d_tofits
from astropy.io import fits, ascii as asci
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from pypeit import utils
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter



sci  = sys.argv[1]
tell = sys.argv[2]

wave, wave_grid_mid, flux, ivar, mask, meta_spec, header = general_spec_reader(sci, ret_flam=False)

#Read the telluric object
TelObj = fits.open(tell)[1].data

tel_wave = TelObj['WAVE'][0]
tel_tel  = TelObj['TELLURIC'][0]
# print(tel_wave)

#interpolate the telluric spectrum to the observed wavelength
tell_interp = scipy.interpolate.interp1d(tel_wave, tel_tel, bounds_error=False, fill_value=0.0)
telluric = tell_interp(wave)

flux_corr = flux*utils.inverse(telluric)
ivar_corr = (telluric > 0.0) * ivar * telluric * telluric
mask_corr = (telluric > 0.0) * mask
sig_corr = np.sqrt(utils.inverse(ivar_corr))


out_name = sci.split('_')[0]+"_"+tell.split('_')[0]+"_telluric_corrected.txt"

#naming things
hdr = fits.open(sci)

tel = hdr[0].header['TELESCOP']
inst = hdr[0].header['INSTRUME']
date_obj = Time(hdr[0].header['MJD'], format = 'mjd')
date = date_obj.isot[0:10]
timeobs = date_obj.isot[11:]
obj = (hdr[0].header['TARGET']).strip()
sn_coord = SkyCoord(ra = hdr[0].header['RA']*u.deg, dec = hdr[0].header['DEC']*u.deg)
 
out_name2 = obj+'_'+tel+'_'+inst+'_'+date+'.flm'
out_plot = obj+'_'+tel+'_'+inst+'_'+date+'.png'

#Define tables to write out

spec = Table(np.array([wave,flux_corr, sig_corr]).T, names = ['wavelength', 'flux', 'fluxerr'])

spec.meta['comments'] = ['wavelength flux fluxerr', 'GROUPS UCSC,YSE,KITS', 'SNID %s'%(obj), 'OBS_GROUP UCSC', 
							'OBS_DATE %s %s'%(date, timeobs), 'INSTRUMENT %s'%inst,
							'RA %f'%sn_coord.ra.deg, 'DEC %f'%sn_coord.dec.deg]

asci.write(np.array([wave,flux_corr, sig_corr]).T, out_name, \
    names = ['wavelength', 'flux', 'fluxerr'], format = 'csv', overwrite = True)

asci.write(spec , out_name2, format = 'no_header', delimiter = ' ', overwrite = True)


######PLOT
YJ = np.logical_and(spec['wavelength'] > 10000, spec['wavelength'] < 13000)
H  = np.logical_and(spec['wavelength'] > 14200, spec['wavelength'] < 17500)
K  = np.logical_and(spec['wavelength'] > 20000, spec['wavelength'] < 24000)

# good_snr = np.logical_and(good, spec[0,:] , 13500)
good_snr = np.logical_or(YJ, H)
good_snr = np.logical_or(good_snr, K)

minwl = None
if minwl is not None:
    good_snr = np.logical_and(good_snr, spec['wavelength'] > minwl)

max_y = np.nanmax(median_filter(spec['flux'][good_snr], size= 5))
min_y = np.nanmin(median_filter(spec['flux'][good_snr], size= 5))
# print(fn)
fig, ax = plt.subplots(1,1,figsize = (15,8))

#smooth with median 
smooth_size = None
if smooth_size is None:
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
fig.savefig(out_plot, bbox_inches = 'tight')


# plt.step(wave, flux, ':', where = 'mid', lw = 1)
# plt.step(wave, flux_corr, where = 'mid', lw = 1)
# plt.show()