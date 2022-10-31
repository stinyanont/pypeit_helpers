import sys
from pypeit.core.telluric import general_spec_reader, save_coadd1d_tofits
from astropy.io import fits
from pypeit import utils
import scipy
import numpy as np
import matplotlib.pyplot as plt

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


plt.step(wave, flux, ':', where = 'mid', lw = 1)
plt.step(wave, flux_corr, where = 'mid', lw = 1)
plt.show()