import astropy.io.fits as fits
import glob

all_fn = glob.glob('s*fits') #All NIRES spectroscopy FITS

for fn in all_fn:
	hdr = fits.open(fn)
	if 'HIP' in hdr[0].header['OBJECT']:
		print(fn+': reset obstype from "%s"'%hdr[0].header['OBSTYPE']+' to "telluric"')
		hdr[0].header['OBSTYPE'] = 'telluric'
		hdr.writeto(fn, overwrite = True)