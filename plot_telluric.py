import matplotlib.pyplot as plt 
import numpy as np
import astropy.io.fits as fits
import scipy
import glob

#This script is to recreate the telluric diagnostic plots. 

all_tel = glob.glob('HIP*tellcorr.fits')

for file in all_tel:
	std_telcor = fits.open(file)

	fig, ax = plt.subplots(3,1,figsize = (12,12))

	good = std_telcor[1].data['wave']>0

	#Full plot
	ax[0].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['flux'][good], lw = 1)
	ax[0].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['flux'][good] * std_telcor[1].data['telluric'][good], lw = 1)
	ax[0].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['obj_model'][good], color = 'k', lw = 1)

	ax[0].set_ylim([np.min(std_telcor[1].data['obj_model'][good]), np.max(std_telcor[1].data['obj_model'][good])])

	#K band CO2
	ax[1].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['flux'][good], lw = 1)
	ax[1].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['flux'][good] * std_telcor[1].data['telluric'][good], lw = 1)
	ax[1].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['obj_model'][good], color = 'k', lw = 1)

	# y_min = std_telcor[1].data['flux'][good][np.argmin(np.abs(std_telcor[1].data['wave'][good] - 21500))]
	wl = (std_telcor[1].data['wave'][good] > 19000) & (std_telcor[1].data['wave'][good] < 21500)
	y_min = np.min(std_telcor[1].data['flux'][good][wl] * std_telcor[1].data['telluric'][good][wl]) * 0.9
	y_max = std_telcor[1].data['obj_model'][good][np.argmin(np.abs(std_telcor[1].data['wave'][good] - 19000))]*1.1


	ax[1].set_xlim([19000, 21500])
	ax[1].set_ylim([y_min, y_max])

	#J band
	ax[2].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['flux'][good], lw = 1)
	ax[2].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['flux'][good] * std_telcor[1].data['telluric'][good], lw = 1)
	ax[2].plot(std_telcor[1].data['wave'][good], std_telcor[1].data['obj_model'][good], color = 'k', lw = 1)

	# y_min = std_telcor[1].data['flux'][good][np.argmin(np.abs(std_telcor[1].data['wave'][good] - 21500))]
	wl = (std_telcor[1].data['wave'][good] > 10600) & (std_telcor[1].data['wave'][good] < 13000)
	y_min = np.min(std_telcor[1].data['flux'][good] * std_telcor[1].data['telluric'][good]) * 0.9
	y_max = std_telcor[1].data['obj_model'][good][np.argmin(np.abs(std_telcor[1].data['wave'][good] - 10600))]*1.1


	ax[2].set_xlim([10600, 13000])
	ax[2].set_ylim([y_min, y_max])



	ax[0].set_xlabel(r'Wavelength ($\rm \AA$)')
	ax[1].set_xlabel(r'Wavelength ($\rm \AA$)')
	ax[2].set_xlabel(r'Wavelength ($\rm \AA$)')

	ax[0].set_ylabel(r"$F_\lambda$")
	ax[1].set_ylabel(r"$F_\lambda$")
	ax[2].set_ylabel(r"$F_\lambda$")

	fig.suptitle(std_telcor[0].header['TARGET'])

	fig.tight_layout()

	fig.savefig('%s_telluric_QA.pdf'%std_telcor[0].header['TARGET'], bbox_inches = 'tight')