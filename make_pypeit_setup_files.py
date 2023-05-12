import astropy.io.ascii as asci
import numpy as np
import sys, glob, os
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy import table
from astroquery.simbad import Simbad
from astropy.io import fits
"""
This script takes the .pypeit file generated by pypeit_setup and create all 
pypeit setup files necessary to run sensitivity function, coadditions, and telluric correction. 
"""

files = sys.argv[1]

#Find where data start
f = open(files, 'r')
header_line = None
valid_line = 0
for line in f:
	if len(line.strip()) > 0:
		# print(line.strip())
		if line.strip()[0] != '#':
			valid_line += 1
			if '|' in line and header_line is None:
				header_line = valid_line
	# else:
	# 	valid_line += 1
# print(header_line)

giant_table = asci.read(files, format = 'fixed_width', header_start=header_line-1, data_start=header_line)
#remove "data end"
giant_table = giant_table[:-1]
# print(giant_table)

#STANDARD STARS
std_stars = list(set(giant_table['target'][giant_table['frametype'] == 'standard']))


# for std in std_stars:
Simbad.add_votable_fields("flux(V)")
Simbad.add_votable_fields("sptype")
add_mag = Simbad.query_objects([x.split('-')[0] for x in std_stars])

print("Here are all the standard stars used this night: ")

for ind, std in enumerate(std_stars):
	print('%s \t %.2f \t %s'%(std, add_mag['FLUX_V'][ind],  add_mag['SP_TYPE'][ind]))
	fn = '%s_sen_input.sen'%std 
	f = open(fn, 'w')
	# f.write("[sensfunc]\n\talgorithm = IR\n\tstar_mag = %.2f\n\tstar_type = A0"%Vmag[ind])
	f.write("[sensfunc]\n\talgorithm = IR\n\tstar_mag = %.2f\n\tstar_type = %s"%(add_mag['FLUX_V'][ind],  add_mag['SP_TYPE'][ind]))
	f.close()

	tfn = '%s_telluric_input.tel'%std 
	tf = open(tfn, 'w')
	# f.write("[sensfunc]\n\talgorithm = IR\n\tstar_mag = %.2f\n\tstar_type = A0"%Vmag[ind])
	tf.write("[telluric]\n\tobjmodel = star\n\tstar_mag = %.2f\n\tstar_type = %s"%(add_mag['FLUX_V'][ind],  add_mag['SP_TYPE'][ind]))
	tf.close()

############################################GENERATE FLUX CALIBRATION AND COADD FILES.########################################################## 
#Matching up science to standard
unique_science = table.unique(giant_table[(giant_table['frametype'] == 'arc,science,tilt')|(giant_table['frametype'] =='science')], 'target')
sci_coord = SkyCoord(ra = unique_science['ra']*u.deg, dec = unique_science['dec']*u.deg)
unique_standard = table.unique(giant_table[giant_table['frametype'] == 'standard'], 'target')
std_coord = SkyCoord(ra = unique_standard['ra']*u.deg, dec = unique_standard['dec']*u.deg)

print(unique_science['target','ra','dec','mjd'])
print(unique_standard['target','ra','dec','mjd'])

print("Check if the following science - telluric association is correct.")

res = []
for ind, i in enumerate(unique_science):
	sep = sci_coord[ind].separation(std_coord)
	am_diff = unique_science[ind]['airmass'] - unique_standard['airmass']
	time_sep = (Time(unique_standard['mjd'], format = 'mjd') - Time(i['mjd'], format = 'mjd')).to(u.min)
	print('for %s'%i['target'])
	for j in range(len(sep)):
		print(unique_standard[j]['target'], sep[j].to(u.deg), time_sep[j], am_diff[j])
	sep[np.abs(time_sep) > 55*u.min] = np.nan
	am_diff[np.abs(time_sep) > 55*u.min] = np.nan
	print("min separation: "+unique_standard[np.nanargmin(np.abs(sep))]['target'])
	print("min time: "+unique_standard[np.nanargmin(np.abs(time_sep))]['target'])
	print("min airmass diff: "+unique_standard[np.nanargmin(np.abs(am_diff))]['target'])
	res += [np.nanargmin(np.abs(sep))]


for ind, i in enumerate(unique_science):
    print(
        unique_science[ind]["target"],
        # sci_coord[ind].to_string('hmsdms'),
        unique_standard[res[ind]]["target"],
        # std_coord[res[ind]].to_string('hmsdms'),
    )

# unique_science = table.unique(giant_table[giant_table['frametype'] == 'arc,science,tilt'], 'target')
# unique_standard = table.unique(giant_table[giant_table['frametype'] == 'standard'], 'target')
# res = []
# for i in unique_science:
# 	std_idx = np.where(unique_standard['calib']==i['calib'])[0][0]
# 	print(std_idx)
# 	res+=[std_idx]

# print("Check if the following science - telluric association is correct.")
# for ind, i in enumerate(unique_science):
#     print(
#         unique_science[ind]["target"],
#         # sci_coord[ind].to_string('hmsdms'),
#         unique_standard[res[ind]]["target"],
#         # std_coord[res[ind]].to_string('hmsdms'),
#     )

#Flux Calibration file
#Now make flux calibration files
f = open('fluxcal.flux', 'w')
f.write("[fluxcalib]\n\textrap_sens = True\n")
f.write('flux read\nfilename    | sensfile\n')
for ind, i in enumerate(unique_science):
	sci_name = unique_science[ind]["target"]
	std_name = unique_standard[res[ind]]["target"]
	sci_fn = glob.glob('spec1d*%s*.fits'%sci_name)
	sci_fn.sort()
	sen_file = std_name+"_sensitivity.fits"
	for fn in sci_fn:
		f.write(fn+' | '+sen_file+'\n')
		# obj_id_full = fits.getheader(fn)['EXT0000']
		# obj_id = obj_id_full.split('-')[0]+'-'+obj_id_full.split('-')[1]
		# f.write('%s \t|\t%s\n'%(fn, obj_id))
for ind, i in enumerate(unique_standard):
	sci_name = unique_standard[ind]["target"]
	std_name = unique_standard[ind]["target"]	
	sci_fn = glob.glob('spec1d*%s*.fits'%sci_name)
	sci_fn.sort()
	sen_file = std_name+"_sensitivity.fits"
	for fn in sci_fn:
		f.write(fn+' | '+sen_file+'\n')	
f.write('flux end')
f.close()

#Now make coadd files
for ind, i in enumerate(unique_science):
	sci_name = unique_science[ind]["target"]
	std_name = unique_standard[res[ind]]["target"]
	f = open('%s_%s.coadd'%(sci_name, std_name), 'w')
	f.write("[coadd1d]\n\tcoaddfile = '%s_corrected_coadd.fits'\n\tsensfuncfile = '%s_sensitivity.fits'\n"%(sci_name, std_name))
	f.write("coadd1d read\npath .\nfilename | obj_id\n")
	#Now loop through all the 1d file from the given object
	sci_fn = glob.glob('spec1d*%s*.fits'%sci_name)
	sci_fn.sort()
	for fn in sci_fn:
		obj_id_full = fits.getheader(fn)['EXT0000']
		obj_id = obj_id_full.split('-')[0]+'-'+obj_id_full.split('-')[1]
		f.write('%s \t|\t%s\n'%(fn, obj_id))
	f.write('coadd1d end')
	f.close()

for ind, i in enumerate(unique_standard):
	sci_name = unique_standard[ind]["target"]
	std_name = unique_standard[ind]["target"]
	f = open('%s_%s.coadd'%(sci_name, std_name), 'w')
	f.write("[coadd1d]\n\tcoaddfile = '%s_corrected_coadd.fits'\n\tsensfuncfile = '%s_sensitivity.fits'\n"%(sci_name, std_name))
	f.write("coadd1d read\npath .\nfilename | obj_id\n")
	#Now loop through all the 1d file from the given object
	sci_fn = glob.glob('spec1d*%s*.fits'%sci_name)
	sci_fn.sort()
	for fn in sci_fn:
		obj_id_full = fits.getheader(fn)['EXT0000']
		obj_id = obj_id_full.split('-')[0]+'-'+obj_id_full.split('-')[1]
		f.write('%s \t|\t%s\n'%(fn, obj_id))
	f.write('coadd1d end')
	f.close()


########BASH SCRIPT TO RUN THESE#############
############Script to run sensitivity function.
f = open('fluxcals.bash', 'w')
f.write('#!/bin/bash\n')
all_sen_func = glob.glob('*sen_input.sen')

for fn in all_sen_func:
	std_name = fn.split('_')[0]
	std_spec = glob.glob('spec1d_*%s*.fits'%std_name)
	std_spec.sort()
	f.write('pypeit_sensfunc %s -o %s_sensitivity.fits -s %s\n'%(std_spec[0], std_name, fn))
f.close()


############Script to do coaddition. 
f = open('coadds.bash', 'w')
f.write('#!/bin/bash\n')
all_coadd_file = glob.glob('*.coadd')

for fn in all_coadd_file:
	f.write('pypeit_coadd_1dspec %s\n'%fn)

############Script to do telluric
f = open('telluric.bash', 'w')
f.write('#!/bin/bash\n')
for ind, i in enumerate(unique_science):
	f.write('pypeit_tellfit %s_corrected_coadd.fits --objmodel poly \n'%i["target"])
for ind, i in enumerate(unique_standard):
	f.write('pypeit_tellfit %s_corrected_coadd.fits -t %s_telluric_input.tel \n'%(i["target"], i["target"]))
f.close()

############Apply A0V telluric to science
f = open('telluric_from_A0V.bash', 'w')
f.write('#!/bin/bash\n')

#Now make coadd files
# path_to_script = '/Users/kaew/work/useful_scripts/pypeit_helpers/apply_telluric.py'
path_to_script = os.path.dirname(os.path.realpath(__file__))+'/apply_telluric.py'
for ind, i in enumerate(unique_science):
	sci_name = unique_science[ind]["target"]
	std_name = unique_standard[res[ind]]["target"]
	sci_file = sci_name+'_corrected_coadd.fits'
	std_model = std_name+'_corrected_coadd_tellmodel.fits'
	f.write("python %s %s %s \n"%(path_to_script, sci_file, std_model))
f.close()

#Lastly, convert to ascii and png

############Lastly, make a script to run everything!
f = open('do_all.bash', 'w')
f.write('#!/bin/bash\nset -e\n')

flux = open('fluxcals.bash', 'r')
coadd = open('coadds.bash', 'r')
tell = open('telluric.bash', 'r')
tellstd =open('telluric_from_A0V.bash', 'r')

for line in flux:
	if line[0] != '#':
		f.write(line)
f.write('pypeit_flux_calib fluxcal.flux\n')
for line in coadd:
	if line[0] != '#':
		f.write(line)
for line in tell:
	if line[0] != '#':
		f.write(line)
for line in tellstd:
	if line[0] != '#':
		f.write(line)


