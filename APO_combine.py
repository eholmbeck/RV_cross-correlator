"""
Author:		E.M. Holmbeck
Last edit:	22 Feb '23
Output:		Writes a file (currently "APO_stack.cl") in the called directory. Open IRAF and run `cl < APO_stack.cl` to actually combine the spectra.

How to use:
Call `python APO_combine.py`. When prompted, enter the path to the directory where you have exposures to combine.
Note that the spectra need to be RV-shifted first (or otherwise at similar Doppler shifts).
Adjust the "t_tol" variable if you only want to combine spectra within some amount of observed time (e.g., if exposures are taken within two days).

Warning:	This is all written for APO right now. Adjust bands, labels, etc., accordingly for other data.

"""

from glob import glob
from astropy.io import fits
from astropy.time import Time
from numpy import inf
from spectrum import FITS
import os.path as path


# Tolerance in days for combining spectra. If the time difference is greater than t_tol, we won't combine them.
t_tol = inf

# Enter path to directory for which you want to combine files
dir = input('Directory: ')
#dir = path.abspath(dir)
dir = dir.replace('\ ', ' ').strip()

# Find all files in root directory and sub-directories that have already been shifted
# NOTE: I haven't really implemented recursive searching yet...
files = glob(path.join(dir,'**/*_rv.fits'), recursive=True)

# Function for displaying all exposures of a given object
# Choose beam/aperture to show. It's on H-beta now
def show_spectra(files, beam=117):
	import matplotlib.pyplot as plt
	fig,ax = plt.subplots(1,1, figsize=(9,6))
	
	lines,handles = [[],[]]
	for i,file in enumerate(files):
		# Skip displaying the combined spectrum, if it exists
		if 'comb' in file: file = file.replace('_rv.fits', '.fits')
		spec = FITS(file)
		try:
			dopcor = spec.header['DOPCOR']
			try: dopcor = float(dopcor)
			except: dopcor = float(dopcor.split()[0])
			#print(dopcor, float(spec.header['VHELIO']))
		except: pass
	
		tax = ax.twinx()
		tax.plot(spec.wavelength[beam], spec.data[beam], label=file.split('/')[-1], lw=0.7, color='C%s' %i),

		l,h = tax.get_legend_handles_labels()
		lines.append(l[0])
		handles.append(h[0])

	ax.legend(lines,handles)
	plt.xlim(4850, 4870)
	plt.show(block=False)
	keep = input('Is this okay? (Use "n" or "x" if not): ')
	plt.close()
	
	if keep == '':
		return True
		
	if keep in 'nNxX':
		return False
	
	return True


# Getting ducks in a row for combining.
names = {}
for f in files:
	name = f.split('/')[-1].split('_APO')[0]
	
	hdul = fits.open(f)
	header = hdul[0].header
	try:
		ncomb = header['NCOMBINE']
	except KeyError:
		ncomb = 1
	
	tobs = Time(header['DATE-OBS'], format='isot', scale='tai')
	bands = len(hdul[0].shape)
	#names = find_night(name, names, tobs)
	
	if name in names:
		TDelta = []

		# I don't think I *need* to iterate over all of them... If the first two are within 1 day,
		# and the third one is within a day of the first, it should be fine???
		#for file,time in names[name]:
		file, time, _, _ = names[name][0]
		TD = tobs - time
		if abs(TD.jd) < t_tol:
			print('Time difference is less than %s days. I will combine these two spectra.' %t_tol)
			names[name].append([f,tobs,ncomb,bands])
		
		else:
			night = 2
			name_nextnight = name+'_n%s' %night
			found = False
			while name_nextnight in names:
				file, time, _, _ = names[name_nextnight][0]
				TD = tobs - time
				print(TD.jd)
				if abs(TD.jd) < t_tol:
					print('Time difference is less than %s days. I will combine these two spectra.' %t_tol)
					names[name_nextnight].append([f,tobs,ncomb,bands])
					found = True
					break
				
				night += 1
				name_nextnight = name+'_n%s' %night
			
			if not found:
				names[name_nextnight] = [[f,tobs,ncomb,bands]]
	
	else: names[name] = [[f,tobs,ncomb,bands]]



# Write the iraf script.
iraf_script = open(path.join(dir, 'APO_stack.cl'), 'w')

for ni,name in enumerate(names):
	# If there is only one exposure, no need to combine
	if len(names[name]) < 2: continue
	
	# If the combined spectrum exists, continue
	find_comb = glob(path.join(dir, '%s_APO35_comb.fits' %name))
	if len(find_comb)>0:
		print('Found ', find_comb, '. Skipping this one.')
		continue
	
	# Check if the data are 3-D (multiband, which includes error spectra) or not
	bands = 3

	# Create a file list of each band
	files = []
	file_list1 = open(path.join(dir,'files%s_band1' %ni), 'w')
	file_list2 = open(path.join(dir,'files%s_band2' %ni), 'w')
	file_list3 = open(path.join(dir,'files%s_band3' %ni), 'w')
	file_list4 = open(path.join(dir,'files%s_band4' %ni), 'w')
	# Error spectrum needs to be squared, then added, so we need an intermediate squared step
	file_listsqr = open(path.join(dir,'files%s_band4_sqr' %ni), 'w')
	
	for i,(file,_,weight,band) in enumerate(names[name]):
		new_fname = file.split('/')[-1]
		#weights.write('%s\n' %weight)
		files.append(file)
		for fi,file_list in enumerate([file_list1,file_list2,file_list3,file_list4]):
			# Duplicate for as many NCOMBINE went into that "exposure"
			file_list.write(('%s[*,*,%s]\n' %(new_fname,fi+1))*weight)

		file_listsqr.write(('%s_squared.fits\n' %new_fname.replace('.fits',''))*weight)
		
		if band<bands: bands=band
	
	
	for file_list in [file_list1,file_list2,file_list3,file_list4]:
		file_list.close()
	
	# Reassign variances to point to the file name
	file_list1 = path.join(dir,'files%s_band1' %ni)
	file_list2 = path.join(dir,'files%s_band2' %ni)
	file_list3 = path.join(dir,'files%s_band3' %ni)
	file_list4 = path.join(dir,'files%s_band4' %ni)
	file_listsqr = path.join(dir,'files%s_band4_sqr' %ni)
	
	# Show exposures and skip bad ones
	if not show_spectra(files):
		print('%s flagged as possible binary.' %name)
		binaries = open(path.join(dir,'binaries_and_baddies.txt'), 'a')
		binaries.write(name+'\n')
		binaries.close()
		iraf_script.write('delete %s\n' %(','.join([file_list1,file_list2,file_list3,file_list4,file_listsqr])))
		continue
	
	# Adding the INSTR back to the star_ID
	name = name.replace('.fits', '')
	if '_n' in name:
		name = name.replace('_n', '_APO35_n')
	else:
		name += '_APO35'
	
	# Options for scombine routine.
	opts = 'first=yes reject=avsigclip combine=sum weight=none'

	# If all exposures have the error spectrum, combine all bands
	if bands==3:
		# Intermediate file names
		specname = name+'_speccomb'
		sigadd = name+'_sigadd'
		signame = name+'_sigcomb'
		band2_name = name+'_band2'
		band3_name = name+'_band3'
		
		# Combine each band; for the sky and background, I'm just adding them...
		iraf_script.write("scombine @files%s_band1 output=%s %s\n" %(ni,specname,opts))
		iraf_script.write("scombine @files%s_band2 output=%s %s\n" %(ni,band2_name,opts))
		iraf_script.write("scombine @files%s_band3 output=%s %s\n" %(ni,band3_name,opts))

		# For the error spectrum. Square the individual spectra, add them, then squareroot
		iraf_script.write('sarith @files%s_band4 ^ 2 output=@files%s_band4_sqr\n' %(ni, ni))
		iraf_script.write("scombine @files%s_band4_sqr output=%s %s\n" %(ni,sigadd,opts))
		iraf_script.write('sarith %s sqrt output=%s\n' %(sigadd,signame))
		
		# imstack all bands
		final_name = name+'_comb'
		iraf_script.write('imstack %s,%s,%s,%s %s\n' %(specname,name+'_band2',name+'_band3',signame,final_name))
		
		# Delete all intermediate files and lists
		iraf_script.write('delete @%s\n' %file_listsqr)
		iraf_script.write('delete %s.fits\n' %band2_name)
		iraf_script.write('delete %s.fits\n' %band3_name)
		iraf_script.write('delete %s.fits\n' %sigadd)
		iraf_script.write('delete %s.fits\n' %signame)
		iraf_script.write('delete %s.fits\n' %specname)
		iraf_script.write('hedit %s DOPCOR,VHELIO delete=yes\n' %final_name)
	
	# If the data are only 2-D (no error spectrum), just combine the data
	elif bands==2:
		specname = name+'_comb'
		iraf_script.write("scombine @files%s_band1 output=%s %s\n" %(ni,specname,opts))
		iraf_script.write('hedit %s DOPCOR,VHELIO delete=yes\n' %specname)
	
	# Clean up
	for file in [file_list1,file_list2,file_list3,file_list4,file_listsqr]:
		iraf_script.write('delete %s\n' %file)
		
iraf_script.close()

	
	
