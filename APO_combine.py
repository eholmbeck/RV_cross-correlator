from glob import glob
from astropy.io import fits
from astropy.time import Time
from numpy import inf
from spectrum import FITS
import os.path as path

t_tol = inf

dir = input('Directory: ')
dir = dir.replace('\ ', ' ').strip()
print(dir)
files = glob(dir+'/*_rv.fits')

'''
def find_night(name, names, tobs, night=2):
	# Reinitialize; otherwise, we get "n2_n2"
	#name = f.split('/')[-1].split('_APO')[0]
	if name in names:
		TDelta = []
		# I don't think I *need* to iterate over all of them... If the first two are within 1 day,
		# and the third one is within a day of the first, it should be fine???
		#for file,time in names[name]:
		file, time = names[name][0]
		TD = tobs - time
		#import pdb
		#pdb.set_trace()
		if TD.jd < 1.5:
			print('Time difference is less than 1.5 days. I will combine these spectra.')
			names[name].append([f,tobs])
			return names
		
		# Otherwise, look for a new dictionary key with "_n2" appended.
		# If it exists, compare to those spectra.
		else:
			while name+'_n%s' %night in names:
				print(names.keys())
				names = find_night(name+'_n%s' %night, names, tobs, night)
				night += 1
			
			names[name+'_n%s' %night] = [[f, tobs]]
			return names
			#names = find_night(name+'_n%s' %night, names, tobs)	
			
	else: names[name] = [[f, tobs]]
	return names
'''

def show_spectra(files):
	import matplotlib.pyplot as plt
	fig,ax = plt.subplots(1,1, figsize=(9,6))

	beam = 117
	lines,handles = [[],[]]
	for i,file in enumerate(files):
		if 'comb' in file: file = file.replace('_rv.fits', '.fits')
		spec = FITS(file)
		try:
			dopcor = spec.header['DOPCOR']
			try: dopcor = float(dopcor)
			except: dopcor = float(dopcor.split()[0])
			#print(dopcor, float(spec.header['VHELIO']))
		except: pass
	
		tax = ax.twinx()
		#if i==len(files)-1:
		#	tax.plot(spec.wavelength[beam], spec.data[beam], label=file.split('/')[-1], lw=1.5, color='gray'),
		#else:
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


names = {}
for f in files:
	name = f.split('/')[-1].split('_APO')[0]
	
	hdul = fits.open(f)
	header = hdul[0].header
	try:
		ncomb = header['NCOMBINE']
	except KeyError:
		#with fits.open(f, mode='update') as hdul: hdul[0].header.append(('NCOMBINE', 1))
		ncomb = 1
	
	tobs = Time(header['DATE-OBS'], format='isot', scale='tai')
	bands = len(hdul[0].shape)
	#names = find_night(name, names, tobs)
	
	if name in names:
		TDelta = []
		#import pdb
		#pdb.set_trace()

		# I don't think I *need* to iterate over all of them... If the first two are within 1 day,
		# and the third one is within a day of the first, it should be fine???
		#for file,time in names[name]:
		file, time, _, _ = names[name][0]
		TD = tobs - time
		print(TD.jd)
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

# TO-DO:


iraf_script = open('APO_stack.cl', 'w')

for ni,name in enumerate(names):
	if len(names[name]) < 2: continue
	find_comb = glob(path.join(dir, '%s_APO35_comb.fits' %name))
	if len(find_comb)>0:
		print('Found ', find_comb, '. Skipping this one.')
		continue

	weights = open(dir+'weights%s' %ni, 'w')
	
	bands = 3

	files = []
	file_list1 = open(dir+'files%s_band1' %ni, 'w')
	file_list2 = open(dir+'files%s_band2' %ni, 'w')
	file_list3 = open(dir+'files%s_band3' %ni, 'w')
	file_list4 = open(dir+'files%s_band4' %ni, 'w')
	file_listsqr = open(dir+'files%s_band4_sqr' %ni, 'w')
	
	for i,(file,_,weight,band) in enumerate(names[name]):
		#new_fname = file.replace('/Volumes/GoogleDrive/Shared drives', '/Users/eholmbeck/Desktop/Drive')
		new_fname = file.split('/')[-1]
		weights.write('%s\n' %weight)
		files.append(file)
		for fi,file_list in enumerate([file_list1,file_list2,file_list3,file_list4]):
			file_list.write('%s[*,*,%s]\n' %(new_fname,fi+1))

		file_listsqr.write('%s_squared\n' %new_fname.replace('.fits',''))
		
		if band<bands: bands=band
	
		
	for file_list in [file_list1,file_list2,file_list3,file_list4]:
		file_list.close()
	
	weights.close()
	weights = dir+'weights%s' %ni
	file_list1 = dir+'files%s_band1' %ni
	file_list2 = dir+'files%s_band2' %ni
	file_list3 = dir+'files%s_band3' %ni
	file_list4 = dir+'files%s_band4' %ni
	file_listsqr = dir+'files%s_band4_sqr' %ni
	
	if not show_spectra(files):
		print('%s possibly flagged as binary.' %name)
		binaries = open(dir+'/binaries_and_baddies.txt', 'a')
		binaries.write(name+'\n')
		binaries.close()
		iraf_script.write('delete %s\n' %(','.join([weights,file_list1,file_list2,file_list3,file_list4,file_listsqr])))
		continue
	
	name = name.replace('.fits', '')
	if '_n' in name:
		name = name.replace('_n', '_APO35_n')
	else:
		name += '_APO35'

	if bands==3:
		specname = name+'_speccomb'
		string = "scombine @files%s_band1 output=%s first=yes reject=avsigclip combine=average weight=@%s\n" %(ni,specname,weights)
		if len(string)<=132:
			print('String too long!')
			iraf_script.write(string)
		else:
			binaries = open(dir+'/binaries_and_baddies.txt', 'a')
			binaries.write(name+' input string too long\n')
			binaries.close()
			continue
	
		sigadd = name+'_sigadd'
		signame = name+'_sigcomb'
		band2_name = name+'_band2'
		band3_name = name+'_band3'

		# Square the individual spectra, add them, then squareroot
		iraf_script.write("scombine @files%s_band2 output=%s first=yes reject=avsigclip combine=average weight=@%s\n" %(ni,band2_name,weights))
		iraf_script.write("scombine @files%s_band3 output=%s first=yes reject=avsigclip combine=average weight=@%s\n" %(ni,band3_name,weights))
		iraf_script.write('sarith @files%s_band4 ^ 2 output=@files%s_band4_sqr\n' %(ni, ni))
		iraf_script.write("scombine @files%s_band4_sqr output=%s first=yes reject=avsigclip combine=average weight=@%s\n" %(ni,sigadd,weights))
		iraf_script.write('sarith %s sqrt output=%s\n' %(sigadd,signame))
	
		final_name = name+'_comb'
		iraf_script.write('imstack %s,%s,%s,%s %s\n' %(specname,name+'_band2',name+'_band3',signame,final_name))
	
		iraf_script.write('delete @%s\n' %file_listsqr)
		iraf_script.write('delete %s.fits\n' %band2_name)
		iraf_script.write('delete %s.fits\n' %band3_name)
		iraf_script.write('delete %s.fits\n' %sigadd)
		iraf_script.write('delete %s.fits\n' %signame)
		iraf_script.write('delete %s.fits\n' %specname)
		iraf_script.write('hedit %s DOPCOR,VHELIO delete=yes\n' %final_name)

	elif bands==2:
		specname = name+'_comb'
		iraf_script.write("scombine %s output=%s first=yes reject=avsigclip combine=average weight=@%s\n" %(list_band1.replace('[*,*,1]',''),specname,weights))
		iraf_script.write('hedit %s DOPCOR,VHELIO delete=yes\n' %specname)
	
	iraf_script.write('delete %s\n' %(','.join([weights,file_list1,file_list2,file_list3,file_list4,file_listsqr])))
		
		
iraf_script.close()

	
	
