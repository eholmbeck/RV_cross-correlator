import numpy as np
from astropy.io import fits


class FITS:
	def __init__(self,file):
		self.file_name = file
		# Put all this data reassignment to a subroutine; I think the file is being overwritten.
		
		for method in ["self.__read_dispersion__()", "self.__read_neid__()"]:
			try:
				self.wavelength = eval(method)
				break
			except: continue
		
	def __read_neid__(self):
		self.hdu = fits.open(self.file_name)
		self.first_beam = 0
		self.dispersion = None

		self.header = self.hdu[0].header
		'''
		if 'SCIWAVE' not in self.hdu:
			image_header = hdu[0].header
		else:
			image_header = hdu['SCIWAVE'].header
		'''
		
		self.data = self.hdu['SCIFLUX'].data
		self.apertures = self.data.shape[0]
		wavelength = self.hdu['SCIWAVE'].data
		
		'''
		metadata = {}
		for key, value in hdu.items():
			if key in metadata:
				metadata[key] += value
			else:
				metadata[key] = value
		
		self.dispersion = {}
		wavelength = {aperture:np.zeros(metadata['NAXIS1']) for aperture in range(metadata['NAXIS2'])}
		for order in range(1,metadata['NAXIS2']+1):
			crval = metadata[f'CRVAL{order:.0f}']
			cdelt = metadata[f'CDELT{order:.0f}']
			poly_type = metadata[f'PS{order:.0f}_0']
			poly_param_order = metadata[f'PS{order:.0f}_1']
			#pixels = np.arange(*eval(metadata[f'PS{order:.0f}_2']))
			echelle_order = metadata[f'PS{order:.0f}_3']
			leg_range = eval(metadata[f'PS{order:.0f}_4'])
			#poly_params = [metadata[f'PV{order:.0f}_{param:.0f}'] for param in range(poly_param_order)]
			poly_params = []
			for param in range(poly_param_order):
				value = metadata[f'PV{order:.0f}_{param:.0f}']
				if value is None: poly_params.append(np.nan)
				else: poly_params.append(value)
			
			x = np.linspace(leg_range[0], leg_range[1], metadata['NAXIS1'])
			if poly_type=='Legendre' or 'legendre' in poly_type:
				from scipy.special import legendre
				for li in range(poly_param_order):
					wavelength[order-1] += poly_params[li]*legendre(li)(x)
			
			else: print(poly_type)
			
			self.dispersion[order-1] = self.data[order-1]
		
		self.data = self.dispersion
		'''
		aps = range(len(wavelength))
		self.first_beam = min(aps)
		
		data_copy = self.data
		self.data = {}
		for ai,a in enumerate(aps):
			self.data[a] = data_copy[ai]
		
		self.ra_key = 'QRA'
		self.dec_key = 'QDEC'
		return wavelength

	def __read_dispersion__(self):
		self.hdu = fits.open(self.file_name)
		self.header = self.hdu[0].header

		if len(self.hdu[0].shape) == 3:
			self.data = self.hdu[0].data[1]
		else:
			self.data = self.hdu[0].data

		self.apertures = self.data.shape[0]
		self.dispersion = None
		
		fits.conf.strip_header_whitespace = False
		dispersion = ['']*self.apertures

		# Every time the counter is even, it means we are between the double quote,
		# and we should keep that information. 
		# The number of double quotes in a line will be len(split('"'))-1
		# Therefore, at the end of each line, reduce the counter by 1

		counter = 0
		for key in sorted(self.header.keys()):
			if "WAT2_" in key:
				line = self.header[key].split('"')
				for l in line:
					counter += 1
					if len(l)==0: continue
					if counter%2==0:
						dispersion[int(counter/2)-1]+=l
						
				counter -= 1
		
		self.dispersion = dispersion
		
		wavelength = {}
		beam_previous = 1
		rescale = np.float(np.shape(self.data)[1])
		for ap_num in range(self.apertures):
			disp_data = list(map(float,dispersion[ap_num].split()))
			aperture, beam, dispersion_type, dispersion_start, \
				mean_dispersion_delta, num_pixels = disp_data[0:6]

			self.dopshift = disp_data[6]
			self.dopcor = np.sqrt((1.0+self.dopshift)/(1.0-self.dopshift))
			#dispersion_start *= (1.0-dopshift)
			
			#num_pixels = int(num_pixels)
			# Most of the time, this is just mean_dispersion_delta. Mostly.
			disp_delta = mean_dispersion_delta * num_pixels / rescale
			
			if dispersion_type == 0:
				if beam == beam_previous:
					beam_previous = beam
					#wavelength[aperture] = dispersion_start + np.arange(num_pixels) * mean_dispersion_delta
					wavelength[aperture] = dispersion_start + np.arange(int(rescale)) * disp_delta
					wavelength[aperture] /= self.dopcor
				else:
					#wavelength[beam] = dispersion_start + np.arange(num_pixels) * mean_dispersion_delta
					wavelength[beam] = dispersion_start + np.arange(int(rescale)) * disp_delta
					wavelength[beam] /= self.dopcor

				
			elif dispersion_type == 2:
				num_pixels = int(num_pixels)
				coefficients = disp_data[11:]
				# Chebyshev
				# https://iraf.net/irafhelp.php?val=aptrace&help=Help+Page
				if coefficients[0] == 1:
					order,xmin,xmax = coefficients[1:4]
					xmin -= 1 #since python indexes at 0
					xmax -= 1
					x = np.linspace(xmin, xmax, num_pixels)
					n = (2.0*x - (xmax + xmin)) / (xmax - xmin)
					z = []
					z.append([1.0]*num_pixels)
					z.append(n)
					
					for i in range(2, int(order)):
						z.append(2.0*n*z[i-1] - z[i-2])
					
					z = np.array(z)
					c = np.array(coefficients[4:])
					wavelength[beam] = np.dot(c,z)
					
					wavelength[beam] /= self.dopcor
	
					#6570.60 # What it's reading
					#6570.97 # What it should be
				else:
					print("Function type not yet supported. Poke Erika!")
					exit()
			
			else:
				print("Log-linear binned not yet supported. Poke Erika to implement this!")
				exit()
			
		aps = list(map(int, wavelength.keys()))
		self.first_beam = min(aps)
		
		data_copy = self.data
		self.data = {}
		#for a in range(self.apertures):
		#	self.data[a+self.first_beam] = data_copy[a]
		# Using beam instead of aperture
		for ai,a in enumerate(aps):
			#try: self.data[a] = data_copy[a-self.first_beam-new_start]
			self.data[a] = data_copy[ai]
		
		self.ra_key = 'RA'
		self.dec_key = 'DEC'
		return wavelength
		
	'''	
	def wavelength(self, ap_num):
		disp_data = map(float,self.dispersion[ap_num].split())
		
		aperture, beam, dispersion_type, dispersion_start, \
			mean_dispersion_delta, num_pixels = disp_data[0:6]
		
		if aperture not in self.apertures:
			# Assume linear dispersion for now...
			self.apertures[aperture] = dispersion_start + np.arange(num_pixels) * mean_dispersion_delta
	
		return self.apertures[aperture]
	'''

def show_spectra(files):
	import matplotlib.pyplot as plt
	fig,ax = plt.subplots(1,1, figsize=(9,6))

	beam = 117
	lines,handles = [[],[]]
	for i,file in enumerate(files):
		#if 'comb' in file: file = file.replace('_rv.fits', '.fits')
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
	plt.show()
	
if False:
	import glob
	files = glob.glob('/Volumes/GoogleDrive/Shared drives/RPA/data/Reduced_rvcorrected/202008_APO/*_rv.fits')

	for i in range(int(len(files)/5)):
		if (len(files)-(5*i)) < 5:
			show_spectra(files[5*i:])
		else:
			show_spectra(files[5*i:5*(i+1)])
	
	#show_spectra(files[5*(i+1):])
	show_spectra(files[-5:])
