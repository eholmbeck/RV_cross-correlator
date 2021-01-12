import numpy as np
from astropy.io import fits


class FITS:
	def __init__(self,file):
		self.file_name = file
		self.hdu = fits.open(self.file_name)
		self.header = self.hdu[0].header
		if len(self.hdu[0].shape) == 3:
			self.data = self.hdu[0].data[1]
		else:
			self.data = self.hdu[0].data
		
		self.apertures = self.data.shape[0]
		self.first_beam = 0
		self.dispersion = None
		self.wavelength = self.__read_dispersion__()
		
	def __read_dispersion__(self):
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
					if counter%2==0:
						dispersion[int(counter/2)-1]+=l
		
				counter -= 1
		
		self.dispersion = dispersion
		
		wavelength = {}
		beam_previous = 1
		for ap_num in range(self.apertures):
			disp_data = list(map(float,dispersion[ap_num].split()))
			aperture, beam, dispersion_type, dispersion_start, \
				mean_dispersion_delta, num_pixels = disp_data[0:6]
			
			if dispersion_type == 0:
				if beam == beam_previous:
					beam_previous = beam
					wavelength[aperture] = dispersion_start + np.arange(num_pixels) * mean_dispersion_delta
				else:
					wavelength[beam] = dispersion_start + np.arange(num_pixels) * mean_dispersion_delta

			else:
				coefficients = disp_data[6:]
				# Chebyshev
				# https://iraf.net/irafhelp.php?val=aptrace&help=Help+Page
				if coefficients[5] == 1:
					order,xmin,xmax = coefficients[6:9]
					x = np.linspace(xmin, xmax, num_pixels)
					n = (2.0*x - (xmax + xmin)) / (xmax - xmin)
					z = [0.0]*(2+int(order))
					z[0] = 1.0
					z[1] = n
					for i in range(2, int(order)+2, 2):
						z[i] = 2.0*n*z[i-1] - z[i-2]
					z = np.array(z[0:int(order)])
					c = np.array(coefficients[9:])
					wavelength[beam] = np.sum(c*z)
				else:
					print "Function type not yet supported. Poke Erika!"
					exit()
			
		aps = list(map(int, wavelength.keys()))
		self.first_beam = min(aps)
		data_copy = self.data
		self.data = {}
		#for a in range(self.apertures):
		#	self.data[a+self.first_beam] = data_copy[a]
		# Using beam instead of aperture
		for a in aps:
			self.data[a] = data_copy[a-self.first_beam]

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