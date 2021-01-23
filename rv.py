import scipy.signal
from scipy.interpolate import interpolate
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
from scipy import stats

from astropy.convolution import convolve, Box1DKernel
import astropy.units as u
from astropy.constants import c
ckms = c.to(u.km/u.s).value

def gauss(x, sigma, a, mu, b):
	g = np.exp(-0.5*((x-mu)/sigma)**2)
	return a*g + b
	
def line(x, m, b):
	return m*x + b
	
def triangle(x, centroid, m1, m2, b1):
	return np.piecewise(x, [x<centroid],\
            			[lambda x:m1*x + b1 - m1*centroid,\
            			 lambda x:m2*x + b1 - m2*centroid])


# ==========================================================
# ==========================================================
#template has 61, science has 58

def cross_correlate(template, science, aperture, log=None):
	#ap_difference = science.apertures - template.apertures
	#template_aperture = aperture
	over_zero = np.where(science.data[aperture] > 0.0)[0]
	science.data[aperture] = science.data[aperture][over_zero]
	
	try:
		template.wavelength[aperture]
		science.wavelength[aperture]
	except:
		return None
		
	science.wavelength[aperture] = science.wavelength[aperture][over_zero]

	# Find over-lapping regions
	# First, flip the data if necessary
	descending = True
	'''
	# this works
	# If the science is descending and the template is ascending
	if science.wavelength[aperture][1] < science.wavelength[aperture][0]:
		descending = True
		if template.wavelength[aperture][1] > template.wavelength[aperture][0]:
			rev_arr = template.wavelength[aperture][::-1]
			template.wavelength[aperture] = rev_arr
			rev_arr = template.data[aperture][::-1]
			template.data[aperture] = rev_arr
	
	# If the science is ascending and the template is descending, reverse the template
	elif template.wavelength[aperture][1] < template.wavelength[aperture][0]:
		rev_arr = template.wavelength[aperture][::-1]
		template.wavelength[aperture] = rev_arr
		rev_arr = template.data[aperture][::-1]
		template.data[aperture] = rev_arr

	'''
	# This doesn't??
	if template.wavelength[aperture][-1] > template.wavelength[aperture][0]:
		rev_arr = template.wavelength[aperture][::-1]
		template.wavelength[aperture] = rev_arr
		rev_arr = template.data[aperture][::-1]
		template.data[aperture] = rev_arr
		#descending = True

	if science.wavelength[aperture][-1] > science.wavelength[aperture][0]:
		rev_arr = science.wavelength[aperture][::-1]
		science.wavelength[aperture] = rev_arr
		rev_arr = science.data[aperture][::-1]
		science.data[aperture] = rev_arr
		#descending = True
	
	
	if descending:
		overlap_range = [np.where(science.wavelength[aperture]<=min([max(science.wavelength[aperture]),\
							max(template.wavelength[aperture])]))[0][0],
						 np.where(science.wavelength[aperture]>=max([min(science.wavelength[aperture]),\
							min(template.wavelength[aperture])]))[0][-1]]
	else:
		overlap_range = [np.where(science.wavelength[aperture]>=max([min(science.wavelength[aperture]),\
							min(template.wavelength[aperture])]))[0][0],
						 np.where(science.wavelength[aperture]<=min([max(science.wavelength[aperture]),\
							max(template.wavelength[aperture])]))[0][-1]]

	try:
		wavelengths = science.wavelength[aperture][overlap_range[0]:overlap_range[1]+1]
		scidata = science.data[aperture][overlap_range[0]:overlap_range[1]+1]
	except:
		wavelengths = science.wavelength[aperture][overlap_range[0]:-1]
		scidata = science.data[aperture][overlap_range[0]:-1]

	
	low_SNR = False
	signal = np.average(scidata)
	SNR = np.sqrt(signal)
	if SNR<10:
		knots = 20
	else:
		knots = 40
	
	# Take out clear zeros in data at beginning and end of aperture
	# This is mostly unnecessary
	zeros_flag = scidata[0]==0.0
	while zeros_flag:
		scidata = scidata[1:]
		wavelengths = wavelengths[1:]
		zeros_flag = scidata[0]==0.0
	
	zeros_flag = scidata[-1]==0.0
	while zeros_flag:
		scidata = scidata[:-1]
		wavelengths = wavelengths[:-1]
		zeros_flag = scidata[-1]==0.0
	
	
	if descending:
		overlap_range = [np.where(template.wavelength[aperture]>=wavelengths[0])[0][-1],
						 np.where(template.wavelength[aperture]<=wavelengths[-1])[0][0]]
	
	else:
		overlap_range = [np.where(template.wavelength[aperture]<=wavelengths[0])[0][-1],
						 np.where(template.wavelength[aperture]>=wavelengths[-1])[0][0]]

	# This will cut off one at the end since python doesn't include the last index
	# This will catch the case the the overlap_range extends to the last element
	try:
		tempdata = template.data[aperture][overlap_range[0]:overlap_range[1]+1]
	except:
		tempdata = template.data[aperture][overlap_range[0]:-1]
		
	
	# Interpolate the one with a smaller wavelength spread.
	dw_science = abs((wavelengths[-1]-wavelengths[0])/float(len(wavelengths)))
	dw_template = abs((template.wavelength[aperture][overlap_range[0]] - \
					   template.wavelength[aperture][overlap_range[1]]) / \
					   float(len(tempdata)))
	
	# Rebin the science
	global ckms

	if dw_science < dw_template:
		if wavelengths[0] != template.wavelength[aperture][overlap_range[0]]:
			overlap_range[0] += 1
		if wavelengths[-1] != template.wavelength[aperture][overlap_range[1]]:
			overlap_range[1] -= 1
		
		try:
			tempdata = template.data[aperture][overlap_range[0]:overlap_range[1]+1]
			wavelengths = template.wavelength[aperture][overlap_range[0]:overlap_range[1]+1]
		except:
			tempdata = template.data[aperture][overlap_range[0]:-1]
			wavelengths = template.wavelength[aperture][overlap_range[0]:-1]

	'''		
	else:
		fit = interpolate.interp1d(template.wavelength[aperture], template.data[aperture])
		tempdata = fit(wavelengths)
	'''
	dv_average = ckms* (np.power(wavelengths[-1]/wavelengths[0], 1.0/float(len(wavelengths))) - 1.0)

	new_wavelengths = [wavelengths[0]]
	for i in range(1,len(wavelengths)):
		new_wavelengths.append(new_wavelengths[i-1]*dv_average/ckms + new_wavelengths[i-1])
	
	new_wavelengths = np.array(new_wavelengths)
	
	fit = interpolate.interp1d(science.wavelength[aperture], science.data[aperture])
	scidata = fit(new_wavelengths)
	fit = interpolate.interp1d(template.wavelength[aperture], template.data[aperture])
	tempdata = fit(new_wavelengths)

	
	allx = np.arange(len(new_wavelengths))
	# Knot spacing
	n_knots=40
	knots = len(allx)/n_knots #Denominator is approx number of knots
	#allx_knots = [knots/2. + np.mean(allx[i*knots:(i+1)*knots]) for i in range(len(allx)/knots)]
	allx_knots = [np.mean(allx[i*knots:(i+1)*knots+1]) for i in range((len(allx)-1)/knots)]
	# Trim off the last if they aren't associated with a knot.
	
	trim = abs(len(allx)%len(allx_knots) -1)
	if trim != 0:
		allx = allx[0:-trim]
		tempdata = tempdata[0:-trim]
		scidata = scidata[0:-trim]
	
	
	# Bin the data, then shift the "allx" by half a knot
	binned_template_knots = [np.mean(tempdata[i*knots:(i+1)*knots+1]) \
							 for i in range(len(allx_knots))]
							 #for i in range(len(tempdata)/knots)]
	
	# Tie the ends down so it doesn't read the maximum as the last point
	binned_spline = CubicSpline([allx[0]]+allx_knots+[allx[-1]],\
					 [tempdata[0]]+binned_template_knots+[tempdata[-1]])
	
	template_norm = tempdata / binned_spline(allx)
	
	# Repeat for the science spectrum
	science_knots = [np.mean(scidata[i*knots:(i+1)*knots]+1) \
					 for i in range(len(allx_knots))]
					 #for i in range(len(scidata)/knots)]
					 
	# Tie the ends down so it doesn't read the maximum as the last point
	binned_spline = CubicSpline([allx[0]]+allx_knots+[allx[-1]],\
					 [scidata[0]]+science_knots+[scidata[-1]])
	
	science_norm = scidata / binned_spline(allx)
	
	no_nans = list(set(np.where(~np.isnan(template_norm))[0]).intersection(set(np.where(~np.isnan(science_norm))[0])))
	science_norm = science_norm[no_nans]
	template_norm = template_norm[no_nans]
	wavelengths = new_wavelengths[no_nans]
		
	#soln = scipy.signal.correlate(template_norm, science_norm, mode="full")
	soln = np.correlate(template_norm, science_norm, mode="full")
	
	x = np.array(range(len(soln)), dtype=float)
		
	fit = curve_fit(triangle, x, soln,\
					p0=(len(soln)/2., 1.0, -1.0, len(soln)/2.)
					)[0]
	
	# TRY THIS
	soln_shrink = soln - triangle(range(len(soln)), *fit)

	# Find maximum in the correlation function
	# Trim 10% off the front and back just cuz
	solmax = np.where(soln_shrink==max(soln_shrink[len(soln)/10:-len(soln)/10]))[0][0]
	
	# Center it so that negative and positive correspond to pixel shifts
	# Have to add one since it goes from -(N-1) to 0 to +(N-1)
	# solmax is the location in the CCF of the maximum
	
	ccf_wavelengths = list(wavelengths-wavelengths[0])
	for w in wavelengths[1:]:
		ccf_wavelengths.insert(0,-1.*(w-wavelengths[0]))

	xmax = float(solmax) - float(len(science_norm) - 1.0)
	
	# r is the normalized centroid
	# s is the sigma from a gaussian fit to the centroid	
	# Fit a gaussian to the peak	
	center_peak = [soln_shrink[solmax]]
	
	center_x = [0.0]
	i = 0

	# only fit the center before it starts oscillating
	deriv = -1.
	while deriv <= 0:
		try:
			deriv = soln_shrink[solmax+i+1]-soln_shrink[solmax+i]
		except:
			import matplotlib.pyplot as plt
			plt.plot(soln)
			plt.plot(x, triangle(x, *fit))
			plt.plot(x, soln_shrink)
			plt.show()
		
		center_peak.append(soln_shrink[solmax+i+1])
		center_x.append(float(i+1))
		i += 1
	
	j = 0
	#deriv = soln_shrink[solmax-j]-soln_shrink[solmax-(j+1)]
	deriv = 1.0
	while deriv >= 0:
		deriv = soln_shrink[solmax-j]-soln_shrink[solmax-(j+1)]
		center_peak.insert(0, soln_shrink[solmax-(j+1)])
		center_x.insert(0, -float(j+1))
		j += 1
	
	if len(center_peak) <= 3:
		#return None, None, None, 0
		return None
	
	
	center_x = center_x[0:-1]
	center_peak = center_peak[0:-1]
	
	try:
		# Fit the centroid of the correlation peak with a gaussian
		fit = curve_fit(gauss, center_x, center_peak,\
						p0=(0.5, np.max(soln_shrink)-np.min(soln_shrink), 0.0, \
							np.mean([center_peak[0],center_peak[-1]])
							)
						)[0]
	except:
		#return None, None, None, 0
		return None
		
	s, a, mu, b = fit
	npoints = len(center_x)
	#chi2 = stats.chisquare(np.array(center_peak), np.array(gauss(center_x, *fit)))[0]
	
	# True peak height over average peak height
	siga_squared = [(soln_shrink[n+int(mu)] - soln_shrink[int(mu)-n])**2. \
					for n in range(int(len(wavelengths)-abs(mu)-1))]
	
	r = abs(gauss(mu, *fit)/np.sqrt(np.mean(siga_squared)))

	# Convert STDEV to FWHM
	w = 2.355*abs(s)
	# From RVSAO documentation:
	# https://iopscience.iop.org/article/10.1086/316207/pdf
	# adapted from Tonry & Davis (1979)
	sigma = (3./8.)*w/(1.+r)
		
	# There is an error here about "above interp range"
	centroid = xmax+mu
	
	# converts pixel location to wavelength
	winterp = interpolate.interp1d(range(len(wavelengths)), wavelengths)
	
	# If the centroid is negative, shift to blue, so the shift
	# should also be negative.
	s_centroid = np.sign(centroid)

	# winterp(centroid) asks where the centroid is in wavelength space
	shift = s_centroid*(wavelengths[0]-winterp(abs(centroid)))
	#centroid = abs(centroid)
		
	'''
	total_shift = []
	for i,w in enumerate(wavelengths):
		try:
			wrest = winterp(float(i)-centroid)
			dw = -1.0*(w-wrest)/wrest
			total_shift.append(ckms*dw)
		except:
			continue
	
	RV = sum(total_shift)/float(len(total_shift))
	'''
	#RV = ckms*shift/center_wavelength
	RV = ckms*shift/wavelengths[0]
	
	error = []
	
	total_high = []
	for i,w in enumerate(wavelengths):
		try:
			dw = (w - winterp(float(i)-(centroid+sigma)))/w
			total_high.append(ckms*dw)
		except:
			continue

	high = abs(sum(total_high)/float(len(total_high)))
	error.append(abs(high-abs(RV)))

	total_low = []
	for i,w in enumerate(wavelengths):
		try:
			dw = (w - winterp(float(i)-(centroid-sigma)))/w
			total_low.append(ckms*dw)
		except:
			continue

	low = abs(sum(total_low)/float(len(total_low)))
	error.append(abs(low-abs(RV)))
	
	'''
	try:
		# Convert centroid pixel+sigma to wavelength
		high = winterp(centroid+sigma)
		# Then shift 
		high_shift = abs(wavelengths[0] - high)
		print centroid, sigma, wavelengths[0], winterp(centroid+sigma)
		error.append(abs((ckms*high_shift/wavelengths[0])-RV))
	except:
		pass

	try:
		low = winterp(centroid-sigma)
		low_shift = abs(wavelengths[0] - low)
		error.append(abs(RV-ckms*low_shift/wavelengths[0]))
	except:
		pass
	'''

	error = np.average(error)
	
	if log!=None:
		log.write('Aperture %s\n\n' %aperture)
		log.write('\tGaussian fit ({:} points):\n'.format(npoints))
		log.write('\t\tcentroid  : {:>4.2f}\n'.format(centroid))
		log.write('\t\tsigma     : {:>4.2f}\n'.format(s))
		log.write('\t\tamplitude : {:>4.2f}\n'.format(a))
		log.write('\t\ty-offset  : {:>4.2f}\n'.format(b))

		log.write('\tT&D r-value : {:>4.2f}\n'.format(r))
		log.write('\tFWHM        : {:>4.2f}\n'.format(w))
		log.write('\tsigma (pix) : {:>4.2f}\n'.format(sigma))
		log.write('\tshift       : {:>4.2f}\n'.format(shift))
	
		log.write('\tDelta RV    : {:4.2f}\n'.format(RV))
		log.write('\tccf sigma   : {:4.2f}\n\n'.format(error))
		#log.write('-'*80 + '\n')
	
	if False:
		#x = center_x
		x = np.linspace(center_x[0]+xmax, center_x[-1]+xmax, 50)
		y = gauss(np.linspace(center_x[0], center_x[-1], 50), *fit)
		fig,axes = plt.subplots(2,1)
		print s_centroid*centroid, mu
		print gauss(mu, *fit), gauss(0.0, *fit), r, w, sigma, shift, RV
		
		axes[0].plot(soln_shrink)
		axes[0].plot(center_x+solmax,gauss(center_x,*fit))
		axes[0].set_xlim(solmax-10*s, 10*s+solmax)

		velocities = [-((wavelengths[p]/wavelengths[0])-1.0) for p in np.arange(len(wavelengths)-1,-1,-1)]
		velocities += list(np.array(velocities[::-1][1:])*-1)
		velocities = ckms*np.array(velocities)
		
		axes[1].plot(velocities, soln_shrink, color='gray', lw=1)
		
		# Gaussian x is still in pixelspace.
		fit = interpolate.interp1d(range(1-len(wavelengths),len(wavelengths)),\
								   velocities)
		
		axes[1].plot(fit(x), y, color='C0', ls=':')
		xmin = fit(min(x))
		xmax = fit(max(x))
		xwidth = xmax-xmin
		axes[1].set_xlim(xmin-(2*xwidth), xmax+(2*xwidth))
		
		axes[1].axvline(0, color='gray', lw=0.8, ls='--')
		#fit(x[np.where(y==max(y))[0][0]])
		axes[1].axvline(-RV,\
							color='C0', lw=0.8, ls='-')	
	
		plt.show()
		plt.close()
	
	velocities = [-((wavelengths[p]/wavelengths[0])-1.0) for p in np.arange(len(wavelengths)-1,-1,-1)]
	velocities += list(np.array(velocities[::-1][1:])*-1.0)
	velocities = ckms*np.array(velocities)
		
	# Gaussian x is still in pixelspace.
	wvfit = interpolate.interp1d(range(1-len(wavelengths),len(wavelengths)),\
							   velocities)
	
	x = np.linspace(center_x[0]+xmax, center_x[-1]+xmax, 50)
	xg= np.linspace(center_x[0], center_x[-1], 50)
	
	ccf = [velocities, soln_shrink, wvfit(x), gauss(xg, *fit)]
	comparison = [wavelengths, science_norm, template_norm, allx]
	
	return [aperture, RV, error, 1], ccf, comparison


# ==========================================================
# ==========================================================

def rv_by_aperture(template, science, aperture_list):
	#aperture_list = range(20)
	log = open('rv_ccf.log', 'a')
	
	log.write('='*80 + '\n')
	log.write('{:^80}\n'.format('RV cross-correlator'))

	import getpass,socket
	from datetime import datetime

	hostname = socket.gethostname()
	username = getpass.getuser()
	log.write('{:^80}\n\n'.format(username+'@'+hostname+' '+datetime.now().strftime("%d/%m/%Y %H:%M")))

	log.write('Object      : '+science.file_name+'\n')
	log.write('Template    : '+template.file_name+'\n')
	log.write('RV template : '+str(template.header['VHELIO'])+' km/s\n\n')
	
	'''
	rv_data = []
	for a in aperture_list:
		shift, RV, error, keep, ccf = cross_correlate(template, science, a, log=log)
		if RV!=None and error!=None:
			rv_data.append([a, RV, error, keep])
		else:
			print 'Not converged for aperture %s' %a

	rv_data = np.array(rv_data).T

	log.close()
	return rv_data
	'''	
	rv_data = []
	ccf = []
	comparison = []
	for a in aperture_list:
		result = cross_correlate(template, science, a, log=log)
		if result != None:
		#if RV!=None and error!=None:
			rv_data.append(result[0])
			ccf.append(result[1])
			comparison.append(result[2])
		else:
			print 'Not converged for aperture %s' %a

	rv_data = np.array(rv_data).T

	log.close()

	return rv_data, ccf, comparison
	
	
# ==========================================================
# ==========================================================

# array of ap, rv, err
def rv_average(data_list, sig_clip=2.0):
	median = round(np.median(data_list[1]),4)
	variances = []
	num = 0.
	keep = np.where(data_list[3]==1.0)[0]
	'''
	# I don't remember what this was for...
	for d in data_list[1,keep]:
		d = round(d,4)
		if d != median:
			variances.append((median-d)**2.)
			num += 1.
	
	if np.isnan(np.sqrt(np.sum(variances)/num)):
		total_variance = 500.
		print "  xxx  "
	else:
		total_variance = np.sqrt(np.sum(variances)/num)
	'''
	
	# Don't clip
	if sig_clip == 0 or sig_clip == None:
		avg = np.average(data_list[1,keep], weights=1./(data_list[2,keep]**2.0))

		# Average aperture uncertainty
		std1 = np.sqrt(sum([(err)**2. for err in data_list[2,keep]]))
		std1 = std1/np.sqrt(float(len(data_list[1,keep])))
	
		# Standard error of the mean
		std = np.sqrt(np.sum([(d-avg)**2. for d in data_list[1,keep]])/len(data_list[1,keep]))
		std = std/np.sqrt(len(data_list[1,keep])-1)
	
		#return data_list, avg, np.sqrt(std1**2. + std**2.), len(data_list[1])
		return data_list, avg, std1, len(data_list[1])
	
	for i,d in enumerate(data_list.T):
		if abs(d[1]-median) > sig_clip*d[2]:
			data_list[3,i] = 0
	
	# How many are kept now?
	keep = np.where(data_list[3]==1.0)[0]
	if len(keep)==0:
		return data_list, median, np.nan, 0
	
	keep = np.where(data_list[3]==1.0)[0]
	if len(keep)==1:
		return data_list, median, data_list[2,0], 1
	
	keep = np.where(data_list[3]==1.0)[0]
	if len(keep) != len(data_list[0]):
		data_list, avg, std, num = rv_average(data_list[:,keep], sig_clip)

	keep = np.where(data_list[3]==1.0)[0]
	weights = np.array([1.0/w**2. for w in data_list[2,keep]])
	avg = np.average(data_list[1,keep], 
					 weights=weights)
	
	# Average aperture uncertainty; sum the square of the error, divide
	# by the number of measurements, take the square root
	std1 = np.sqrt(np.sum(weights)/float(len(weights)))
	
	# Standard error of the mean
	std = np.sqrt(np.sum([(d-avg)**2. for d in data_list[1,keep]])/len(data_list[1,keep]))
	std = std/np.sqrt(len(data_list[1,keep])-1)
	
	#return data_list, avg, np.sqrt(std1**2. + std**2.), len(data_list[1,np.where(data_list[3]==1)[0]])
	return data_list, avg, std1, len(data_list[1,keep])

# ==========================================================
# ==========================================================

def logger():
	log = open('rv.log', 'a')
	return
	log.close()
	

def update_log(data):
	log.write(data)
	return