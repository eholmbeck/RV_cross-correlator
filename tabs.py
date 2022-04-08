try:
	import Tkinter as tk
except:
	import tkinter as tk
#import Tkinter as tk

from tkinter import ttk
from tkinter import messagebox

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
from matplotlib.lines import Line2D

import numpy as np
import spectrum
from tkinter.filedialog import askopenfilename,asksaveasfile
import rv

from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.time import Time
from astropy.constants import c


# =========================================================================

class Session:
	def __init__(self, root_window):	
		self.tab_parent = ttk.Notebook(root_window)
		tab1 = ttk.Frame(self.tab_parent)
		tab2 = ttk.Frame(self.tab_parent)
		tab3 = ttk.Frame(self.tab_parent)
		tab4 = ttk.Frame(self.tab_parent)
		self.tab_parent.add(tab1, text="Spectra")
		self.tab_parent.add(tab2, text="Cross-Correlate")
		self.tab_parent.add(tab3, text="Telluric Offset")
		self.tab_parent.add(tab4, text="Doppler Shift")

		self.start = Start(root_window, tab1)
		self.tellurics = Tellurics(root_window, tab3)
		self.corrections = Corrections(root_window, tab2, self.tellurics)
		self.doppler = Doppler(root_window, tab4)

		menubar_fmt = "{:<30}{:>2}"
		menubar = tk.Menu(root_window)
		filemenu = tk.Menu(menubar, tearoff=0)
		filemenu.add_command(label=menubar_fmt.format("New spectrum", u"\u2318"+"S"), command=self.open_science)
		root_window.bind('<Command-s>', self.open_science)
		filemenu.add_command(label=menubar_fmt.format("New template", u"\u2318"+"T"), command=self.open_template)
		root_window.bind('<Command-t>', self.open_template)
		filemenu.add_command(label=menubar_fmt.format("Export shifted spectrum", u"\u2318"+"E"), command=self.save_shifted_spectrum)
		root_window.bind('<Command-e>', self.save_shifted_spectrum)
		filemenu.add_separator()
		filemenu.add_command(label=menubar_fmt.format("Exit",""), command=root_window.quit)
		menubar.add_cascade(label="File", menu=filemenu)
		
		exportmenu = tk.Menu(menubar, tearoff=1)
		exportmenu.add_command(label=menubar_fmt.format("Export RV", u"\u2318"+"E"), command=self.write_log)
		root_window.bind('<Command-e>', self.write_log)
		menubar.add_cascade(label="Export", menu=exportmenu)

		root_window.config(menu=menubar)
		
	def open_science(self, event=None):
		filename = self.start.open_science_file()
		#self.start.science = spectrum.FITS(filename)
		self.corrections.science = self.start.science
		self.tellurics.science = self.start.science
		self.doppler.science = self.start.science
		return filename

	def open_template(self, event=None):
		filename = self.start.open_template_file()
		#self.start.template = spectrum.FITS(filename)
		self.corrections.template = self.start.template
		return filename
	
	# I think this is right?
	def save_shifted_spectrum(self, event=None):
		self.doppler.save_shifted_spectrum()

	def pack(self):
		self.tab_parent.pack(expand=1, fill="both")
	

	def write_log(self, event=None):		
		try:
			data_path = '/'.join(self.corrections.science.file_name.split('/')[:-1])
		except:
			print('ERROR: No measurements available to export!')
			return
		
		log_filename = data_path+'/RV_measurements.txt'
		from os.path import isfile
		if isfile(log_filename):
			log = open(log_filename, 'a')
		else:
			log = open(log_filename, 'w')
			log.write('Filename\tBCV\tRV\tRV_err\tTelluric offset\tAps\tAps_RV\tAps_err\n')
		
		keep = np.where(self.corrections.rv_data[3]==1.)[0]
		
		log.write(self.corrections.science.file_name.split('/')[-1]+'\t')
		log.write('{:>+.3f}'.format(self.corrections.bcv)+'\t')
		log.write('{:>+.3f}'.format(self.corrections.vhelio[0])+'\t')
		log.write('{:>.3f}'.format(self.corrections.vhelio[1])+'\t')
		log.write('{:>+.3f}'.format(self.tellurics.telluric_shift[0])+'\t')
		log.write(','.join(['{:>.0f}'.format(v) for v in self.corrections.rv_data[0][keep]])+'\t')
		log.write(','.join(['{:>.3f}'.format(v) for v in self.corrections.rv_data[1][keep]])+'\t')
		log.write(','.join(['{:>.3f}'.format(v) for v in self.corrections.rv_data[2][keep]])+'\n')
		log.close()
		
		messagebox.showinfo(title="Measurements Saved!", message="Measurements have been written to %s" %log_filename)

		return
	

# =========================================================================

class Start:
	def __init__(self, root_window, window):
		self.window = window
		self.root_window = root_window
		self.ap_num = None
		self.science = None
		self.template = None
		
		ws = self.window.winfo_screenwidth() # width of the screen
		hs = self.window.winfo_screenheight() # height of the screen
		dpi = self.window.winfo_fpixels('1i')
		
		self.fig = Figure(figsize=(hs*0.8/dpi, ws*0.5*0.8/dpi))
		#self.fig = Figure()
		self.ax = self.fig.add_subplot(111)
		self.template_ax = self.ax.twinx()
		
		self.spec_plot = FigureCanvasTkAgg(self.fig, self.window)
		self.spec_plot.get_tk_widget().grid(column=0, rowspan=4)

		toolbar_frame = tk.Frame(self.window)
		toolbar_frame.grid(column=0,rowspan=4,sticky="w")
		self.toolbar = NavigationToolbar2Tk(self.spec_plot, toolbar_frame)		

		self.aplabel = tk.Label(self.window, text="Aperture")
		self.aplabel.grid(column=1, row=2, sticky="e")
		
		self.ap_text = tk.StringVar()
		self.ap_text.set(1)
		self.ap_num = self.ap_text.get()
		#self.ap_text.trace("w", self.new_ap)
		self.ap_entry = tk.Entry(self.window, textvariable=self.ap_text)
		self.ap_entry.grid(column=2, row=2, sticky="w")
		self.new_ap
		self.ap_entry.bind("<Return>", self.new_ap)
		
		
		self.quit = tk.Button(self.window, text="Quit", command=quit)
		self.quit.grid(column=1, row=5)

		self.root_window.bind('<Right>', self.right)
		self.root_window.bind('<Left>', self.left)

	
	def new_ap(self, *args):
		self.ap_num = int(self.ap_text.get())
		self.ax.cla()
		self.template_ax.cla()
		self.plot_ap(self.science, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.template, self.ap_num, ax=self.template_ax, lw=1, zorder=2, color="gray")
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		return
	

	def plot_ap(self, spectrum, ap_num, ax=None, **kwargs):
		if ax is None: ax = self.ax
		if spectrum != None:
			try:
				ax.plot(spectrum.wavelength[ap_num], spectrum.data[ap_num], **kwargs)
			except:
				print("Aperture %s does not exist!" %ap_num)
		return

	def open_science_file(self):
		try:
			self.science_plot.remove()
		except:
			pass

		filename = askopenfilename(filetypes=(("fits files", "*.fits"), ("All files", "*,*")))
		self.science = spectrum.FITS(filename)
		self.ap_num = self.science.first_beam
		self.science_plot = self.plot_ap(self.science, self.science.first_beam, lw=1)
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.science.first_beam),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		
		#self.ap_num = self.science.first_beam
		self.ap_text.set(self.ap_num)
		return filename

	def open_template_file(self):
		try:
			self.template_plot.remove()
		except:
			pass

		filename = askopenfilename(filetypes=(("fits files", "*.fits"), ("All files", "*,*")))
		self.template = spectrum.FITS(filename)
		self.vhelio_template_exists()
		if self.template.first_beam > self.science.first_beam:
			self.ap_num = self.template.first_beam
			self.science_plot = self.plot_ap(self.science, self.ap_num, lw=1)
			self.template_plot = self.plot_ap(self.template, self.ap_num, ax=self.template_ax, lw=1, color="gray")

		else:
			self.template_plot = self.plot_ap(self.template, self.science.first_beam, ax=self.template_ax, lw=1, color="gray")

		self.spec_plot.draw()
		return filename
	

	def vhelio_template_exists(self):
		vh = None
		try:
			vh = float(self.template.header["VHELIO"])
			return
			
		except:
			answer = messagebox.askyesno(title="Warning", message="Template file has no VHELIO in fits header. Would you like to add one now?")
			if answer == True:
				try: import tkSimpleDialog as simpledialog
				except ModuleNotFoundError: import tkinter.simpledialog as simpledialog
				vh = simpledialog.askfloat("", "Enter the heliocentric velocity:",
                               parent=self.root_window)

				from astropy.io import fits
				#newtemplate = self.template.file_name.replace(".fits", "")+"_tmp.fits"
				data, hdr = fits.getdata(self.template.file_name, header=True)
				hdul = fits.HDUList()
				new_header = self.template.header
				new_header["VHELIO"] = '{:.3f}'.format(vh)
		
				fits.writeto(self.template.file_name, data, new_header, overwrite=True)
				hdul.close()
		
		print("Using VHELIO =", "{:.3f}".format(vh), "for template spectrum.")
		return
				

	def right(self,event):
		if self.ap_num == self.science.first_beam + self.science.apertures:
			return
		
		self.ap_num += 1
		self.ap_text.set(self.ap_num)
		self.ax.cla()
		self.template_ax.cla()
		self.plot_ap(self.science, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.template, self.ap_num, ax=self.template_ax, lw=1, zorder=2, color="gray")
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		
	def left(self,event):
		if self.ap_num == self.science.first_beam:
			return
		
		self.ap_num -= 1
		self.ap_text.set(self.ap_num)
		self.ax.cla()
		self.template_ax.cla()
		self.plot_ap(self.science, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.template, self.ap_num, ax=self.template_ax, lw=1, zorder=2, color="gray")
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()

# =========================================================================
#rv_data = rv.rv_by_aperture(template, science, range(science.apertures))
#new_data_list, avg, std, n = rv.rv_average(rv_data, 0.2)


class Corrections:
	def __init__(self, root_window, window, telluric_tab):
		self.window = window
		self.root_window = root_window
		self.science = None
		self.template = None
		self.ccf = []
		
		ws = self.window.winfo_screenwidth() # width of the screen
		hs = self.window.winfo_screenheight() # height of the screen
		dpi = self.window.winfo_fpixels('1i')

		self.fig = Figure(figsize=(hs*0.8/dpi, ws*0.5*0.8/dpi))
		self.rv_ax = self.fig.add_subplot(2,1,1,picker=True)
		self.ccf_ax = self.fig.add_subplot(2,1,2,picker=True)
		self.ccf_ax.set_xlabel('Velocity Shift (km/s)')
		self.gaussian_point_chosen = None
		self.chosen_row = None
		
		self.rv_plot = FigureCanvasTkAgg(self.fig, self.window)
		self.rv_plot.get_tk_widget().grid(column=0, row=0, rowspan=10)

		toolbar_frame = tk.Frame(self.window)
		toolbar_frame.grid(column=0,rowspan=10, sticky="w")
		self.toolbar = NavigationToolbar2Tk(self.rv_plot, toolbar_frame)		
		
		self.cbutton = tk.Button(self.window, text="Run cross-correlation", 
				command=self.calculate_rvs)
		self.cbutton.grid(column=1, row=0, columnspan=2)

		self.zbutton = tk.Button(self.window, text="Set HRV to 0", 
				command=self.set_to_zero)
		self.zbutton.grid(column=1, row=1, columnspan=2)
		
		self.sglabel = tk.Label(self.window, text="Outlier sigma:")
		self.sglabel.grid(column=1, row=2, sticky="e")
		self.sg_text = tk.StringVar()
		self.sg_text.set(0)
		self.sig_clip = float(self.sg_text.get())
		
		self.sg_entry = tk.Entry(self.window, textvariable=self.sg_text, width=8)
		self.sg_entry.grid(column=2, row=2, sticky="w")
		self.new_sg
		self.sg_entry.bind("<Return>", self.new_sg) 
				
		self.rv_data = None
		self.keep = None
		self.new_data_list = None
		self.vhelio = None
		self.vm_template = None
		self.telluric_tab = telluric_tab
		
		self.vh = tk.StringVar()
		self.vhlabel = tk.Label(self.window, text="RV helio:")
		self.vhlabel.grid(column=1, row=3, sticky="e")
		self.vhbox = tk.Label(self.window, textvariable=self.vh)
		self.vhbox.grid(column=2, row=3, sticky="w")
			
		self.bcv_text = tk.StringVar()
		self.bcvlabel = tk.Label(self.window, text="BCV:")
		self.bcvlabel.grid(column=1, row=4, sticky="e")
		'''
		self.bcvbox = tk.Label(self.window, textvariable=self.bcv_text)
		self.bcvbox.grid(column=2, row=4, sticky="w")
		'''
		#self.bcv_float = float(self.bcv_text.get())
		self.bcv_entry = tk.Entry(self.window, textvariable=self.bcv_text, width=8)
		self.bcv_entry.grid(column=2, row=4, sticky="w")
		self.bcv_entry.bind("<Return>", self.set_bcv) 


		self.fig.canvas.mpl_connect('pick_event', self.on_click)
		self.fig.canvas.mpl_connect('key_press_event', self.figure_key_press)
		
		self.rv_table_box = tk.Frame(self.window)
		self.rv_table_box.grid(row=5, column=1, columnspan=2, sticky="nw")
		self.rv_checks = None

		self.rv_frame = tk.Canvas(self.rv_table_box, borderwidth=2)
		self.rv_frame.bind_all("<MouseWheel>", self.on_mousewheel)


	def set_bcv(self, event):
		self.bcv = float(self.bcv_text.get())
		self.bcv_text.set("{:.2f}".format(self.bcv))
		
		vm = self.vm_template + self.delta_rv[1]
		ckms = c.to(u.km/u.s).value
		vmvb_c = vm*self.bcv/ckms
		self.rv_data[1] = vm + self.bcv + vmvb_c
		
		#self.rv_plot.focus()
		self.rv_table()
		self.rv_average(sig_clip=0)
		

	def bcv_corr(self, spectrum):
		if False:
		#if "BCV" in spectrum.header.keys():
			return spectrum.header["BCV"]
		
		else:
			try:
				site_name = spectrum.header["OBSERVAT"].lower()
			except:
				site_name = spectrum.header["SITENAME"].lower()
			
			site = EarthLocation.of_site(site_name)
			
			if 'lco' in site_name:
				RA = spectrum.header["RA-D"]
				DEC = spectrum.header["DEC-D"]
				sc = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg)
				UT_DATE = spectrum.header['UT-DATE']
				UT_TIME = spectrum.header["UT-START"]
				ut_start = Time(UT_DATE + ' ' + UT_TIME, format='iso', scale="utc")
				ut_mid = ut_start + (spectrum.header['EXPTIME']/2.)*u.s

			elif 'mcdonald' in site_name:
				RA = spectrum.header['RA']
				DEC = spectrum.header['DEC']
				sc = SkyCoord(ra=RA, dec=DEC, unit=(u.hourangle, u.deg))
				UT_DATE = spectrum.header['DATE-OBS']
				UT_TIME = spectrum.header["UT"]
				ut_start = Time(UT_DATE + ' ' + UT_TIME, format='iso', scale="utc")
				ut_mid = ut_start + (spectrum.header['EXPTIME']/2.)*u.s

			elif 'apo' in site_name:
				RA = spectrum.header['RA']
				DEC = spectrum.header['DEC']
				sc = SkyCoord(ra=RA, dec=DEC, unit=(u.hourangle, u.deg))
				UT = spectrum.header['DATE-OBS']
				ut_start = Time(UT, format='isot', scale="tai")
				ut_mid = ut_start + (spectrum.header['EXPTIME']/2.)*u.s
			
			#TODO: Generalize this!
			# -------------------------------------------------------------
			else:
				try:
					site = EarthLocation.of_site(spectrum.header["OBSERVAT"]) 
				
				except:
					site = EarthLocation.from_geodetic(lat=spectrum.header["SITELAT"]*u.deg, 
							lon=spectrum.header["SITELONG"]*u.deg,
							height=spectrum.header["SITEALT"]*u.m)

				try:
					RA = spectrum.header["RA-D"]
					DEC = spectrum.header["DEC-D"]
					sc = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg)
				except:
					RA = spectrum.header['RA']
					DEC = spectrum.header['DEC']
					sc = SkyCoord(ra=RA, dec=DEC, unit=(u.hourangle, u.deg))

				# TO-DO: CHANGE TO UT-MID
				try:
					UT_DATE = spectrum.header['DATE-OBS']
				except:
					UT_DATE = spectrum.header['UT-DATE']
			
				try:
					UT_TIME = spectrum.header["UT"]
				except TypeError:
					UT_TIME = spectrum.header["UT-TIME"]
				
				ut_start = Time(UT_DATE + ' ' + UT_TIME, format='iso', scale="utc")
				ut_mid = ut_start + (spectrum.header['EXPTIME']/2.)*u.s
			# -------------------------------------------------------------			
			
			barycorr = sc.radial_velocity_correction(obstime=ut_mid, location=site) 
			bcv = barycorr.to(u.km/u.s).value
			return bcv

	
	def calculate_rvs(self):
		rv_data, telluric_data = rv.rv_by_aperture(\
										self.template, self.science, \
										self.science.data.keys(), \
										tellurics=True)
		self.delta_rv, self.ccf = rv_data
		self.rv_data = np.copy(self.delta_rv)

		self.telluric_tab.set_telluric_data(telluric_data)

		self.correct_rvs()
		
		vm = self.vm_template + self.delta_rv[1]
		ckms = c.to(u.km/u.s).value
		vmvb_c = vm*self.bcv/ckms
		self.rv_data[1] = vm + self.bcv + vmvb_c
		
		self.rv_table()
		self.rv_average(sig_clip=0)
	
	
	def set_to_zero(self):
		#self.vhelio = [0.0,0.0]
		# Set DRV=0; then test.fits (combined) is about vm off from v_helio.
		self.delta_rv[1] *= 0.0
		vm = self.vm_template + self.delta_rv[1]
		ckms = c.to(u.km/u.s).value
		vmvb_c = vm*self.bcv/ckms
		#self.rv_data[1] = vm + self.bcv + vmvb_c
		#self.rv_data[1] = vm + vmvb_c
		self.rv_data[1] *= 0.0
		
		self.rv_table()
		self.rv_average(sig_clip=0)
		
		return
	
	
	def rv_average(self, sig_clip=None):
		if sig_clip == None:
			sig_clip = self.sig_clip
		rv_data, avg, std, n = rv.rv_average(self.rv_data, sig_clip)
		
		for i,r in enumerate(rv_data.T):
			aploc = np.where(self.rv_data[0]==r[0])[0][0]
			if r[3] == 1:
				self.rv_data[3,aploc] = 1
				self.delta_rv[3,aploc] = 1
			else:
				self.rv_data[3,aploc] = 0
				self.delta_rv[3,aploc] = 0
		
		for i in range(len(self.rv_data[3])):
			self.rv_checks[i].set(int(self.rv_data[3,i]))
			
		self.vhelio = [avg,std]
		self.keep = np.where(self.rv_data[3]==1)[0]
		
		ckms = c.to(u.km/u.s).value
		self.science.header["DOPCOR"] = (self.vhelio[0] - self.bcv)/(1.+self.bcv/ckms)
		self.science.header["VHELIO"] = self.vhelio[0]
		self.vh.set("{:.2f}".format(self.vhelio[0])+u" \u00B1 "+"{:.2f}".format(self.vhelio[1]))
		
		# Add table
		self.plot_results()
		return self.rv_ax
	
	
	def correct_rvs(self):
		# RV + correction
		self.bcv = self.bcv_corr(self.science)
		self.bcv_text.set("{:.2f}".format(self.bcv))
		
		bcv_template = self.bcv_corr(self.template)
		ckms = c.to(u.km/u.s).value
		self.vm_template = (float(self.template.header["VHELIO"]) - bcv_template)/(1.+(bcv_template/ckms))
		#correction = float(self.template.header["VHELIO"]) - bcv_template \
		#		+ bcv_science
		return
		
		
	def plot_results(self):
		avg,std = self.vhelio
		self.rv_ax.cla()
		self.rv_ax.text(0.05, 0.95, "%s apertures used" %len(self.keep),\
					va="top", transform=self.rv_ax.transAxes)
		
		self.rv_ax.errorbar(self.rv_data[0], self.rv_data[1], yerr=self.rv_data[2],\
						ls="", markerfacecolor="white", marker="o", color="C0",\
						alpha=0.7, picker=5)
		
		self.rv_ax.scatter(self.rv_data[0,self.keep], \
						self.rv_data[1,self.keep], color="C0",\
						marker="o", zorder=4, picker=5)
		
		self.rv_ax.plot([self.science.first_beam,self.science.first_beam+self.science.apertures], \
					 [avg, avg], color="C0", ls="--")
		
		self.rv_ax.set_ylim(np.min(self.rv_data[1,self.keep]-2.5*(self.rv_data[2,self.keep])),\
						 np.max(self.rv_data[1,self.keep]+2.5*(self.rv_data[2,self.keep]))) 

		self.rv_plot.draw()

		return
		
		
	def new_sg(self, *args):
		self.sig_clip = float(self.sg_text.get())
		self.rv_average()
		return


	def on_click(self,event):
		if isinstance(event.artist, Line2D):
			ind = event.ind[0]
			
			x_move = event.artist.get_xdata()[ind]
			loc = np.where(self.rv_data[0]==x_move)[0][0]
			self.plot_ccf(loc)
		
		return			

	def on_mousewheel(self, event):
		self.rv_frame.yview_scroll(-1*(event.delta), "units")
			
	def rv_table(self):
		#rv_frame.pack(side="left", fill="both", padx=10, pady=10)
		self.rv_frame.grid(row=0, column=0, sticky="news")
		buttons_frame = tk.Frame(self.rv_frame, borderwidth=2, relief="groove")
		buttons_frame.grid(row=0, column=0, sticky="news")

		ROWS = len(self.rv_data[0])
		ROWS_DISP = 24
		
		ttk.Label(buttons_frame, text="Use").grid(row=0,column=0, sticky="nwe")
		ttk.Label(buttons_frame, text="Aperture", width=9).grid(row=0,column=1, sticky="new")
		ttk.Label(buttons_frame, text="RV", width=8).grid(row=0,column=2, sticky="nwe")
		ttk.Label(buttons_frame, text="RV err", width=8).grid(row=0,column=3, sticky="nwe")
		
		self.rv_checks = [tk.IntVar(value=int(r[3])) for r in self.rv_data.T]
		for i,r in enumerate(self.rv_data.T):
			cb = ttk.Checkbutton(buttons_frame, variable=self.rv_checks[i], 
								 command=lambda i=i: self.select_aps(i),
								 takefocus=False)
			cb.grid(column=0, row=i+1)
			#self.rv_checks[i].set(int(r[3]))
			
			ttk.Label(buttons_frame, text="{:>5.0f}".format(r[0])).grid(column=1, row=i+1, sticky='news')
			ttk.Label(buttons_frame, text="{:>+6.2f}".format(r[1])).grid(column=2, row=i+1, sticky='news')
			ttk.Label(buttons_frame, text="{:>6.2f}".format(r[2])).grid(column=3, row=i+1, sticky='news')
			#boxes[x].append(Checkbutton(master, variable = boxVars[x][y], command = lambda x = x: checkRow(x)))
			#boxes[x][y].grid(row=x+1, column=y+1)
		
		rv_scroll = ttk.Scrollbar(self.rv_table_box, orient="vertical", command=self.rv_frame.yview)
		#rv_scroll.pack(side="right", fill="y")
		rv_scroll.grid(row=0, column=0, sticky="nes")
		self.rv_frame.configure(yscrollcommand=rv_scroll.set)#, scrollregion=rv_frame.bbox("all"))
		
		self.rv_frame.create_window((0,0), window=buttons_frame, anchor=tk.NW)
		buttons_frame.update_idletasks()
		#bbox = self.rv_table_box.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.
		bbox = self.rv_frame.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.

		# Define the scrollable region as entire canvas with only the desired
		# number of rows and columns displayed.
		w, h = bbox[2]-bbox[1], bbox[3]-bbox[1]
		dw, dh = int((w/4) * 4), int((h/ROWS) * ROWS_DISP)
		self.rv_frame.configure(scrollregion=bbox, width=w, height=dh)
	
		return
	
	def select_aps(self, row):
		self.rv_checks[row].set(int(abs(self.rv_data[3,row]%2-1)))
		self.rv_data[3,row] = abs(self.rv_data[3,row]%2-1)
		
		self.keep = np.where(self.rv_data[3]==1)[0]
		self.rv_average(sig_clip=0)

	def plot_ccf(self, row):
		row = int(row)
		self.ccf_ax.cla()
		if row != self.chosen_row:
			self.gaussian_point_chosen = None
			self.chosen_row = row
		
		self.ccf_ax.plot(self.ccf[row][0], self.ccf[row][1], color='gray', lw=1)
		self.ccf_ax.plot(self.ccf[row][2], self.ccf[row][3], color='C0', ls=':')
		xmin = min(self.ccf[row][2])
		xmax = max(self.ccf[row][2])
		xwidth = xmax-xmin
		self.ccf_ax.set_xlim(xmin-max(xwidth,80), xmax+max(xwidth,80))
		#self.ccf_ax.set_xlim(-350, 350)
		
		self.ccf_ax.text(0.05, 0.95, "Aperture %s" %int(self.rv_data[0][row]),\
					va="top", transform=self.ccf_ax.transAxes)
		self.ccf_ax.axvline(0, color='k', lw=2, alpha=0.5)

		# Recover old RV
		ckms = c.to(u.km/u.s).value
		vbc1 = 1.0+(self.bcv/ckms)
		#doppler_shift = self.rv_data[1][row] - self.vm_template*vbc1 - self.bcv
		doppler_shift = self.delta_rv[1][row]
		#self.ccf_ax.axvline(-doppler_shift/vbc1,\
		self.ccf_ax.axvline(-doppler_shift,\
							color='C0', lw=0.8, ls='-')
		
		self.rv_plot.draw()

	
	def figure_key_press(self, event):
		if event.key in "gk":
			x_select = event.xdata #event.artist.get_xdata()[ind]
			if self.gaussian_point_chosen == None:
				self.gaussian_point_chosen = x_select
				print('Ap. {:.0f}: Point 1 selected for refitting at {:.1f}.'.format(int(self.rv_data[0][self.chosen_row]), x_select))
			else:
				self.refit_CCF(x_select)
				self.gaussian_point_chosen = None
		
		return

	
	def refit_CCF(self, lambda2):
		from scipy.optimize import curve_fit
		
		lambda1 = self.gaussian_point_chosen
		if lambda1 > lambda2:
			x2temp = lambda1
			lambda1 = lambda2
			lambda2 = x2temp
		
		# x chosen by velocity space, must we go back to pixels?
		# If the data are decreasing...
		if self.ccf[self.chosen_row][0][0] > self.ccf[self.chosen_row][0][1]:
			x1 = np.where(self.ccf[self.chosen_row][0]>=lambda2)[0][-1]
			x2 = np.where(self.ccf[self.chosen_row][0]<=lambda1)[0][0]

		else:
			x1 = np.where(self.ccf[self.chosen_row][0]<=lambda1)[0][-1]
			x2 = np.where(self.ccf[self.chosen_row][0]>=lambda2)[0][0]
		
		x = self.ccf[self.chosen_row][0][x1:x2]
		y = self.ccf[self.chosen_row][1][x1:x2]

		if len(x) < 4:
			self.gaussian_point_chosen = None
			print('Error! Try choosing points again.')
			return

		def gauss(x, sigma, a, mu, b):
			g = np.exp(-0.5*((x-mu)/sigma)**2)
			return a*g + b
		
		p0 = (5.0, np.max(y)-np.min(y), np.mean([lambda1,lambda2]), np.min(y))
		fit = curve_fit(gauss, x, y, p0=p0)[0]
		
		s, a, mu, b = fit
		siga_squared = [(self.ccf[self.chosen_row][1][n+int(mu)] - self.ccf[self.chosen_row][1][int(mu)-n])**2. \
						for n in range(int(len(self.ccf[self.chosen_row][1])-abs(mu)-1))]
	
		r = abs(gauss(mu, *fit)/np.sqrt(np.mean(siga_squared)))
		w = 2.355*abs(s)
		sigma = (3./8.)*w/(1.+r)
		
		ckms = c.to(u.km/u.s).value
		#current_dRV = self.rv_data[1,self.chosen_row]
		#vb_c1 = 1.0 + (self.bcv/ckms)
		vmvb_c = (self.vm_template - mu)*self.bcv/ckms
		
		self.rv_data[1,self.chosen_row] = -mu + self.vm_template + self.bcv + vmvb_c
		self.rv_data[2,self.chosen_row] = sigma
		self.rv_data[3,self.chosen_row] = 1.0

		self.delta_rv[1,self.chosen_row] = -mu
		self.delta_rv[2,self.chosen_row] = sigma
		self.delta_rv[3,self.chosen_row] = 1.0
		
		fitx = np.linspace(lambda1,lambda2,50)
		self.ccf[self.chosen_row][2] = fitx
		self.ccf[self.chosen_row][3] = gauss(fitx, *fit)
		
		self.plot_ccf(self.chosen_row)
		self.rv_average(sig_clip=0)
		# Update value in table
		# TODO: Make this update, not rebuild the table from scratch!
		#self.rv_table()
		
		return
	
	
# =========================================================================
# This tab find the Telluric offset

class Tellurics:
	def __init__(self, root_window, window):
		self.window = window
		self.root_window = root_window
		self.science = None
		self.template = None
		self.ccf = []
		
		ws = self.window.winfo_screenwidth() # width of the screen
		hs = self.window.winfo_screenheight() # height of the screen
		dpi = self.window.winfo_fpixels('1i')

		self.fig = Figure(figsize=(hs*0.8/dpi, ws*0.5*0.8/dpi))
		self.rv_ax = self.fig.add_subplot(2,1,1,picker=True)
		self.ccf_ax = self.fig.add_subplot(2,1,2,picker=True)
		self.ccf_ax.set_xlabel('Velocity Shift (km/s)')
		self.gaussian_point_chosen = None
		self.chosen_row = None
		
		self.rv_plot = FigureCanvasTkAgg(self.fig, self.window)
		self.rv_plot.get_tk_widget().grid(column=0, row=0, rowspan=6)

		toolbar_frame = tk.Frame(self.window)
		toolbar_frame.grid(column=0,rowspan=6, sticky="w")
		self.toolbar = NavigationToolbar2Tk(self.rv_plot, toolbar_frame)		
		
		self.cbutton = tk.Button(self.window, text="Apply telluric offset", 
				command=self.telluric_offset)
		self.cbutton.grid(column=1, row=0, columnspan=2)
		
		self.keep = None
		self.new_data_list = None
		self.telluric_shift = None
		
		self.vh = tk.StringVar()
		self.vhlabel = tk.Label(self.window, text="Offset:")
		self.vhlabel.grid(column=1, row=2, sticky="e")
		self.vhbox = tk.Label(self.window, textvariable=self.vh)
		self.vhbox.grid(column=2, row=2, sticky="w")

		self.vth = tk.StringVar()
		self.vthlabel = tk.Label(self.window, text="Total RV:")
		self.vthlabel.grid(column=1, row=3, sticky="e")
		self.vthbox = tk.Label(self.window, textvariable=self.vth)
		self.vthbox.grid(column=2, row=3, sticky="w")
				
		self.fig.canvas.mpl_connect('pick_event', self.on_click)
		self.fig.canvas.mpl_connect('key_press_event', self.figure_key_press)
		
		self.rv_table_box = tk.Frame(self.window)
		self.rv_table_box.grid(row=4, column=1, columnspan=2, sticky="nw")
		self.rv_checks = None

		self.rv_frame = tk.Canvas(self.rv_table_box, borderwidth=2)
		self.rv_frame.bind_all("<MouseWheel>", self.on_mousewheel)

	
	def telluric_offset(self):
		ckms = c.to(u.km/u.s).value
		self.science.header["DOPCOR"] -= self.telluric_shift[0]
		self.science.header["VHELIO"] -= self.telluric_shift[0]
		self.vth.set("{:.2f}".format(self.science.header["VHELIO"]))

		
	def set_telluric_data(self, telluric_data):
		self.telluric_data, self.ccf = telluric_data
		self.rv_table()
		self.rv_average(sig_clip=0)
		
	
	def rv_average(self, sig_clip=None):
		if sig_clip == None:
			sig_clip = self.sig_clip
		rv_data, avg, std, n = rv.rv_average(self.telluric_data, sig_clip)
		
		for i,r in enumerate(rv_data.T):
			aploc = np.where(self.telluric_data[0]==r[0])[0][0]
			if r[3] == 1.0:
				self.telluric_data[3,aploc] = 1.0
			else:
				self.telluric_data[3,aploc] = 0.0
		
		for i in range(len(self.telluric_data[3])):
			self.rv_checks[i].set(int(self.telluric_data[3,i]))
			
		self.telluric_shift = [avg,std]
		self.keep = np.where(self.telluric_data[3]==1)[0]
		
		self.vh.set("{:.2f}".format(self.telluric_shift[0])+u" \u00B1 "+"{:.2f}".format(self.telluric_shift[1]))
		try: self.vth.set("{:.2f}".format(self.science.header["VHELIO"]  - self.telluric_shift[0]))
		except: pass
		
		# Add table
		self.plot_results()
		return self.rv_ax
	
		
	def plot_results(self):
		avg,std = self.telluric_shift
		self.rv_ax.cla()
		self.rv_ax.text(0.05, 0.95, "%s apertures used" %len(self.keep),\
					va="top", transform=self.rv_ax.transAxes)
		
		self.rv_ax.errorbar(self.telluric_data[0], self.telluric_data[1], yerr=self.telluric_data[2],\
						ls="", markerfacecolor="white", marker="o", color="C0",\
						alpha=0.7, picker=5)
		
		self.rv_ax.scatter(self.telluric_data[0,self.keep], \
						self.telluric_data[1,self.keep], color="C0",\
						marker="o", zorder=4, picker=5)
		
		self.rv_ax.plot([self.science.first_beam,self.science.first_beam+self.science.apertures], \
					 [avg, avg], color="C0", ls="--")
		
		self.rv_ax.set_ylim(np.min(self.telluric_data[1,self.keep]-2.5*(self.telluric_data[2,self.keep])),\
						 np.max(self.telluric_data[1,self.keep]+2.5*(self.telluric_data[2,self.keep]))) 

		self.rv_plot.draw()

		return
		

	def on_click(self,event):
		if isinstance(event.artist, Line2D):
			ind = event.ind[0]
			
			x_move = event.artist.get_xdata()[ind]
			loc = np.where(self.telluric_data[0]==x_move)[0][0]
			self.plot_ccf(loc)
		
		return			
			
	def on_mousewheel(self, event):
		self.rv_frame.yview_scroll(-1*(event.delta), "units")

	def rv_table(self):
		self.rv_frame.grid(row=0, column=0, sticky="news")
		self.rv_frame.grid(row=0, column=0, sticky="news")
		buttons_frame = tk.Frame(self.rv_frame, borderwidth=2, relief="groove")
		buttons_frame.grid(row=1, column=0, sticky="news")

		ROWS = len(self.telluric_data[0])
		ROWS_DISP = 24
		
		ttk.Label(buttons_frame, text="Use").grid(row=0,column=0, sticky="nwe")
		ttk.Label(buttons_frame, text="Aperture", width=9).grid(row=0,column=1, sticky="new")
		ttk.Label(buttons_frame, text="RV", width=8).grid(row=0,column=2, sticky="nwe")
		ttk.Label(buttons_frame, text="RV err", width=8).grid(row=0,column=3, sticky="nwe")
		
		self.rv_checks = [tk.IntVar(value=int(r[3])) for r in self.telluric_data.T]
			
		for i,r in enumerate(self.telluric_data.T):
			cb = ttk.Checkbutton(buttons_frame, variable=self.rv_checks[i], 
								 command=lambda i=i: self.select_aps(i),
								 takefocus=False)
			cb.grid(column=0, row=i+1)
			#self.rv_checks[i].set(int(r[3]))
			
			ttk.Label(buttons_frame, text="{:>5.0f}".format(r[0])).grid(column=1, row=i+1, sticky='news')
			ttk.Label(buttons_frame, text="{:>+6.2f}".format(r[1])).grid(column=2, row=i+1, sticky='news')
			ttk.Label(buttons_frame, text="{:>6.2f}".format(r[2])).grid(column=3, row=i+1, sticky='news')
			#boxes[x].append(Checkbutton(master, variable = boxVars[x][y], command = lambda x = x: checkRow(x)))
			#boxes[x][y].grid(row=x+1, column=y+1)
		
		rv_scroll = ttk.Scrollbar(self.rv_table_box, orient="vertical", command=self.rv_frame.yview)
		#rv_scroll.pack(side="right", fill="y")
		rv_scroll.grid(row=0, column=0, sticky="nes")
		self.rv_frame.configure(yscrollcommand=rv_scroll.set)#, scrollregion=rv_frame.bbox("all"))
		
		self.rv_frame.create_window((0,0), window=buttons_frame, anchor=tk.NW)
		buttons_frame.update_idletasks()
		#bbox = self.rv_table_box.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.
		bbox = self.rv_frame.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.

		# Define the scrollable region as entire canvas with only the desired
		# number of rows and columns displayed.
		w, h = bbox[2]-bbox[1], bbox[3]-bbox[1]
		dw, dh = int((w/4) * 4), int((h/ROWS) * ROWS_DISP)
		self.rv_frame.configure(scrollregion=bbox, width=w, height=dh)

		return

	
	def select_aps(self, row):
		self.rv_checks[row].set(int(abs(self.telluric_data[3,row]%2-1)))
		self.telluric_data[3,row] = abs(self.telluric_data[3,row]%2-1)
		
		self.keep = np.where(self.telluric_data[3]==1)[0]
		self.rv_average(sig_clip=0)

	def plot_ccf(self, row):
		row = int(row)
		self.ccf_ax.cla()
		if row != self.chosen_row:
			self.gaussian_point_chosen = None
			self.chosen_row = row
		
		self.ccf_ax.plot(self.ccf[row][0], self.ccf[row][1], color='gray', lw=1)
		self.ccf_ax.plot(self.ccf[row][2], self.ccf[row][3], color='C0', ls=':')
		xmin = min(self.ccf[row][2])
		xmax = max(self.ccf[row][2])
		xwidth = xmax-xmin
		self.ccf_ax.set_xlim(xmin-max(xwidth,80), xmax+max(xwidth,80))
		#self.ccf_ax.set_xlim(-350, 350)
		
		self.ccf_ax.text(0.05, 0.95, "Aperture %s" %int(self.telluric_data[0][row]),\
					va="top", transform=self.ccf_ax.transAxes)
		self.ccf_ax.axvline(0, color='k', lw=2, alpha=0.5)

		# Recover old RV
		self.ccf_ax.axvline(-self.telluric_data[1][row], color='C0', lw=0.8, ls='-')

		self.rv_plot.draw()

	
	def figure_key_press(self, event):
		if event.key in "gk":
			x_select = event.xdata #event.artist.get_xdata()[ind]
			if self.gaussian_point_chosen == None:
				self.gaussian_point_chosen = x_select
				print('Ap. {:.0f}: Point 1 selected for refitting at {:.1f}.'.format(int(self.telluric_data[0][self.chosen_row]), x_select))
			else:
				self.refit_CCF(x_select)
				self.gaussian_point_chosen = None
		
		return

	
	def refit_CCF(self, lambda2):
		from scipy.optimize import curve_fit
		
		lambda1 = self.gaussian_point_chosen
		if lambda1 > lambda2:
			x2temp = lambda1
			lambda1 = lambda2
			lambda2 = x2temp
		
		# x chosen by velocity space, must we go back to pixels?
		# If the data are decreasing...
		if self.ccf[self.chosen_row][0][0] > self.ccf[self.chosen_row][0][1]:
			x1 = np.where(self.ccf[self.chosen_row][0]>=lambda2)[0][-1]
			x2 = np.where(self.ccf[self.chosen_row][0]<=lambda1)[0][0]

		else:
			x1 = np.where(self.ccf[self.chosen_row][0]<=lambda1)[0][-1]
			x2 = np.where(self.ccf[self.chosen_row][0]>=lambda2)[0][0]
		
		x = self.ccf[self.chosen_row][0][x1:x2]
		y = self.ccf[self.chosen_row][1][x1:x2]

		if len(x) < 4:
			self.gaussian_point_chosen = None
			print('Error! Try choosing points again.')
			return

		def gauss(x, sigma, a, mu, b):
			g = np.exp(-0.5*((x-mu)/sigma)**2)
			return a*g + b
		
		p0 = (5.0, np.max(y)-np.min(y), np.mean([lambda1,lambda2]), np.min(y))
		fit = curve_fit(gauss, x, y, p0=p0)[0]
		
		s, a, mu, b = fit
		siga_squared = [(self.ccf[self.chosen_row][1][n+int(mu)] - self.ccf[self.chosen_row][1][int(mu)-n])**2. \
						for n in range(int(len(self.ccf[self.chosen_row][1])-abs(mu)-1))]
	
		r = abs(gauss(mu, *fit)/np.sqrt(np.mean(siga_squared)))
		w = 2.355*abs(s)
		sigma = (3./8.)*w/(1.+r)
				
		self.telluric_data[1,self.chosen_row] = -mu
		self.telluric_data[2,self.chosen_row] = sigma
		self.telluric_data[3,self.chosen_row] = 1.0
		
		fitx = np.linspace(lambda1,lambda2,50)
		self.ccf[self.chosen_row][2] = fitx
		self.ccf[self.chosen_row][3] = gauss(fitx, *fit)
		
		self.plot_ccf(self.chosen_row)
		self.rv_average(sig_clip=0)
		# Update value in table
		# TODO: Make this update, not rebuild the table from scratch!
		#self.rv_table()
		
		return

# =========================================================================
# This tab applies the Doppler shift to the science spectrum, and saves a 
# new file

class Doppler:
	def __init__(self, root_window, window):
		self.window = window
		self.root_window = root_window
		self.ap_num = None
		self.science = None
		self.template = None
		self.shifted = None
		self.vhelio = None
		
		self.fig = Figure()
		self.ax = self.fig.add_subplot(111)
		self.spec_plot = FigureCanvasTkAgg(self.fig, self.window)
		self.spec_plot.get_tk_widget().grid(column=0, rowspan=6)

		toolbar_frame = tk.Frame(self.window)
		toolbar_frame.grid(column=0,rowspan=6,sticky="w")
		self.toolbar = NavigationToolbar2Tk(self.spec_plot, toolbar_frame)		

		self.button = tk.Button(self.window, text="Shift to rest", 
				command=self.shift_spectrum)
		self.button.grid(column=1, row=0, columnspan=2)

		self.aplabel = tk.Label(self.window, text="Aperture")
		self.aplabel.grid(column=1, row=1)
		self.ap_entry = tk.Entry(self.window)
		self.ap_entry.grid(column=1, row=1)
		self.ap_entry.insert(0, 0)
		self.ap_num = int(self.ap_entry.get())
		
		self.quit = tk.Button(self.window, text="Quit", command=quit)
		self.quit.grid(column=1, row=3)
		
		self.window.bind('<Right>', self.right)
		self.window.bind('<Left>', self.left)

	
	def plot_ap(self, spectrum, ap_num, **kwargs):
		if spectrum != None:
			try:
				self.ax.plot(spectrum.wavelength[ap_num], spectrum.data[ap_num], **kwargs)
			except:
				print("Aperture %s does not exist!" %ap_num)
		return

	def shift_spectrum(self):
		from astropy.io import fits
		c = 299792458e-3 # km/s
		data, hdr = fits.getdata(self.science.file_name, header=True)
		# hdu is a numpy array
		hdul = fits.HDUList()
		new_header = self.science.header
		# Check for empty wavelength headers:
		to_remove = []
		for key in new_header:
			if "WAT2" in key: to_remove.append(key)

		for key in to_remove: del new_header[key]
		
		try:
			doppler_velocity = float(self.science.header["DOPCOR"])/c
		except:
			doppler_velocity = float(self.science.header["DOPCOR"].split()[0])/c
		
		# see https://iraf.net/irafhelp.php?val=dopcor&help=Help+Page
		one_plus_z_old = self.science.dopcor
		one_plus_z_add = np.sqrt((1.0+doppler_velocity)/(1.0-doppler_velocity))
		one_plus_z_new = one_plus_z_old*one_plus_z_add
		doppler_velocity = (one_plus_z_new**2 - 1)/(one_plus_z_new**2 + 1.0)
		
		dispstring = 'wtype=multispec'
		#dispstring = 'wtype=linear'
		#for l,line in enumerate(self.science.dispersion,int(self.science.first_beam)):
		for l,line in enumerate(self.science.dispersion,1):
			dispstring += ' spec%s = "' %l
			linedata = line.split()
			dispstring += " ".join(linedata[0:6]+['{:.13e}'.format(doppler_velocity)]+\
								linedata[7:])
			dispstring += '"'
		
		wat_counter = 1
		char_counter = 0
		while wat_counter < (len(dispstring)/68)+1:
			new_header["WAT2_{:03.0f}".format(wat_counter)] = dispstring[char_counter:char_counter+68]
			char_counter += 68
			wat_counter += 1
		
		new_header["WAT2_{:03.0f}".format(wat_counter)] = dispstring[char_counter:]
		#new_header["WAT3_001"] = 'wtype=linear'
		#new_header["CTYPE3"] = 'LINEAR'
		new_header["DOPCOR"] = '{:.2f} all'.format(doppler_velocity*c)
		
		newfile = self.science.file_name.replace(".fits", "")+"_rv.fits"
		#newfile = "temp_rv.fits"
		fits.writeto(newfile, data, new_header, overwrite=True)

		#hdul.close()
		self.shifted = spectrum.FITS(newfile)
		
		self.plot_ap(self.science, self.science.first_beam, color="C0", alpha=0.5, lw=1, zorder=3)
		self.plot_ap(self.shifted, self.shifted.first_beam, color="C0", lw=1, zorder=4)
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.shifted.first_beam),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		self.ap_num = self.shifted.first_beam
	

	def save_shifted_spectrum(self):
		# Create save file button too
		newfile = self.science.file_name.replace(".fits", "")+"_rv.fits"
		answer = False
		answer = messagebox.askyesnocancel(title="Save as...", \
					message="Save file as %s (Yes), enter a new name (No), or Cancel?" %newfile)
		
		if answer:
			fits.writeto(newfile, data, new_header, overwrite=True)
		elif not answer:
			newfile = asksaveasfilename(filetypes=("Fits files", "*.fits"))
			# Try this is if the above doesn't work..
			#import tkSimpleDialog as simpledialog
			#newfile = simpledialog.askfloat("", "Enter the heliocentric velocity:",
			#			   parent=self.root_window)
		else:
			return
		
		fits.writeto(newfile, data, new_header, overwrite=True)
		
		print("Saved file to %s" %newfile)

	
	def right(self,event):
		if self.ap_num == self.shifted.first_beam + self.shifted.apertures:
			return
		
		self.ap_num += 1
		#self.ap_text.set(self.ap_num)
		self.ax.cla()
		self.plot_ap(self.shifted, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.science, self.ap_num, color="C0", alpha=0.5, lw=1, zorder=3)
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		
	def left(self,event):
		if self.ap_num == self.shifted.first_beam:
			return
		
		self.ap_num -= 1
		#self.ap_text.set(self.ap_num)
		self.ax.cla()
		self.plot_ap(self.shifted, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.science, self.ap_num, color="C0", alpha=0.5, lw=1, zorder=3)
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		
		
		
		
		
