import Tkinter as tk
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
		self.tab_parent.add(tab1, text="Spectra")
		self.tab_parent.add(tab2, text="Cross-Correlate")
		self.tab_parent.add(tab3, text="Doppler Shift")

		self.start = Start(root_window, tab1)
		self.corrections = Corrections(root_window, tab2)
		self.doppler = Doppler(root_window, tab3)

		menubar_fmt = "{:<30}{:>2}"
		menubar = tk.Menu(root_window)
		filemenu = tk.Menu(menubar, tearoff=0)
		filemenu.add_command(label=menubar_fmt.format("New spectrum", u"\u2318".encode('utf-8')+"S"), command=self.open_science)
		root_window.bind('<Command-s>', self.open_science)
		filemenu.add_command(label=menubar_fmt.format("New template", u"\u2318".encode('utf-8')+"T"), command=self.open_template)
		root_window.bind('<Command-t>', self.open_template)
		filemenu.add_command(label=menubar_fmt.format("Export shifted spectrum", u"\u2318".encode('utf-8')+"E"), command=self.save_shifted_spectrum)
		root_window.bind('<Command-e>', self.save_shifted_spectrum)
		filemenu.add_separator()
		filemenu.add_command(label=menubar_fmt.format("Exit",""), command=root_window.quit)
		menubar.add_cascade(label="File", menu=filemenu)
		
		exportmenu = tk.Menu(menubar, tearoff=1)
		exportmenu.add_command(label=menubar_fmt.format("Export RV", u"\u2318".encode('utf-8')+"E"), command=self.write_log)
		root_window.bind('<Command-e>', self.write_log)
		menubar.add_cascade(label="Export", menu=exportmenu)

		root_window.config(menu=menubar)
		
	def open_science(self, event=None):
		filename = self.start.open_science_file()
		#self.start.science = spectrum.FITS(filename)
		self.corrections.science = self.start.science
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
			log.write('Filename\tBCV\tRV\tRV_err\tAps\tAps_RV\tAps_err\n')
		
		keep = np.where(self.corrections.rv_data[3]==1.)[0]
		
		log.write(self.corrections.science.file_name.split('/')[-1]+'\t')
		log.write('{:>+.3f}'.format(self.corrections.bcv)+'\t')
		log.write('{:>+.3f}'.format(self.corrections.vhelio[0])+'\t')
		log.write('{:>.3f}'.format(self.corrections.vhelio[1])+'\t')
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
		
		self.fig = Figure()
		self.ax = self.fig.add_subplot(111)
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
		self.plot_ap(self.science, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.template, self.ap_num, lw=1, zorder=2, color="gray")
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		return
	

	def plot_ap(self, spectrum, ap_num, **kwargs):
		if spectrum != None:
			try:
				self.ax.plot(spectrum.wavelength[ap_num], spectrum.data[ap_num], **kwargs)
			except:
				print "Aperture %s does not exist!" %ap_num
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
			self.template_plot = self.plot_ap(self.template, self.ap_num, lw=1, color="gray")

		else:
			self.template_plot = self.plot_ap(self.template, self.science.first_beam, lw=1, color="gray")

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
				import tkSimpleDialog as simpledialog
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
		
		print "Using VHELIO =", "{:.3f}".format(vh), "for template spectrum."	
		return
				

	def right(self,event):
		if self.ap_num == self.science.first_beam + self.science.apertures:
			return
		
		self.ap_num += 1
		self.ap_text.set(self.ap_num)
		self.ax.cla()
		self.plot_ap(self.science, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.template, self.ap_num, lw=1, zorder=2, color="gray")
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()
		
	def left(self,event):
		if self.ap_num == self.science.first_beam:
			return
		
		self.ap_num -= 1
		self.ap_text.set(self.ap_num)
		self.ax.cla()
		self.plot_ap(self.science, self.ap_num, lw=1, zorder=3)
		self.plot_ap(self.template, self.ap_num, lw=1, zorder=2, color="gray")
		self.ax.text(0.05, 0.95, "Aperture %s" %int(self.ap_num),\
					va="top", transform=self.ax.transAxes)
		self.spec_plot.draw()

# =========================================================================
#rv_data = rv.rv_by_aperture(template, science, range(science.apertures))
#new_data_list, avg, std, n = rv.rv_average(rv_data, 0.2)


class Corrections:
	def __init__(self, root_window, window):
		self.window = window
		self.root_window = root_window
		self.science = None
		self.template = None
		self.ccf = []
		
		self.fig = Figure()
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
		
		self.cbutton = tk.Button(self.window, text="Run cross-correlation", 
				command=self.calculate_rvs)
		self.cbutton.grid(column=1, row=0, columnspan=2)
		
		self.sglabel = tk.Label(self.window, text="Outlier sigma:")
		self.sglabel.grid(column=1, row=1, sticky="e")
		self.sg_text = tk.StringVar()
		self.sg_text.set(0)
		self.sig_clip = float(self.sg_text.get())
		
		self.sg_entry = tk.Entry(self.window, textvariable=self.sg_text, width=8)
		self.sg_entry.grid(column=2, row=1, sticky="w")
		self.new_sg
		self.sg_entry.bind("<Return>", self.new_sg) 
				
		self.rv_data = None
		self.keep = None
		self.new_data_list = None
		self.vhelio = None
		self.vm_template = None
		
		self.vh = tk.StringVar()
		self.vhlabel = tk.Label(self.window, text="RV helio:")
		self.vhlabel.grid(column=1, row=2, sticky="e")
		self.vhbox = tk.Label(self.window, textvariable=self.vh)
		self.vhbox.grid(column=2, row=2, sticky="w")
			
		self.bcv_text = tk.StringVar()
		self.bcvlabel = tk.Label(self.window, text="BCV:")
		self.bcvlabel.grid(column=1, row=3, sticky="e")
		self.bcvbox = tk.Label(self.window, textvariable=self.bcv_text)
		self.bcvbox.grid(column=2, row=3, sticky="w")
		
		self.fig.canvas.mpl_connect('pick_event', self.on_click)
		self.fig.canvas.mpl_connect('key_press_event', self.figure_key_press)
		
		self.rv_table_box = tk.Frame(self.window)
		self.rv_table_box.grid(row=4, column=1, columnspan=2, sticky="nw")
		self.rv_checks = None
		
	
	def bcv_corr(self, spectrum):
		if False:
		#if "BCV" in spectrum.header.keys():
			return spectrum.header["BCV"]
		
		else:
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
			'''
			try:
				UT_MID = spectrum.header["UT-DATE"]+" "+spectrum.header["UT-START"]
			except TypeError:
				UT_MID = spectrum.header["UT-DATE"]+" "+spectrum.header["UT-TIME"]
			'''
			UT_MID = spectrum.header['DATE']
			barycorr = sc.radial_velocity_correction(obstime=Time(UT_MID, format='isot', scale="utc"), location=site) 
			bcv = barycorr.to(u.km/u.s).value
			return bcv

	
	def calculate_rvs(self):
		self.rv_data, self.ccf, self.comparison = rv.rv_by_aperture(\
										self.template, self.science, \
										self.science.data.keys())
		
		self.correct_rvs()
		vm = self.vm_template + self.rv_data[1]
		ckms = c.to(u.km/u.s).value
		vmvb_c = vm*self.bcv/ckms
		self.rv_data[1] = vm + self.bcv + vmvb_c
		
		self.rv_table()
		self.rv_average(sig_clip=0)
		
	
	def rv_average(self, sig_clip=None):
		if sig_clip == None:
			sig_clip = self.sig_clip
		rv_data, avg, std, n = rv.rv_average(self.rv_data, sig_clip)
		
		#import pdb
		#pdb.set_trace()
		for i,r in enumerate(rv_data.T):
			aploc = np.where(self.rv_data[0]==r[0])[0][0]
			if r[3] == 1.0:
				self.rv_data[3,aploc] = 1.0
			else:
				self.rv_data[3,aploc] = 0.0
		
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
			'''
			if x_move in self.rv_data[0,self.keep]:
				loc = np.where(self.rv_data[0,self.keep]==x_move)[0][0]
				self.rv_data[3,loc] = abs(self.rv_data[3,loc]%2-1)
				self.rv_data = self.rv_data[:,self.rv_data[0]!=x_move]
				self.keep = np.where(self.rv_data[3]==1)[0]
			
			# NEED TO FIX THESE NEW_DATA_LISTs 
			else:
				loc = np.where(self.rv_data[0]==x_move)[0][0]
				self.new_data_list = np.append(self.new_data_list, np.array([self.rv_data[:,loc]]).T, axis=-1)
				self.new_data_list = self.new_data_list[:,self.new_data_list.argsort()[0]]
			
			self.new_data_list, avg, std, n = rv.rv_average(self.new_data_list, 2.0)
			'''
			loc = np.where(self.rv_data[0]==x_move)[0][0]
			#self.select_aps(loc)
			# Used to be that a click unselects. Now, let's plot the CCF
			self.plot_ccf(loc)
			
			#self.rv_data[3,loc] = abs(self.rv_data[3,loc]%2-1)
			#self.CheckVar.set(abs(self.rv_data[3,loc]%2-1))
			'''
			self.keep = np.where(self.rv_data[3]==1)[0]
			avg = np.mean(self.rv_data[1,self.keep])
			std = np.sqrt(np.sum([(d-avg)**2. for d in \
						 self.rv_data[1,self.keep]])/len(self.rv_data[1,self.keep]))
			std = std/np.sqrt(len(self.rv_data[1,self.keep])-1)
			n = len(self.keep)
			
			self.vhelio = [avg,std]
			self.doppler_velocity = self.vhelio[0] - self.bcv
			self.ax.cla()
			self.ax.text(0.05, 0.95, "%s apertures used" %int(n),\
						va="top", transform=self.ax.transAxes)
			
			self.plot_results()
			self.vh.set("{:.2f}".format(self.vhelio[0])+u" \u00B1 "+"{:.2f}".format(self.vhelio[1]))
			'''
			# Call select_aps too
			#self.select_aps
			
			
	def rv_table(self):
		rv_frame = tk.Canvas(self.rv_table_box, borderwidth=2)
		#rv_frame.pack(side="left", fill="both", padx=10, pady=10)
		rv_frame.grid(row=0, column=0, sticky="news")
		buttons_frame = tk.Frame(rv_frame, borderwidth=2, relief="groove")
		buttons_frame.grid(row=1, column=0, sticky="news")

		ROWS = len(self.rv_data[0])
		ROWS_DISP = 15
		
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
		
		rv_scroll = ttk.Scrollbar(self.rv_table_box, orient="vertical", command=rv_frame.yview)
		#rv_scroll.pack(side="right", fill="y")
		rv_scroll.grid(row=0, column=0, sticky="nes")
		rv_frame.configure(yscrollcommand=rv_scroll.set)#, scrollregion=rv_frame.bbox("all"))
		
		rv_frame.create_window((0,0), window=buttons_frame, anchor=tk.NW)
		buttons_frame.update_idletasks()
		#bbox = self.rv_table_box.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.
		bbox = rv_frame.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.

		# Define the scrollable region as entire canvas with only the desired
		# number of rows and columns displayed.
		w, h = bbox[2]-bbox[1], bbox[3]-bbox[1]
		dw, dh = int((w/4) * 4), int((h/ROWS) * ROWS_DISP)
		rv_frame.configure(scrollregion=bbox, width=w, height=dh)
		#buttons_frame.configure(width=2*w, height=dh)
		#self.rv_table_bow.configure(width=2*w, height=dh)
		
		'''
		rows = 10
		for i in range(1,rows):
			for j in range(1,6):
				self.rv_table_box = tk.Label(self.window)
				self.rv_table_box.grid(column=2, row=4, stick="news")

				button = Button(rv_frame, padx=7, pady=7, text="[%d,%d]" % (i,j))
				button.grid(row=i, column=j, sticky='news')

		vsbar = Scrollbar(FMas, orient="vertical", command=self.rv_table_box.yview)
		vsbar.grid(row=3, column=1)
		'''
	
	def select_aps(self, row):
		self.rv_checks[row].set(int(abs(self.rv_data[3,row]%2-1)))
		self.rv_data[3,row] = abs(self.rv_data[3,row]%2-1)
		
		self.keep = np.where(self.rv_data[3]==1)[0]
		self.rv_average(sig_clip=0)
		
		'''
		#self.new_data_list = self.rv_data[:,np.where(self.rv_data[3]==1)[0]]
		avg = np.mean(self.rv_data[1,self.keep])
		n = len(self.keep)
		std = np.sqrt(np.sum([(d-avg)**2. for d in \
					 self.rv_data[1,self.keep]])/float(n))
		std = std/np.sqrt(float(n-1))
		
		std1 = np.sqrt(sum([(err)**2. for err in \
					   self.rv_data[2,self.keep]]))/float(n)
		
		total_std = np.sqrt(std**2. + std1**2.)
		
		self.vhelio = [avg,total_std]
		'''
		
		#self.plot_results()
		'''
		self.vh.set("{:.2f}".format(self.vhelio[0])+u" \u00B1 "+"{:.2f}".format(self.vhelio[1]))
		
		ckms = c.to(u.km/u.s).value
		self.science.header["DOPCOR"] = (self.vhelio[0] - self.bcv)/(1.+self.bcv/ckms)
		self.science.header["VHELIO"] = self.vhelio[0]
		
		weights = np.array([1.0/w**2. for w in self.rv_data[2,self.keep]])
		print self.doppler_velocity[self.keep]
		avg = np.average(self.doppler_velocity[self.keep], 
						 weights=weights)
	
		std1 = np.sqrt(np.sum(weights)/float(len(weights)))
		print self.science.header["DOPCOR"], avg, std1
		'''

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
		doppler_shift = self.rv_data[1][row] - self.vm_template*vbc1 - self.bcv
		self.ccf_ax.axvline(-doppler_shift/vbc1,\
							color='C0', lw=0.8, ls='-')

		self.rv_plot.draw()

	
	def figure_key_press(self, event):
		if event.key in "gGkK":
			x_select = event.xdata #event.artist.get_xdata()[ind]
			if self.gaussian_point_chosen == None:
				self.gaussian_point_chosen = x_select
				print('Point 1 selected for refitting. Please choose Point 2.')
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
		x1 = np.where(self.ccf[self.chosen_row][0]<=lambda1)[0][-1]
		x2 = np.where(self.ccf[self.chosen_row][0]>=lambda2)[0][0]
		
		def gauss(x, sigma, a, mu):
			g = np.exp(-0.5*((x-mu)/sigma)**2)
			return a*g

		x = self.ccf[self.chosen_row][0][x1:x2]
		y = self.ccf[self.chosen_row][1][x1:x2]
		
		fit = curve_fit(gauss, x, y,
						p0=(5., np.max(y), np.mean([lambda1,lambda2]))
						)[0]
		
		s, a, mu = fit
		siga_squared = [(self.ccf[self.chosen_row][1][n+mu] - self.ccf[self.chosen_row][1][mu-n])**2. \
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
		
		fitx = np.linspace(lambda1,lambda2,50)
		self.ccf[self.chosen_row][2] = fitx
		self.ccf[self.chosen_row][3] = gauss(fitx, *fit)
		
		self.plot_ccf(self.chosen_row)
		self.rv_average(sig_clip=0)
		# Update value in table
		self.rv_table()
		
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
				print "Aperture %s does not exist!" %ap_num
		return

	def shift_spectrum(self):
		from astropy.io import fits
		c = 299792458e-3 # km/s
		data, hdr = fits.getdata(self.science.file_name, header=True)
		# hdu is a numpy array
		hdul = fits.HDUList()
		new_header = self.science.header
		
		try:
			doppler_velocity = float(self.science.header["DOPCOR"])
		except:
			doppler_velocity = float(self.science.header["DOPCOR"].split()[0])
		
		dispstring = 'wtype=multispec'
		#for l,line in enumerate(self.science.dispersion,int(self.science.first_beam)):
		for l,line in enumerate(self.science.dispersion,1):
			dispstring += ' spec%s = "' %l
			linedata = line.split()
			dispstring += " ".join(linedata[0:6]+['{:.13e}'.format(doppler_velocity/c)]+\
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
		new_header["DOPCOR"] = '{:.2f} all'.format(doppler_velocity)
		
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
		answer = messagebox.askyesnocancel(title="Save as...", \
					message="Save file as %s (Yes), enter a new name (No), or Cancel?" %newfile)
		
		if answer == True:
			fits.writeto(newfile, data, new_header, overwrite=True)
		elif answer == False:
			newfile = asksaveasfilename(filetypes=("Fits files", "*.fits"))
			# Try this is if the above doesn't work..
			#import tkSimpleDialog as simpledialog
			#newfile = simpledialog.askfloat("", "Enter the heliocentric velocity:",
			#			   parent=self.root_window)
		else:
			return
		
		fits.writeto(newfile, data, new_header, overwrite=True)
		
		print "Saved file to %s" %newfile

	
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
		
		
		
		
		
