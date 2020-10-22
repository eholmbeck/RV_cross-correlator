try:
	import Tkinter as tk
except:
	# Detecting python3
	import tk
	
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np
import spectrum
from tkinter.filedialog import askopenfilename

'''
template = spectrum.FITS('/Volumes/Files/Desktop/Spectra/Magellan/hd122563blue_multi_magellan_07_th.fits')
science = spectrum.FITS('/Volumes/Files/Desktop/Spectra/Magellan/201707_Magellan/cd-299441blue_multi.fits')
import rv

#for i in range(int(science.first_beam), int(science.first_beam)+10):
#	rv.cross_correlate(template, science, i)

start = int(science.first_beam)
data = rv.cross_correlate(template, science, 94)
exit()

rv_data, ccf = rv.rv_by_aperture(template, science, range(start,start+science.apertures))
print rv_data[1]

new_data_list, avg, std, n = rv.rv_average(rv_data)

plt.errorbar(rv_data[0], rv_data[1], yerr=rv_data[2], ls="", marker="o")
plt.scatter(new_data_list[0], new_data_list[1], facecolor="None", edgecolor="red", marker="o", zorder=4)
#plt.plot([0,40], [avg, avg])

print std
#plt.plot(soln)
plt.show()

exit()
'''


# To do: run in batch mode; read template to use from header

import tabs

root = tk.Tk()

# get screen width and height
ws = root.winfo_screenwidth() # width of the screen
hs = root.winfo_screenheight() # height of the screen

# set the dimensions of the screen 
# and where it is placed
root.geometry('%dx%d+%d+%d' % (ws*0.8, hs*0.8, ws*0.1, 0.0))

root.title("Radial Velocity Calculator")

session = tabs.Session(root)
session.pack()
root.mainloop()


'''
top = tk.Tk()

# Code to add widgets will go here...
frame = tk.Frame(top)

button = tk.Button(frame, text="Open Science File", command=open_file)

button.pack(side=tk.LEFT)

top.mainloop()

exit()
fig,ax = plt.subplots(1,1)
ax.plot(spec[0], spec[1])
spec_plot = FigureCanvasTkAgg(fig, top)

spec_plot.get_tk_widget().pack()

top.mainloop()
'''