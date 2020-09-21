##############################################################
# Integrate PBay Spectrum with actualy Alluxa profile to get Measurements
# To make spectrum, go to Kavli research folder where PBay is stored
# Outputs: Plots saved to plots folder
# Inputs: Spectra from pbay + their log should be saved to spectra folder
###############################################################

import numpy as np

import matplotlib.pyplot as plt
from astropy.time import Time
from astroplan import FixedTarget, Observer
from astroplan.plots import plot_airmass
from astroplan import moon

import matplotlib,sys,os
from astropy.io import fits
from scipy.signal import gaussian
import ephem
from astropy.table import Table
import datetime

plt.ion()

font = {'weight' : 'normal',
        'size'   : 8}

matplotlib.rc('font', **font)

datapath = './data/'

def read_transits(file, site='Keck', frac=1.0):
	"""
	read transits from exoplanet archive transit planner CSV file

	site options from astroplan: e.g. KPNO, Keck
	"""


	data = Table.read(file,format='ascii')

	startimes = data['targetobsstartcalendar']
	endtimes  = data['targetobsendcalendar']
	midtimes  = data['midpointcalendar']
	duration  = data['transitduration']
	fracs     = data['fractionobservable']

	full_transits =  np.where(fracs>=frac)[0]

	target = FixedTarget.from_name(planet.strip('b'))
	obs    = Observer.at_site(site)

	return startimes, endtimes, midtimes, duration, fracs, full_transits, target, obs




def plot_all(startimes, endtimes, midtimes, duration, fracs, full_transits, target, obs):
	n = len(full_transits)
	w = int(round(np.sqrt(n)))
	h = int(round(n/w))

	half_duration = datetime.timedelta(hours=0.5*duration[0])

	fig, axs = plt.subplots(w,h,figsize=(10,10),sharey=True)
	
	m = 0
	while m<n:
		for i in range(w):
			for j in range(h):
				year = midtimes[m][6:10]
				mid_time = Time(year  + '-' + midtimes[m].replace('/','-').replace('-%s'%year,'') + ':00')
				
				start_time = -1*half_duration + mid_time.datetime
				end_time   = half_duration + mid_time.datetime
	
				ax = plot_airmass(target, apo, start_time+datetime.timedelta(hours=2), ax=axs[i,j],brightness_shading=False, altitude_yaxis=False)
				
				# get moon info
				moon_ill = round(moon.moon_illumination(mid_time),2)
				moon_coord = moon.get_moon(mid_time,location=apo.location) 
				moon_sep = np.abs(moon_coord.ra  - target.ra)

				axs[i,j].axvspan(start_time, end_time,ymin=0, ymax=10, color='plum',alpha=0.8)


				# narrower x-axis
				xlo = datetime.datetime(day=start_time.day, year=start_time.year, month= start_time.month)
				axs[i,j].set_xlim(xlo + datetime.timedelta(hours=4),xlo + datetime.timedelta(hours=17))

				txt_shift = 1 if start_time.hour > 4 else -8

				txt = axs[i,j].annotate('moon:%s\n%sdeg'%(moon_ill,round(moon_sep.value)), \
					(xlo - datetime.timedelta(hours=txt_shift),2.5),fontsize=8)
				txt.set_bbox(dict(facecolor='white', alpha=0.5))

				start = start_time - datetime.timedelta(hours=10)
				twilight1 =obs.twilight_evening_astronomical(Time(start), which='next').datetime
				twilight2 = obs.twilight_morning_astronomical(Time(start), which='next').datetime
				ax.axvspan(twilight1, twilight2, ymin=0, ymax=1, color='lightgrey',zorder=-10)

				m+=1

		if i+j == w+h:
			break




if __name__=='__main__':
	planets = ['WASP-127b', 'WASP-76b']

	for planet in planets:
		filename = datapath + 'transits_%s.txt' %planet

		args = read_transits(filename,site='Keck')

		plot_all(*args)

		plt.savefig('observing_chart_%s.eps'%planet)
		plt.savefig('observing_charts_%s.pdf'%planet)

###################### Generate data for 


