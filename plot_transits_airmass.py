##############################################################
# 1. download transit data from exoplanet archive 
# 2. put file in data/planet/ folder
# 3. update planets list (to do, make it an input)
# 4. run scripts
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

	return startimes[full_transits], endtimes[full_transits], midtimes[full_transits], duration[full_transits], target, obs


def plot_all(startimes, endtimes, midtimes, duration, target, obs):
	n = len(startimes)
	w = int(np.ceil(np.sqrt(n)))
	h = int(np.ceil(n/w))

	half_duration = datetime.timedelta(hours=0.5*duration[0])

	fig, axs = plt.subplots(w,h,figsize=(3*h,3*w),sharey=True)
	
	m = 0
	for i in np.arange(w):
		for j in np.arange(h):
			if m>=n:
				axs[i,j].set_visible(False)
				continue

			year = midtimes[m][6:10]
			mid_time = Time(year  + '-' + midtimes[m].replace('/','-').replace('-%s'%year,'') + ':00')
			
			start_time = -1*half_duration + mid_time.datetime
			end_time   = half_duration + mid_time.datetime

			# add airmass plot - sample airmass time window from UT 0 so x axis reads correct data - created start of grid by dumb hack of taking start_time and subtracting start_time.hour
			dt_array = np.arange(0,30,0.1)
			time_array = [start_time - datetime.timedelta(hours=start_time.hour) + datetime.timedelta(hours=q) for q in dt_array]
			ax = plot_airmass(target, obs, time_array, ax=axs[i,j],brightness_shading=False, altitude_yaxis=False)
			
			# get moon info
			moon_ill = round(moon.moon_illumination(mid_time),2)
			moon_coord = moon.get_moon(mid_time,location=obs.location) 
			moon_sep = np.abs(moon_coord.ra  - target.ra)

			axs[i,j].axvspan(start_time, end_time,ymin=0, ymax=10, color='plum',alpha=0.8)

			# narrower x-axis
			xlo = datetime.datetime(day=start_time.day, year=start_time.year, month= start_time.month)
			axs[i,j].set_xlim(xlo + datetime.timedelta(hours=4),xlo + datetime.timedelta(hours=17))

			# Put moon info in each panel
			txt_shift = 5 if start_time.hour > 8 else 12
			txt = axs[i,j].annotate('moon:%s\n%sdeg'%(moon_ill,round(moon_sep.value)),\
				(xlo + datetime.timedelta(hours=txt_shift),2.5),\
				fontsize=8)
			txt.set_bbox(dict(facecolor='white', alpha=0.5))

			start = start_time - datetime.timedelta(hours=10)
			twilight1 =obs.twilight_evening_astronomical(Time(start), which='next').datetime
			twilight2 = obs.twilight_morning_astronomical(Time(start), which='next').datetime
			ax.axvspan(twilight1, twilight2, ymin=0, ymax=1, color='lightgrey',zorder=-10)

			fig.suptitle(planet + ' from ' + site, fontsize=18)
			m+=1

	plt.subplots_adjust(hspace=0.6,top=0.92)


def run():
	# goal: run from command line: plot_transits airmass -site Keck -target WASP-76b -color coral -return_coord -ext pdf
	pass

if __name__=='__main__':
	planets = ['WASP-12b', 'WASP-127b']#,'WASP-31b','WASP-103b','WASP-80b','WASP-76b', 'WASP-127b']
	site = 'Keck' #'Palomar', 'KPNO'

	for planet in planets:
		filename = datapath + 'transits_%s.txt' %planet

		args = read_transits(filename,site=site)

		plot_all(*args)

		plt.savefig('./plots/observing_chart_%s_%s.eps'%(planet,site)
		plt.savefig('./plots/observing_charts_%s_%s.pdf'%(planet,site)

###################### 


