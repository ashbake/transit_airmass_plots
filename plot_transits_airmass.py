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
planet = 'WASP-76b'
filename= datapath + 'transits_%s.txt' %planet

data = np.loadtxt(filename,dtype='str',delimiter=',')

startimes = data[:,12]
endtimes  = data[:,13]
midtimes  = data[:,9]
duration  = data[:,7].astype(float)
fracs     = data[:,14].astype(float)

full_transits =  np.where(fracs==1.0)[0]


target = FixedTarget.from_name(planet.strip('b'))
apo    = Observer.at_site('KPNO')
observer = apo

fig, axs = plt.subplots(4,4,figsize=(10,10),sharey=True)
j = 0
for m in range(4):
	for n in range(3):
		#start_time = Time('2020-' + startimes[full_transits[j]].replace('/','-').replace('-2020','') + ':00')
		#end_time = Time('2020-' + endtimes[full_transits[j]].replace('/','-').replace('-2020','') + ':00')
		# calc transit
		half_duration = datetime.timedelta(hours=0.5*duration[0])
		if midtimes[full_transits[j]].startswith('01'):
			mid_time = Time('2021-' + midtimes[full_transits[j]].replace('/','-').replace('-2021','') + ':00')
		else:
			mid_time = Time('2020-' + midtimes[full_transits[j]].replace('/','-').replace('-2020','') + ':00')
		start_time = -1*half_duration + mid_time.datetime
		end_time   = half_duration + mid_time.datetime

		ax = plot_airmass(target, apo, start_time +datetime.timedelta(hours=2) , ax=axs[m,n],brightness_shading=False, altitude_yaxis=False)

		# get moon info
		moon_ill = round(moon.moon_illumination(mid_time),2)
		moon_coord = moon.get_moon(mid_time,location=apo.location) 
		moon_sep = moon_coord.ra  - target.ra

		#plot transit
		ax.axvspan(start_time, end_time,ymin=0, ymax=10, color='lavender')
		
		# narrower x-axis
		xlo = datetime.datetime(day=start_time.day, year=start_time.year, month= start_time.month)
		axs[m,n].set_xlim(xlo - datetime.timedelta(hours=0),xlo + datetime.timedelta(hours=14))

		txt_shift = 1 if start_time.hour > 4 else -8

		txt = axs[m,n].annotate('moon:%s\n%sdeg'%(moon_ill,round(moon_sep.value)), \
			(xlo - datetime.timedelta(hours=txt_shift),2.5),fontsize=8)
		txt.set_bbox(dict(facecolor='white', alpha=0.5))

		#start = start_time - datetime.timedelta(hours=10)
		#twilight1 =observer.twilight_evening_astronomical(Time(start), which='next').datetime
		#twilight2 = observer.twilight_morning_astronomical(Time(start), which='next').datetime
		#ax.axvspan(twilight1, twilight2, ymin=0, ymax=1, color='lightgrey',zorder=-10)

		j+=1

plt.subplots_adjust(hspace=0.7,wspace=0.2)







planet = 'WASP-127b'
filename= datapath + 'transits_%s.txt' %planet

data = np.loadtxt(filename,dtype='str',delimiter=',')

startimes = data[:,13]
endtimes  = data[:,14]
midtimes  = data[:,10]
duration  = data[:,7].astype(float)
fracs     = data[:,15].astype(float)

full_transits =  np.where(fracs>0.6)[0]


target = FixedTarget.from_name(planet.strip('b'))
apo    = Observer.at_site('KPNO')


n=3# fill in 4th row
for m in range(4):
	j= m
	half_duration = datetime.timedelta(hours=0.5*duration[0])
	if midtimes[full_transits[j]].startswith('01'):
		mid_time = Time('2021-' + midtimes[full_transits[j]].replace('/','-').replace('-2021','') + ':00')
	else:
		mid_time = Time('2020-' + midtimes[full_transits[j]].replace('/','-').replace('-2020','') + ':00')
	start_time = -1*half_duration + mid_time.datetime
	end_time   = half_duration + mid_time.datetime

	ax = plot_airmass(target, apo, start_time+datetime.timedelta(hours=2), ax=axs[m,n],brightness_shading=False, altitude_yaxis=False)

	# get moon info
	moon_ill = round(moon.moon_illumination(mid_time),2)
	moon_coord = moon.get_moon(mid_time,location=apo.location) 
	moon_sep = np.abs(moon_coord.ra  - target.ra)

	#plot transit
	ax.axvspan(start_time, end_time,ymin=0, ymax=10, color='plum',alpha=0.8)
	
	# narrower x-axis
	xlo = datetime.datetime(day=start_time.day, year=start_time.year, month= start_time.month)
	axs[m,n].set_xlim(xlo - datetime.timedelta(hours=2),xlo + datetime.timedelta(hours=14))

	txt_shift = 1 if start_time.hour > 4 else -8

	txt = axs[m,n].annotate('moon:%s\n%sdeg'%(moon_ill,round(moon_sep.value)), \
		(xlo - datetime.timedelta(hours=txt_shift),2.5),fontsize=8)
	txt.set_bbox(dict(facecolor='white', alpha=0.5))

	start = start_time - datetime.timedelta(hours=10)
	twilight1 =observer.twilight_evening_astronomical(Time(start), which='next').datetime
	twilight2 = observer.twilight_morning_astronomical(Time(start), which='next').datetime
	ax.axvspan(twilight1, twilight2, ymin=0, ymax=1, color='lightgrey',zorder=-10)

	j+=1


axs[0,3].set_title('WASP-127b',color='plum',fontsize=12)

axs[0,1].set_title("WASP-76b",color='steelblue',fontsize=12)


plt.savefig('observing_charts.eps')
plt.savefig('observing_charts.pdf')

###################### Generate data for 


