##############################################################
# purpose: plot best planets to do atmospheric characterization
# 1. download hot jupiter list from https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PS
# 2. run script
###############################################################

import numpy as np
import matplotlib.pyplot as plt

from astropy.time import Time

import matplotlib,sys,os
from astropy.io import fits
from astropy.table import Table

plt.ion()

font = {'weight' : 'normal',
        'size'   : 8}

matplotlib.rc('font', **font)

data_dir = './data/'

def read_targets(file,file2=None):
	"""
	read exoplanet archive confirmed planet csv file
	"""
	hj1 = Table.read(file,format='ascii')

	if file2!=None:
			#https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=transitspec
		hj2 = Table.read(file2,format='ascii')

		hj = vstack([hj1,hj2],join_type='outer')
	else:
		hj = hj1

	vmag  = hj['sy_vmag']
	jmag  = hj['sy_jmag']
	teff  = hj['st_teff'] 
	rstar = hj['st_rad']
	a     = hj['pl_orbsmax']
	rp    = hj['pl_radj']
	mass  = hj['pl_bmassj']
	name  = hj['pl_name']
	dec   = hj['dec']
	
	return vmag, jmag, teff, rstar, a, rp, mass, name, dec

def get_scale_height(teff, rstar, a, rp, mass, names):
	"""
	calculate planet scale heights
	"""

	kb = 1.38064852e-23 #m2 kg s-2 K-1
	mu = 2.3 * 1.6605390e-27 # kg
	G = 6.67408e-11 #m3 kg-1 s-2

	# estimate Temp planet in kelvin
	Teq = (1/4.)**(1/4.) * teff* np.sqrt(0.00465047 * rstar/a) # 0.00465047 AU per Rsun
	gravity = G * mass * 1.898e27 / (rp * 69911000.0)**2 #kg from m_jup

	# get H -> kb * T/ mu / g
	H = kb * Teq / mu / gravity # meters

	# calculate A like Sing did
	A = 2* rp*0.10049 * H*1.4374e-9 / rstar**2 #0.10049 rsun/rjupiter, 1.4374e-9 rsun/meters

	return A

def calc_s(mode, tot_eff=0.2,ntransits=3,target_snr=3,tel_diam=5.1,dl_l=0.001):
	"""
	calc  snr 3sigma

	mode: NIR or Optical (J band or V band, respectively)

	"""
	tel_area = np.pi * (tel_diam/2.)**2
	exp_one_transit = 2.3 * 3600. # approx half of a hot jupiter transit
	mags = np.arange(5,17,0.2)
	s_3sig=  np.zeros(len(mags))
	# optical
	if mode=='V':
		base_jansky = 3640
	elif mode =='J':
		base_jansky = 1600
	for i in range(len(mags)):
		signal = tot_eff * 10**(-0.4*mags[i]) * base_jansky*1.51e7 * tel_area * exp_one_transit *dl_l # photons in band
		noise = np.sqrt(2/signal) # sqrt 2 bc combining in and out of transit
		s_3sig[i] = target_snr*noise/np.sqrt(ntransits)


	return mags,s_3sig

def plot_planets(mags, tele, A, mag, names, dec, scale_fac=2, decmin=-20):
	"""
	plot planets: V mag vs A (transit %) 
	"""
	fig, ax = plt.subplots(figsize=(10,8))

	isub = np.where(dec > decmin)[0]
	plt.semilogx(scale_fac*A[isub]*100,mag[isub],'ko',ms=3)
	plt.ylim(15,7.2)
	plt.xlim(8e-4,0.45)
	plt.plot(tele*100,mags,'r--',lw=.6)

	for i in range(len(isub)): # labelplanets
		if np.isnan(A[isub][i]):
			pass
		else:
			if (type(mag[isub][i]) != np.ma.core.MaskedConstant) & (type(A[isub][i]) != np.ma.core.MaskedConstant):
				if (mag[isub][i] < 15) & (2.51*A[isub][i]*100 > 8e-4):
					plt.text(scale_fac*A[isub][i]*100,mag[isub][i],names[isub][i],fontsize=4,color='steelblue')

	plt.grid(color='lightgray',linestyle='--',zorder=-1000)
	ax.tick_params(axis='both',direction='in',width=1.3,which='both',length=3,pad=5)
	plt.xlabel('Transit Signal (%)',fontweight='bold',fontstyle='italic')
	plt.ylabel('%s magnitude'%mode,fontweight='bold',fontstyle='italic')

	# shade for hubble
	plt.text(0.0012, 7.8, '%s Times Scale Height'%scale_fac, fontsize=14,color='gray')

	plt.title('SNR=3, $N_{transits}=%s$'%ntransits,color='red')
	plt.savefig('plots/Sing_inspired_%s.png'%scale_fac)


if __name__=='__main__':
	file = data_dir + 'hotjupiters.csv' # period <10days, radius > 0.8 rjup
	vmag, jmag, teff, rstar, a, rp, mass, names, dec = read_targets(file, file2=None)
	A = get_scale_height(teff, rstar, a, rp, mass, names)

	# choose either 'J' or 'V' band
	mode = 'J'
	mag  = jmag if mode=='J' else vmag

	# define instrument/observatory params
	binning = 20 # chose these to match instr. i care about, not J/V bands
	res = 90000.0      # resolution of instrument
	ntransits = 3      # no. of transits assuming to have observed
	dl_l = binning/res # cal delta lambda over lambda
	tel_diam = 4.5     # telescope diameter in meters

	mags, tele = calc_s(mode=mode, tot_eff=0.1,  ntransits=ntransits, target_snr=3, tel_diam=tel_diam, dl_l=dl_l)

	plot_planets(mags, tele, A, mag, names, dec, scale_fac=3, decmin=-20)


###################### 

# goal: run from command line: input tel diameter, add NIR option using J band mag, DEC filter
