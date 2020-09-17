##############################################################
# Integrate PBay Spectrum with actualy Alluxa profile to get Measurements
# To make spectrum, go to Kavli research folder where PBay is stored
# Outputs: Plots saved to plots folder
# Inputs: Spectra from pbay + their log should be saved to spectra folder
###############################################################

import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import matplotlib,sys,os
from astropy.io import fits
from scipy.signal import gaussian

#plt.ion()


font = {'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

def gaussian(x, shift, sig):
    ' Return normalized gaussian with mean shift and var = sig^2 '
    gaus = np.exp(-.5*((x - shift)/sig)**2)/(sig * np.sqrt(2*np.pi))
    return gaus/max(gaus)


def bindat(x,y,nbins):
    """
    Bin Data
    
    Inputs:
    ------
    x, y, nbins
    
    Returns:
    --------
    arrays: bins, mean, std  [nbins]
    """
    # Create bins (nbins + 1)?
    n, bins = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    
    # Calculate bin centers, mean, and std in each bin
    bins = (bins[1:] + bins[:-1])/2
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    
    # Return bin x, bin y, bin std
    return bins, mean, std

######################### INPUTS #####################
t_exp   = 900.0  # seconds
res     = 90000 # resolution l/dl
mode    = 'gaus' # pick 'gaus' for gaussians matching/scaled from HARPS or 'model' for Goyal model

######################### LOAD THINGS #####################
# Load spectrum 
def calc_spec(D1, D2, snr_1, nspec, nvis, res, NaK_rat, nbins=4500,mode='gaus'):
	"""

	"""
	if mode=='spec':
		f = np.loadtxt('../data/WASP-076/trans-eqpt_WASP-076_1.00_+1.0_1.50_model.txt')
		f = np.loadtxt('../data/WASP-127/trans-eqpt_WASP-127_1.00_+1.0_1.50_model.txt')
		wl, trans = f[:,0][1000:4000],f[:,1][1000:4000]

		f2    = interp1d(wl*1000,trans, kind='linear',fill_value=0,bounds_error=False)
		xneid = np.arange(380, 920, 650./res) # split the resolution between Na and K wavelengths
		yneid = f2(xneid)
	elif mode =='gaus':
		xneid = np.arange(380, 920, 650./res)
		yneid = np.zeros(len(xneid))
		# add gaussians 
		yneid += D2 * gaussian(xneid,588.9930,0.0619) #D2 Na
		yneid += D1 * gaussian(xneid,589.5903,0.0680) #D1

		yneid += (1/NaK_rat) * D2 * gaussian(xneid,766.48991,0.0619) #D2 K
		yneid += (1/NaK_rat) * D1 * gaussian(xneid,769.89645,0.0680) #D1

	# resample to NEID resolution
	sigma      = 1/(snr_1/np.sqrt(2))/np.sqrt(nspec*nvis)
	yneid_data = yneid + sigma * np.random.standard_normal(len(yneid))

	# Bin data
	bins, means, stds = bindat(xneid,yneid_data,nbins)

	return xneid, yneid, yneid_data, bins, means, sigma/np.sqrt(len(xneid)/nbins), sigma


def plot_simspec(planets,dat):
	fig, ax = plt.subplots(3,2,figsize=(9,6), sharex=True, sharey=False)
	for i in range(3):
		for j, planet in enumerate(planets):
			subdat = dat[planet]
			iplot = np.where((subdat['x_%s'%i] > 765) & (subdat['x_%s'%i] < 772))[0]
			ax[i,j].errorbar(subdat['x_%s' %i][iplot],100*subdat['ydat_%s' %i][iplot],100*subdat['sigma'],fmt='.',
				c='lightgray',zorder=0,
				label=subdat['NaK_rat_%s'%i])
			iplotbin = np.where((subdat['xbin_%s'%i] > 765) & (subdat['xbin_%s'%i] < 772))[0]
			ax[i,j].errorbar(subdat['xbin_%s'%i][iplotbin],100*subdat['ybin_%s'%i][iplotbin],100*subdat['sbin_%s'%i],\
				fmt='.',c='k',zorder=10) 
			ax[i,j].plot(subdat['x_%s'%i][iplot],100*subdat['y_%s'%i][iplot], 'r', zorder=9)
			if j==0:
				ax[i,j].set_ylim(-1.5,0.8)
			elif j==1:
				ax[i,j].set_ylim(-2.5,1.0)

	fig.subplots_adjust(hspace=0,wspace=0.15,left=0.15,bottom=0.14)
	ax[2,0].set_xlabel('Wavelength (nm)')
	ax[2,1].set_xlabel('Wavelength (nm)')
	ax[1,0].set_ylabel(r'$F_{in}/F_{out}$ -1 (%)')
	ax[0,0].set_xlim(766,770.5)

	ax[0,0].set_title(planets[0])
	ax[0,1].set_title(planets[1])

	ax[0,0].text(766.2,-1.1,'Na/K = %s' %subdat['NaK_rat_0'])
	ax[1,0].text(766.2,-1.1,'Na/K = %s' %subdat['NaK_rat_1'])
	ax[2,0].text(766.2,-1.1,'Na/K = %s' %subdat['NaK_rat_2'])

	plt.savefig('./simulated_trans_spec.pdf')
	plt.savefig('./simulated_trans_spec.eps')



def calc_nexp(tdur, texp):
	"""
	calculate number exposures can do in one transit observation
	based on neid overhead calculator, assumine one visit

	Time = (# visits to target) * (180s + (Exposure time in s) * (# exposures per visit) + 30s * (# exposures per visit - 1))

	inputs: tdur: duration of transit  (s)
			texp: exposure time (s)
	"""
	return (tdur + 30) / (180 + texp + 30)


def calc_time(texp,nexp):
	"""
	texp (s)
	nexp

	return time in hours
	"""
	time = (180 + texp) * nexp + 30 * (nexp -1)

	return time/3600.

###################### Generate data for 
if __name__=='__main__':
	planets = ['WASP-076b', 'WASP-127b']
	NaK_rats = [0.5, 1.0,2.0]
	dat = {}
	for planet in planets: 
		dat[planet] = {}
		if planet == 'WASP-076b':  # 
			dat[planet]['D2'], dat[planet]['D1'] =  -0.00373, -0.00508
			dat[planet]['texp']  = 900 #s
			dat[planet]['snr_1'] = 64.73*np.sqrt(5) # 9.5 v band mag, 6250 Teff, 767nm wavelength
			dat[planet]['nvis']  = 4
			dat[planet]['nbins'] = 10000
			dat[planet]['tdur']  = 3.683*3600 #hours
			dat[planet]['nspec'] = int(calc_nexp(dat[planet]['tdur'], dat[planet]['texp']))
		elif planet == 'WASP-127b':
			dat[planet]['D2'], dat[planet]['D1'] =   -0.0114, -0.0073
			dat[planet]['texp']  = 900 #s
			dat[planet]['snr_1'] = 49.6 *np.sqrt(5)# 10.16 v band mag, 5750 Teff
			dat[planet]['nvis']  = 2
			dat[planet]['nbins'] = 10000
			dat[planet]['tdur']  = 4.3*3600#hours to seconds
			dat[planet]['nspec'] = int(calc_nexp(dat[planet]['tdur'], dat[planet]['texp']))
		for i in range(3):
			dat[planet]['NaK_rat_%s' %i] = NaK_rats[i]
			# calc sim
			dat[planet]['x_%s' %i], dat[planet]['y_%s'%i], dat[planet]['ydat_%s' %i],\
					dat[planet]['xbin_%s' %i], dat[planet]['ybin_%s' %i], dat[planet]['sbin_%s' %i], \
					dat[planet]['sigma'] = calc_spec(dat[planet]['D1'], dat[planet]['D2'], dat[planet]['snr_1'], dat[planet]['nspec'], \
									dat[planet]['nvis'], res,dat[planet]['NaK_rat_%s' %i], \
									nbins=dat[planet]['nbins'],mode='gaus')



	###################### Plot
	plot_simspec(planets,dat)
	for planet in planets:
		subdat = dat[planet]
		print(planet)
		print('sigma of %s spectra and %s visits:  %s' %(subdat['nspec'], subdat['nvis'], subdat['sigma']))
		print('sigma after binning %s: ' %subdat['sbin_1'])
		tot_time = subdat['nvis']*calc_time(subdat['texp'],2*subdat['nspec'])
		print('total time in hours: %s  nights: %s ' %(tot_time,tot_time/10.0))



