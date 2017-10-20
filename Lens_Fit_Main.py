import numpy as np
import matplotlib
# import scipy as sc
import scipy.optimize as opti
#from optimize import brentq
#matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 16})

import math as ma


## FOR READING AND BINNING DATA
import Lens_Fit_ReadData
from Lens_Fit_ReadData import *

## MODELS and Err Funcs
import Lens_Fit_FitFuncs
from Lens_Fit_FitFuncs import *


## Dynesty an emcee specific funcs
import Lens_Fit_dynestyFuncs
from Lens_Fit_dynestyFuncs import *

## parellel for dynesty
from multiprocessing import Pool
ppool = Pool(processes=8)
procs = 8



### FITTING OPTIONS
FitIt = True
##
FitMeth = "dynesty"

#FitMeth = "FminFit"
##
full_ecc = False ##fit only the lens flare and feed best fit params into Doppler curve
DopLens_Fit = True #fit enitre Dop+Lens to entire lightcurve

### PLOTTING OPTIONS
plt_res = True #plot best fit
plt_ICs = False


#Physical Constants
G    = 6.67384*10**(-8)
c    = 2.99792458*10**(10)
Msun = 1.989*10**33
day2sec = 24.*3600.
day2yr = 1./365.25
yr2sec = 3.154*10.**7












################################################
###IMPORT DATA
################################################
print "Importing Data"
#MAKE USER INPUT LATER
#Kep KLS AGN
filename="1918+49_lc_corrected.dat"
err_fac = 0.01
err_fac_zm = 0.005
day_intervl = 20.0
day_intervl_zm = 0.5
dats = Get_FitData(filename, err_fac, err_fac_zm, day_intervl, day_intervl_zm)

bkgnd  = dats["bkgnd"]
bkgmean  = dats["bkgmean"]

##All data entire dt of data
t_rdayp  = dats["t_rdayp"]
magp    = dats["magp"]
mag_errs = dats["mag_errs"]
mag_err_flare = dats["mag_err_flare"]  ## diff error bars in flare

##All data, Zoomed in on flare
tzoom       = dats["tzoom"]
flrzoom     = dats["flrzoom"]
mag_err_zm  = dats["mag_err_zm"]



##Binned over entire dt of data
t_avg       = dats["t_avg"]
mag_avg     = dats["mag_avg"]
mag_err_avg = dats["mag_err_avg"]
mag_err_flare_avg = dats["mag_err_flare_avg"]  ## diff error bars in binned flare

##Binned, Zoomed in on flare
t_avg_zm       = dats["t_avg_zm"]
mag_avg_zm     = dats["mag_avg_zm"]
mag_err_avg_zm = dats["mag_err_avg_zm"]
print "Data Imported"
################################################
###^^IMPORT DATA^^
################################################
















######################
###DEBUG/ PLOT ICs
######################

##ecc Dop params
#pDop = [ecc,cosw,T0,KK,Tbin,fs]
Mtst = 8.
qtst = 0.1
ecctst = 0.5
coswtst = 0.0

ecctst1 = 0.2
ecctst2 = 0.5

coswtst1 = 0.5
coswtst2 = 0.0

Tbintst = 450.
CosItst = 0.99
KKtst = 0.1
fstst = 0.5
alptst = 1.5

ph0tst = -2.*np.pi*247/Tbintst

pDop0 = [Mtst, qtst, ecctst, coswtst, Tbintst*0.0, KKtst, Tbintst, fstst, alptst]

pDopLens0 = [Mtst, qtst, ecctst2, coswtst, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]
pDopLens1 = [Mtst, qtst, ecctst2, coswtst1, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]
pDopLens2 = [Mtst, qtst, ecctst2, coswtst2, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]

plens = [Mtst, Tbintst, qtst, CosItst, 0.0]

plens_ecc = [Mtst, Tbintst, qtst, CosItst, 0.0, ecctst, coswtst]
plens_ecc1 = [Mtst, Tbintst, qtst, CosItst, 0.0, ecctst, coswtst1]
plens_ecc2 = [Mtst, Tbintst, qtst, CosItst, 0.0, ecctst, coswtst2]


tt = np.linspace(195.0, 210.0, 1000)
tl = np.linspace(0.0, 1000.0, 1000)


# plt.figure()
# plt.plot(tl, MagPSprm(plens, tl*day2sec))
# plt.plot(tl, MagPSsec(plens, tl*day2sec))
# plt.show()


# plt.figure()
# plt.plot(tl, MagPSprm_ecc(plens_ecc, tl*day2sec))
# plt.plot(tl, MagPSsec_ecc(plens_ecc, tl*day2sec))
# #plt.plot(tl, MagPSsec_ecc(plens_ecc1, tl*day2sec))
# #plt.plot(tl, MagPSsec_ecc(plens_ecc2, tl*day2sec))
# plt.show()

#p0E = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]
#from fmin on flare zoom
p0E = [  7.96595983e+00,   4.48578302e+02,   8.82467081e-02, 9.89532469e-01,  -3.44860788e+00, 6.03603989e-01, 1.63034655e-04, 4.51238859e-01]
	#pDopLens0 = [Mtst, qtst, ecctst2, coswtst, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]

#From fmin DopLens fit
#p0E = [8.2157630055802926, 510.75101131540401, 0.10324063269547901, 0.92412490822557425, -183.45673927598466, 0.52290916488381289, 0.00024182304074598073, 0.50658607995971572]
#Dop0 = [p0E[0], p0E[2], p0E[5], p0E[6], p0E[4]*p0E[1]/2./np.pi, p0E[3], p0E[1], p0E[7], alptst]
	
# Dop0 = [ 7.96465683e+00,   1.01889801e-01,   5.67425637e-01,
#          2.80781763e-05,  -2.41601578e+02,   9.91013118e-01,
#          4.43887801e+02,   4.37456235e-01,   1.53406427e+00]
Dop0 = [ 7.96787162e+00,   1.06195167e-01,   5.68348476e-01,
         6.49502429e-06,  -2.50977976e+02,   9.91452061e-01,
         4.53272107e+02,   4.15220989e-01,   1.47790115e+00]



if (plt_ICs):
	print "PLOTTING ICs"

	plt.figure(figsize=[10,5])

	plt.subplot(211)
	#plt.plot(tl, DopLum(pDop0, tl*day2sec))
	plt.scatter(t_rdayp, magp/bkgmean, s=1, alpha=0.5)
	plt.errorbar(t_rdayp, magp/bkgmean, yerr=mag_err_flare/bkgmean, alpha=0.5,linestyle="none")

	plt.scatter(t_avg, mag_avg/bkgmean, s=1, alpha=0.5, color='black', zorder = 10)
	plt.errorbar(t_avg, mag_avg/bkgmean, yerr=mag_err_avg/bkgmean, alpha=0.5, color='black',linestyle="none", zorder = 10)

	plt.plot(tl, DopLumLens(Dop0, tl), linestyle="--", zorder = 20)








	plt.subplot(212)


	plt.scatter(tzoom, flrzoom, s=1, alpha=0.5)
	plt.errorbar(tzoom, flrzoom, yerr=mag_err_zm, alpha=0.5,linestyle="none")


	plt.scatter(t_avg_zm, mag_avg_zm, s=1, alpha=0.5, color='black', zorder = 10)
	plt.errorbar(t_avg_zm, mag_avg_zm, yerr=mag_err_avg_zm, alpha=0.5, color='black',linestyle="none", zorder = 10)
	
	plt.plot(tl, MagPS_ecc(p0E, tl), linestyle="--", zorder = 20)

	plt.xlim(190,220)

	#plt.plot(tl, DopLumLens(pDopLens1, tl*day2sec), linestyle=":")
	#plt.plot(tl, DopLumLens(pDopLens2, tl*day2sec), linestyle="-.")
	
	#plt.show()
	plt.savefig("DopLensICs.png")
	######################
	###^DEBUG/PLOT ICS^
	######################


























if (FitIt):

	######################
	###FITTING
	######################
	if (FitMeth =="FminFit"):
		if (full_ecc): ##fit only the lens flare and feed best fit params into Doppler curve
			#p =   M,    Prest,   q,    CosI,    ph0, ecc,    cosw, fs = p
			p0E = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]
			print "Fitting Lens Flare Zoom in"
			popt  = sc.optimize.fmin(Mag_ecc_Fmin_Err2,  p0E, args=(t_avg_zm, mag_avg_zm, mag_err_avg_zm), full_output=1, disp=False, ftol=0.001)[0]
			print "DONE Fitting Lens Flare Zoom in"
			DopOpt = [popt[0], popt[2], popt[5], popt[6], popt[4]*popt[1]/2./np.pi, popt[3], popt[1], popt[7], 1.0]
			#pDopLens0 = [Mtst, qtst, ecctst2, coswtst, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]


		elif (DopLens_Fit): #fit enitre Dop+Lens to entire lightcurve
			p0E = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]
			
			#pDopLens0 = [Mtst, qtst, ecctst2, coswtst, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]
			Dop0 = [p0E[0], p0E[2], p0E[5], p0E[6], p0E[4]*p0E[1]/2./np.pi, p0E[3], p0E[1], p0E[7], alptst]
				
			print "Fitting DOPLENS Avgs"
			DopOpt  = sc.optimize.fmin(DopLens_Fmin_Err2,  Dop0, args=(t_avg, mag_avg/bkgmean, mag_err_flare_avg/bkgmean), full_output=1, disp=False, ftol=0.001)[0]
			print "DONE Fitting DOPLENS Avgs"

			#p0E = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]
			popt = [DopOpt[0], DopOpt[6], DopOpt[1], DopOpt[5], DopOpt[4]/DopOpt[6]*2.*np.pi, DopOpt[2], DopOpt[3], DopOpt[7]]
		else:
			print "Specify what we are fitting"
			brk

	elif (FitMeth == "dynesty"):
		from scipy.stats import truncnorm
		# seed the random number generator
		np.random.seed(2)
		import dynesty
		if (full_ecc): ##fit only the lens flare and feed best fit params into Doppler curve
			#p =   M,    Prest,   q,    CosI,    ph0, ecc,    cosw, fs = p
			print "DYNESTY Fitting Lens Flare Zoom in"
			#DYNESTY STUFF HRE
			###setup params
			param_names = [r"$M$", r"$P_{\rm{r}}$", r"$q$",  r"$\cos{I}$", r"$\phi_0$", r"$e$", r"$\cos{w}$", r"$f_2$"]
			p0 = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]

			ndim = len(p0)
			nwalkers = ndim*4

			# initialize dynamic nested sampler
			dsampler = dynesty.DynamicNestedSampler(ln_Fzoom_posterior_dyn, ptform_Zoom, ndim, logl_args=(t_avg_zm, mag_avg_zm, mag_err_avg_zm), sample='rwalk', bound='multi', update_interval=6.*ndim, walks=50, queue_size=procs, pool=ppool)


			#popt  = 
			#print "DONE Fitting Lens Flare Zoom in"
			#DopOpt = [popt[0], popt[2], popt[5], popt[6], popt[4]*popt[1]/2./np.pi, popt[3], popt[1], popt[7], 1.0]


		elif (DopLens_Fit): #fit enitre Dop+Lens to entire lightcurve
			p0E = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]
			#pDopLens0 = [Mtst, qtst, ecctst2, coswtst, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]
			
			param_names = [r"$M$", r"$q$", r"$e$", r"$\cos{w}$", r"$\phi_0$", r"$\cos{I}$", r"$P_{\rm{r}}$", r"$f_2$", r"$\alpha_{\nu}$"]
			p0 = [p0E[0], p0E[2], p0E[5], p0E[6], p0E[4]*p0E[1]/2./np.pi, p0E[3], p0E[1], p0E[7], alptst]
				
			ndim = len(p0)
			nwalkers = ndim*4

			# initialize dynamic nested sampler
			dsampler = dynesty.DynamicNestedSampler(ln_DopLens_posterior_dyn, ptform_DopLens, ndim, logl_args=(t_avg_zm, mag_avg_zm, mag_err_avg_zm), sample='rwalk', bound='multi', update_interval=6.*ndim, walks=50, queue_size=procs, pool=ppool)


			# print "DYNESTY Fitting DOPLENS Avgs"
			

			# #DopOpt  = 
			# print "DONE Fitting DOPLENS Avgs"

			# #p0E = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]
			# popt = [DopOpt[0], DopOpt[6], DopOpt[1], DopOpt[5], DopOpt[4], DopOpt[2], DopOpt[3], DopOpt[7]]
		else:
			print "Specify what we are fitting"
			brk



		# run with 120 initial live points until dlogz=0.01
		# add 100 live points at a time
		# use default weight function with 100% posterior weight
		# use automated stopping criteria (100% posterior weight)
		print "Running dsampler"
		dsampler.run_nested(nlive_init=10*ndim, nlive_batch=100, dlogz_init=0.01, wt_kwargs={'pfrac': 1.0})


		##ANALYSIS
		print "Analysing Dynesty Results"
		dres = dsampler.results  # store results
		with open("LSQ12dlf_dynesty_DopLens%g_Zoom%g.pkl" %(DopLens_Fit, full_ecc), "wb") as out:
			pickle.dump(dres, out)

		    
		# diagnostic summary plot
		from dynesty import plotting as dyplot

		print "Summary Plot"
		fig, ax = dyplot.runplot(dres, logplot=True)
		fig.tight_layout()
		fig.savefig("dyntest_sumry.png")


		# trace plot
		print "Trace Plot"
		fig1, ax1 = dyplot.traceplot(dres, labels=param_names, show_titles=True)
		fig1.tight_layout()
		fig1.savefig("dyntest_trace.png")



		print "Corner Plot"
		fig2, ax2 = dyplot.cornerplot(dres, show_titles=True, labels=param_names)
		fig2.savefig("dyntest_corner.png")

else:
	popt = p0E
	DopOpt = Dop0














if (plt_res):


	######################
	###PLOT RESULTS
	######################
	plt.figure(figsize=[10,7])



	plt.subplot(211)
	plt.scatter(t_rdayp, magp, s=1, alpha=0.5)

	plt.scatter(t_avg, mag_avg, s=1, alpha=0.5, color='black')
	plt.errorbar(t_avg, mag_avg, yerr=mag_err_avg, alpha=0.5, color='black', linestyle="none")

	plt.plot(t_rdayp, bkgnd, color="blue")

	plt.plot(tl, DopLumLens(DopOpt, tl)*bkgmean, color="orange", alpha=1.0, linewidth=3, linestyle="--", zorder=10)

	#plt.xlim(t_rdayp[0], t_rdayp[len(t_rdayp)-1])
	plt.ylabel("counts per second")
	plt.xlabel("rest frame days")








	plt.subplot(212)
	plt.title("Zoom-in on flare")

	plt.scatter(tzoom, flrzoom, s=1, alpha=0.5)
	plt.errorbar(tzoom, flrzoom, yerr=mag_err_zm, alpha=0.5,linestyle="none")


	plt.scatter(t_avg_zm, mag_avg_zm, s=1, alpha=0.5, color='black', zorder = 10)
	plt.errorbar(t_avg_zm, mag_avg_zm, yerr=mag_err_avg_zm, alpha=0.5, color='black',linestyle="none", zorder = 10)


	#	plt.plot(tt, MagPS_ecc(popt, tt), color='orange', alpha=1.0, linewidth=3, linestyle="--", zorder=10)

	# if (full_ecc):
	# 	plt.plot(tt, MagPS_ecc(popt, tt), color='orange', alpha=1.0, linewidth=3, linestyle="--", zorder=10)
	# 	#p0E = [Mtst, Tbintst, qtst, CosItst, ph0tst, ecctst, coswtst, fstst]
	# 	plt.figtext(0.8, 0.47, "Best fit parameters:", fontsize=12)
	# 	plt.figtext(0.8, 0.43, r"$M = 10^{%g} M_{\odot}$" %popt[0], fontsize=12)
	# 	plt.figtext(0.8, 0.39, r"$P = %g$ days"           %popt[1], fontsize=12)
	# 	plt.figtext(0.8, 0.35, r"$q = %g$"                %popt[2], fontsize=12)
	# 	plt.figtext(0.8, 0.31, r"$ecc = %g$"              %popt[5], fontsize=12)
	# 	plt.figtext(0.8, 0.27, r"$cosw = %g$"             %popt[6], fontsize=12)
	# 	plt.figtext(0.8, 0.23, r"$\cos{I} = %g$"          %popt[3], fontsize=12)
	# 	plt.figtext(0.8, 0.19, r"$\phi_0 = %g$"           %popt[4], fontsize=12)
	# 	plt.figtext(0.8, 0.15, r"$f_2 = %g$"              %popt[7], fontsize=12)
	# if (DopLens_Fit):
	plt.plot(tt, DopLumLens(DopOpt, tt), color="orange", alpha=1.0, linewidth=3, linestyle="--", zorder=10)
	#pDopLens0 = [Mtst, qtst, ecctst2, coswtst, Tbintst*0.0, CosItst, Tbintst, fstst, alptst]
	plt.figtext(0.8, 0.47, "Best fit parameters:", fontsize=12)
	plt.figtext(0.8, 0.43, r"$M = 10^{%g} M_{\odot}$" %DopOpt[0], fontsize=12)
	plt.figtext(0.8, 0.39, r"$P = %g$ days"           %DopOpt[6], fontsize=12)
	plt.figtext(0.8, 0.35, r"$q = %g$"                %DopOpt[1], fontsize=12)
	plt.figtext(0.8, 0.31, r"$ecc = %g$"              %DopOpt[2], fontsize=12)
	plt.figtext(0.8, 0.27, r"$cosw = %g$"             %DopOpt[3], fontsize=12)
	plt.figtext(0.8, 0.23, r"$\cos{I} = %g$"          %DopOpt[5], fontsize=12)
	plt.figtext(0.8, 0.19, r"$\phi_0 = %g$"           %DopOpt[4], fontsize=12)
	plt.figtext(0.8, 0.15, r"$f_2 = %g$"              %DopOpt[7], fontsize=12)
	plt.figtext(0.8, 0.11, r"$\alpha = %g$"           %DopOpt[8], fontsize=12)



	plt.ylabel(r"$\mathcal{M}$")
	plt.xlabel("rest frame days")


	##DEBUG
	#plt.plot(tt, MagPS_ecc(popt, tt), color="black")


	plt.xlim(tzoom[0]-10, tzoom[len(tzoom)-1]+10)
	plt.tight_layout()

	if (full_ecc):
		Savename = "BestFits_"+FitMeth+"_Flareonly_eccON.png" 
	# Savename = Savename.replace('.', 'p')
	# Savename = Savename.replace('ppng', '.png')
	else:
		Savename = "BestFits_"+FitMeth+"_DopandFlare_eccON.png" 
	plt.savefig(Savename)
	######################
	###^PLOT RESULTS^
	######################






