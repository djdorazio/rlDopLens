import numpy as np




#Physical Constants
G    = 6.67384*10**(-8)
c    = 2.99792458*10**(10)
Msun = 1.989*10**33
day2sec = 24.*3600.
day2yr = 1./365.25
yr2sec = 3.154*10.**7



#err_fac = 0.01


def Get_FitData(filename, err_fac, err_fac_zm, day_intervl, day_intervl_zm):
	#################################
	### Now import MBHB LC
	#################################
	t_rday   =  np.genfromtxt("1918+49_lc_corrected.dat", comments="#", delimiter=" ", usecols=0)
	mag      =  np.genfromtxt("1918+49_lc_corrected.dat", comments="#", delimiter=" ", usecols=1)

	mag_errs = err_fac*mag

	mag_err_flare = mag_errs
	# for i in range(len(mag_errs)):
	# 	if (t_rday[i] > 199 and t_rday[i] < 205):
	# 		mag_err_flare[i] = mag_errs[i]/5.0



	t_rdayp   =  np.genfromtxt("1918+49_lc_corrected.dat", comments="#", delimiter=" ", usecols=0)
	magp      =  np.genfromtxt("1918+49_lc_corrected.dat", comments="#", delimiter=" ", usecols=1)



	######################################################################
	######################################################################
	print "BINNING DATA"
	## put in time order
	tS   = zip(t_rdayp,magp,mag_errs)
	tS.sort()
	TtS  = np.transpose(tS)
	t_srt   =  np.array(TtS[0])
	mag_srt    =  np.array(TtS[1])
	mag_sigsrt =  np.array(TtS[2])


	###### get average value for each cluster of data points in time
	day_intervlend = day_intervl
	iendd = np.where(t_srt>t_srt[len(t_srt)-1]-day_intervlend)[0][0]
	iseg = []
	iseg.append(-1)
	i=0
	while (i<=iendd-1):
		j=1
		while (j<=len(t_srt)-1-i):
			if (t_srt[i+j] > 198.0 and t_srt[i+j] < 210.0):
				day_intervl_use = day_intervl_zm
			else:
				day_intervl_use = day_intervlend
				#print "long"
			if (abs((t_srt[i+j] - t_srt[i])) > day_intervl_use ):
				iseg.append(j+i)
				i = i+j
				j=len(t_srt)
				
			else:
				j=j+1
				#print j
		print i
	iseg.append(len(t_srt)-1)

	t_avg = []
	mag_avg = []


	#VarMn_test_W1 = []
	#VarMn_test_W2 = []


	for i in range(0 , len(iseg)-1):
		t_avg.append(np.mean(t_srt[iseg[i]+1:iseg[i+1]+1]))

		mag_avg.append(np.mean(mag_srt[iseg[i]+1:iseg[i+1]+1]))
		

		
		#Nseg = len(V_sigsrt[iseg[i]+1:iseg[i+1]]) + 1


		


		# ## see if variance of means or mean of variances is larger
		# VarMn_test_W1.append(np.mean(W1_sig[iseg[i]+1:iseg[i+1]+1])/np.std(W1_mag[iseg[i]+1:iseg[i+1]+1]))
		# VarMn_test_W2.append(np.mean(W2_sig[iseg[i]+1:iseg[i+1]+1])/np.std(W2_mag[iseg[i]+1:iseg[i+1]+1]))


	t_avg = np.array(t_avg)
	mag_avg = np.array(mag_avg)

	mag_err_avg = err_fac*mag_avg

	mag_err_flare_avg = mag_err_avg
	# for i in range(len(mag_err_avg)):
	# 	if (t_avg[i] > 198 and t_avg[i] < 208):
	# 		mag_err_flare_avg[i] = mag_err_avg[i]/5.0


	print "DONE BINNING DATA"
	######################################################################
	######################################################################



	#sav golay
	import scipy as sc
	from scipy import signal


	bkgnd = sc.signal.savgol_filter(mag, len(mag)/3, 3)
	bkgmean = np.mean(bkgnd)

	bkgnd2 = sc.signal.savgol_filter(mag, len(mag)/500, 7)


	flr = (mag)/bkgnd #np.zeros(len(mag))#mag - bkgnd
	tzm = t_rday#np.zeros(len(mag))#t_rday
	for i in range(len(flr)):
		if (t_rday[i]<203.0-5.0 or t_rday[i]>203.+5.0):
			flr[i]=0.0
			tzm[i]=0.0
		# if (t_rday[i]>203.0-5.0 or t_rday[i]<203.+5.0):
		# 	flr[i]=mag[i]
		# 	tzm[i]=t_rday[i]

	flrzoom = np.trim_zeros(flr)
	tzoom   = np.trim_zeros(tzm)







	######################################################################
	######################################################################
	print "BINNING ZOOM DATA"
	## put in time order
	mag_err_zm = err_fac_zm*flrzoom

	tS   = zip(tzoom,flrzoom,mag_err_zm)
	tS.sort()
	TtS  = np.transpose(tS)
	t_srt   =  np.array(TtS[0])
	mag_srt    =  np.array(TtS[1])
	mag_sigsrt =  np.array(TtS[2])


	###### get average value for each cluster of data points in time
	day_intervl = day_intervl_zm
	iendd = np.where(t_srt>t_srt[len(t_srt)-1]-day_intervl)[0][0]
	iseg = []
	iseg.append(-1)
	i=0
	while (i<=iendd-1):
		j=1
		while (j<=len(t_srt)-1-i):
			if (abs((t_srt[i+j] - t_srt[i])) > day_intervl ):
				iseg.append(j+i)
				i = i+j
				j=len(t_srt)
				
			else:
				j=j+1
				#print j
		print i
	iseg.append(len(t_srt)-1)

	t_avg_zm = []
	mag_avg_zm = []


	for i in range(0 , len(iseg)-1):
		t_avg_zm.append(np.mean(t_srt[iseg[i]+1:iseg[i+1]+1]))
		mag_avg_zm.append(np.mean(mag_srt[iseg[i]+1:iseg[i+1]+1]))
		

	t_avg_zm = np.array(t_avg_zm)
	mag_avg_zm = np.array(mag_avg_zm)

	mag_err_avg_zm = err_fac_zm*mag_avg_zm

	print "DONE BINNING ZOOM DATA"
	######################################################################
	######################################################################


	dicout = {"bkgnd":bkgnd, "bkgmean":bkgmean, "tzoom":tzoom, "flrzoom":flrzoom, "t_rdayp":t_rdayp, "magp":magp, "mag_errs":mag_errs, "t_avg":t_avg, "mag_avg":mag_avg, "mag_err_avg":mag_err_avg, "mag_err_flare":mag_err_flare, "mag_err_flare_avg":mag_err_flare_avg, "t_avg_zm":t_avg_zm, "mag_avg_zm":mag_avg_zm, "mag_err_avg_zm":mag_err_avg_zm,  "mag_err_zm":mag_err_zm}
	return dicout




