import numpy as np
from numpy import *

## MODELS
import Lens_Fit_FitFuncs
from Lens_Fit_FitFuncs import *







############################
## MAX LIK FNS for Simul FIT of W1 W2 and V0bands
############################
#PRIORS
def ptform_Zoom(params):
	M, Prest, q, CosI, ph0, ecc, cosw, fs = params
					
	#Transform from unit cube to actual param values

	M_out = 5.0 + M*6.0 ##from 10^(5) to 10^(11)  

	Prest_out = Prest*1000.0 + 350.00 # make repeat no sooner than at end of data = ~200 + 350 = 550 days

	q_out = q 

	CosI_out = CosI*2.0-1.0

	pho_out = ph0 *2.*ma.pi

	ecc_out = ecc

	cosw_out = cosw*2.0-1.0

	fs_out = fs

	x = np.array([M_out, Prest_out, q_out, CosI_out, pho_out, ecc_out, cosw_out, fs_out])		
	return x


#PRIORS
def ptform_DopLens(params):
	#M, q, ecc, cosw, T0, CosI, Prest, fs, alp = params
	M, q, ecc, cosw, T0, CosI, Prest, fs = params

					
	#Transform from unit cube to actual param values

	M_out = 8.0 + M*2.0 ##from 10^(5) to 10^(11)  
	
	q_out = q 

	ecc_out = ecc

	cosw_out = cosw*2.0-1.0

	T0_out = T0*4350.0#(T0*2.0 - 1.0)*3000.0 #days

	CosI_out = CosI*2.0-1.0

	Prest_out = Prest*4000.0 + 350.00 #days

	fs_out = fs

	#alp_out = 4.*(alp*2.0 - 1.0)

	#x = np.array([M_out, q_out, ecc_out, cosw_out, T0_out, CosI_out, Prest_out, fs_out, alp_out])		
	x = np.array([M_out, q_out, ecc_out, cosw_out, T0_out, CosI_out, Prest_out, fs_out])		

	return x





##ERR2 FUNCS (-LogLik now)
def Mag_ecc_Err2(p, t, y, dy):
	#print "EVAL", p
	#M, Prest, q, CosI, ph0, ecc, cosw, fs = p
	chi = 0.5*(y - MagPS_ecc(p, t) )/ dy
	chi2 = np.sum(chi*chi) 
	#print(chi2)
	return chi2




def DopLens_Err2(p, t, y, dy):
	#print "EVAL", p
	#Mbin, qq, ecc, cosw, T0, CosI, Tbin, fs, alp = params
	chi = 0.5*(y - DopLumLens(p, t) )/ dy
	chi2 = np.sum(chi*chi) 
	#print(chi2)
	return chi2

	# ##LogLik i made negative in liklihood function
	# #print(-negLogLik)
	# return negLogLik #pos chi2





##likliehoods (- Err2 Funcs)
def ln_Fzoom_likelihood(p, t, y, dy):
	return -Mag_ecc_Err2(p, t, y, dy)

def ln_DopLens_likelihood(p, t, y, dy):
	return -DopLens_Err2(p, t, y, dy)
		


##POSTERIORS
##for dynasty with ptform() at top of this file
def ln_Fzoom_posterior_dyn(p, t, y, dy):
	return ln_Fzoom_likelihood(p, t, y, dy)
	

def ln_DopLens_posterior_dyn(p, t, y, dy):
	return ln_DopLens_likelihood(p, t, y, dy)



#FMIN
def Mag_ecc_Fmin_Err2(p, t, y, dy):
	print "EVAL", p
	M, Prest, q, CosI, ph0, ecc, cosw, fs = p
	#if (M<8.5 or M>9.5 or q<0.0 or q>1.0 or CosI<-1.0 or CosI>1.0 or cosw<0.0 or cosw>1.0 or ecc<0.0 or ecc>1.0 or fs<0.0 or fs>1.0 or Prest<300.0 or Prest>10000.):
	if (M<0.0 or M>11.0 or q<0.0 or q>1.0 or CosI<-1.0 or CosI>1.0 or cosw<0.0 or cosw>1.0 or ecc<0.0 or ecc>1.0 or fs<0.0 or fs>1.0 or Prest<0.0 or Prest>10000.):
		return np.inf
	else:
		chi = 0.5*(y - MagPS_ecc(p, t) )/ dy
		chi2 = np.sum(chi*chi) 
		print(chi2)
		return chi2


def DopLens_Fmin_Err2(p, t, y, dy):
	print "EVAL", p
	#M, q, ecc, cosw, T0, CosI, Prest, fs, alp = p
	M, q, ecc, cosw, T0, CosI, Prest, fs = p
	#if (M<8.5 or M>9.5 or q<0.0 or q>1.0 or CosI<-1.0 or CosI>1.0 or cosw<0.0 or cosw>1.0 or ecc<0.0 or ecc>1.0 or fs<0.0 or fs>1.0 or Prest<300.0 or Prest>10000.):
	if (M<0.0 or M>11.0 or q<0.0 or q>1.0 or CosI<-1.0 or CosI>1.0 or cosw<0.0 or cosw>1.0 or ecc<0.0 or ecc>1.0 or fs<0.0 or fs>1.0 or Prest<0.0 or Prest>10000.):
		return np.inf
	else:
		chi = 0.5*(y - DopLumLens(p, t) )/ dy
		chi2 = np.sum(chi*chi) 
		print(chi2)
		return chi2



