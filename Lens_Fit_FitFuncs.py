import numpy as np
import math as ma
import scipy as sc



#Physical Constants
G    = 6.67384*10**(-8)
c    = 2.99792458*10**(10)
Msun = 1.989*10**33
day2sec = 24.*3600.
day2yr = 1./365.25
yr2sec = 3.154*10.**7

# #OPTIONS
# Prm_Lens = True
# qtst = 0.05
# NePrb = 1.0
# rhomx = 1.0
# NRsmax = 10.0
# Ng = 200

#################################
### Lensing PRobs and Timescales
#################################

def asep(M,Prest):
	return (Prest/(2.*(np.pi)))**(2./3) * (G*M)**(1./3.)


###ecc
def rEprm_ecc(t, ph0, M, Prest,q, CosI, ecc, cosw):
	#circ
	T0 = ph0/2./np.pi*Prest
	EEu = np.mod(EAn(t,ecc,T0,Prest), 2.*np.pi)
	fot = 2.*np.arctan2(np.sqrt(1.+ecc)*np.tan(EEu/2.), np.sqrt(1.- ecc))  ##arctan2(top, bot)
	ww =  ma.acos(cosw)#np.arccos(cosw) ##careful ##arg of pericenter
	roa = (1.-ecc*ecc)/(1.+ecc*np.cos(fot + ww))
	return np.nan_to_num( np.sqrt(    4.*G*M/(1.+q)/c/c * (asep(M,Prest)*roa * CosI* -np.sin(fot + ww))    ) )


## time dep value of r_e/a  for sec lens
def rEsec_ecc(t, ph0, M, Prest,q, CosI, ecc, cosw):
	#circ
	T0 = ph0/2./np.pi*Prest
	EEu = np.mod(EAn(t,ecc,T0,Prest), 2.*np.pi)
	fot = 2.*np.arctan2(np.sqrt(1.+ecc)*np.tan(EEu/2.), np.sqrt(1.- ecc))  ##arctan2(top, bot)
	ww =  ma.acos(cosw)#np.arccos(cosw) ##careful ##arg of pericenter
	roa = (1.-ecc*ecc)/(1.+ecc*np.cos(fot + ww))
	return np.nan_to_num(  np.sqrt(     4.*G*M/(1.+1./q)/c/c * asep(M,Prest)*roa * CosI* -np.sin( fot + ww + np.pi )    )  )



def usec_ecc(t, ph0, M, Prest,q, CosI, ecc, cosw):
	### Solve for Eccentric Anamoly and then true eanomaly f(t)
	T0 = ph0/2./np.pi*Prest
	EEu = np.mod(EAn(t,ecc,T0,Prest), 2.*np.pi)
	fot = 2.*np.arctan2(np.sqrt(1.+ecc)*np.tan(EEu/2.), np.sqrt(1.- ecc))  ##arctan2(top, bot)
	ww = ma.acos(cosw)#np.arccos(cosw) ##careful ##arg of pericenter
	re = asep(M,Prest)*(1.-ecc*ecc)/(1.+ecc*np.cos(fot + ww))
	us = re/rEprm_ecc(t, ph0, M, Prest,q, CosI, ecc, cosw) * np.sqrt( (np.cos(fot + ww))**2  + (1.0-CosI**2)*(np.sin(fot + ww))**2 )
	# div by primary Einstein radius
	return us



def uprm_ecc(t, ph0, M, Prest,q, CosI, ecc, cosw):
	### Solve for Eccentric Anamoly and then true eanomaly f(t)
	T0 = ph0/2./np.pi*Prest
	EEu = np.mod(EAn(t,ecc,T0,Prest), 2.*np.pi)
	fot = 2.*np.arctan2(np.sqrt(1.+ecc)*np.tan(EEu/2.), np.sqrt(1.- ecc))  ##arctan2(top, bot)
	ww =  ma.acos(cosw)#np.arccos(cosw) ##careful ##arg of pericenter
	re = asep(M,Prest)*(1.-ecc*ecc)/(1.+ecc*np.cos(fot + ww))
	up = re/rEsec_ecc(t, ph0, M, Prest,q, CosI, ecc, cosw) * np.sqrt( (np.cos(fot + ww + np.pi))**2  + (1.0-CosI**2)*(np.sin(fot + ww + np.pi))**2 )
	# div by secondary Einstein radius
	return up









def MagPSsec_ecc(p, t):
	#plens = [Mbin, Tbin, qq, CosI, T0*2.*np.pi/Tbin]
	M, Prest, q, CosI, ph0, ecc, cosw = p
	M_units = 10.**M*Msun
	Prest_units = Prest*day2sec
	t_units = t*day2sec

	u = np.array(usec_ecc(t_units, ph0, M_units, Prest_units, q, CosI, ecc, cosw) )

	Mag = np.nan_to_num(  (u*u + 2.)/(u*np.sqrt(u*u + 4.)) )

	izero = np.where(Mag==0.0)[0]
	Mag[izero] = 1.0

	return Mag

def MagPSprm_ecc(p, t):
	M, Prest, q, CosI, ph0, ecc, cosw = p
	M_units = 10.**M*Msun
	Prest_units = Prest*day2sec
	t_units = t*day2sec
	u = np.array(uprm_ecc(t_units, ph0, M_units, Prest_units, q, CosI, ecc, cosw) )
	Mag = np.nan_to_num(  (u*u + 2.)/(u*np.sqrt(u*u + 4.)) )
	izero = np.where(Mag==0.0)[0]
	Mag[izero] = 1.0

	return Mag


def MagPS_ecc(p, t):
	M, Prest, q, CosI, ph0, ecc, cosw, fs = p
	pcomp = [M, Prest, q, CosI, ph0, ecc, cosw]
	return fs * MagPSsec_ecc(pcomp, t) +  (1.-fs)#+ (1.-fs)*MagPSprm_ecc(pcomp, t)







##############################
#### Beaming Model
##############################

### SOLVE FOR ECCENTRIC ANAMOLY	
def EEfunc(EEa, t, ecc,T0,Tbin):
	return np.absolute(2.*np.pi/Tbin*(t-T0) - EEa + ecc*np.sin(EEa))
	

		
def EAn(t, ecc,T0,Tbin):
	EEres = []
	for t_i in t:
		EEres.append(sc.optimize.fmin(EEfunc, 0., args=(t_i,ecc,T0,Tbin),disp=False ))
#	#print(EEres)
	return np.transpose(EEres)[0]




def DopLum(params, t):
	Mbin, qq, ecc, cosw, T0, KK, Tbin, fs, alp = params
	
	vmean=0.0
	Mbin = 10.**Mbin*Msun
	Tbin = Tbin*day2sec  #(put in sec)
	T0 = T0*day2sec #- 0.5*Tbin
	KK = KK*c       # put in cm/s
	t = t*day2sec
	
	semimaj = (np.sqrt(G*Mbin)*Tbin/(1+0.2784)/2./np.pi)**(2./3.)

	
	#### Inc measured from FACE ON

	### Solve for Eccentric Anamoly and then f(t)
	EEu = np.mod(EAn(t,ecc,T0,Tbin), 2.*np.pi)
	fot = 2.*np.arctan2(np.sqrt(1.+ecc)*np.tan(EEu/2.), np.sqrt(1.- ecc))
	

	#ww = arccos(cosw)

	
	### NOW GET RADIAL VEL FROM KEPLER PROBLEM
	#KK        = q/(1.+q) * nn*sep * sin(Inc)/sqrt(1.-e*e)   #/q for sec
	#vsec       = vmean + KK*( np.cos(ww + fot) + ecc*np.cos(ww) ) ##div by q to get seconday vel
	vRsec       =  vmean + KK*( cosw*np.cos(fot) + np.sqrt(1.-cosw*cosw)*np.sin(fot) + ecc*cosw ) ## cosw*np.cos(fot) +- sqrt(1. sample two diff halves of w space
	vRprm        = qq*vRsec 
		##A note on the rewriting of the above. I have used that cos(w+f) = cosw*cosf - sinw*sinf 
		## and then used sinw = sqrt(1-cos^2w). But sinw has range -1,1 while sqrt(1-cos^2w) has 
		## range 0,1. So we have to use +- the second term to get both top and bottom halves
		## of the trigonometric plane (cosw works the same on both halves).
				

	
	### NOW COMPUTE REL BEAMING FORM RAD VEL
	r_sec      = semimaj * (1. - ecc*ecc) / (1. + ecc* np.cos( fot ))
	vorbS      = 1./(1.+qq) * np.sqrt(G*Mbin *( 2./r_sec - 1./semimaj ))  #1./(1.+q) * sqrt(G*M/a)
	vorbP	   = qq*vorbS
	GamS       = 1./np.sqrt(1. - vorbS*vorbS/c/c)
	GamP       = 1./np.sqrt( 1. - vorbP*vorbP/c/c)
	
	DopS      = (1./(GamS * (1. - vRsec/c)) )**(3.-alp)
	DopP      = (1./(GamP * (1. + vRprm/c)) )**(3.-alp)
	
	DopLum	  = DopS*fs + DopP*(1.-fs)
	
	
	#mags      = -5./2. * np.log10(DopLum)  ## mag - mag0= -2.5 * log10(F(t)/F_0)  =  -2.5 * log10(DopLum) 
	

	# ## put limits on ecc
	# if (ecc<0. or ecc>=1.0):
	# 	DopLum = 99999999999.0*(t+1.)/(t+1.)
	# 	mags = 99999999999.0*(t+1.)/(t+1.)
	# if (KK/c>=1. or KK/c <0.):
	# 	DopLum = 99999999999.0*(t+1.)/(t+1.)
	# 	mags = 99999999999.0*(t+1.)/(t+1.)
	# if (cosw>1. or cosw < -1):
	# 	DopLum = -99999999999.0*(t+1.)/(t+1.)
	# 	mags = -99999999999.0*(t+1.)/(t+1.)
	# if (fs<0. or fs > 1):
	# 	DopLum = -99999999999.0*(t+1.)/(t+1.)
	# 	mags = -99999999999.0*(t+1.)/(t+1.)



	return DopLum #return fluxes here because data is countspersec ###mags






def DopLumLens(params, t):
	#Mbin, qq, ecc, cosw, T0, CosI, Tbin, fs, alp = params
	Mbin, qq, ecc, cosw, T0, CosI, Tbin, fs = params
	alp = 1.88
	#circ	
	#plens = [Mbin, Tbin, qq, CosI, T0*2.*np.pi/Tbin]
	#ecc
	plens = [Mbin, Tbin, qq, CosI, T0*2.*np.pi/Tbin, ecc, cosw]



	vmean=0.0
	# fs = 1.
	Mbin = 10.**Mbin*Msun
	Tbin = Tbin*day2sec  #(put in sec)
	T0 = T0*day2sec #- 0.5*Tbin
	t = t*day2sec
	#KK = KK*c       # put in cm/s

	#zz = in rest frmame alrdy
	
	semimaj = (np.sqrt(G*Mbin)*Tbin/2./np.pi)**(2./3.)

	##Check this
	KK = 1./(1.+qq) * np.sqrt(G*Mbin/semimaj)  *  CosI ## I from faceon here but from edge on in lens calc!

	
	#### Inc measured from FACE ON

	### Solve for Eccentric Anamoly and then f(t)
	EEu = np.mod(EAn(t,ecc,T0,Tbin), 2.*np.pi)
	fot = 2.*np.arctan2(np.sqrt(1.+ecc)*np.tan(EEu/2.), np.sqrt(1.- ecc))
	




	#ww = arccos(cosw)

	ww = ma.acos(cosw)
	### NOW GET RADIAL VEL FROM KEPLER PROBLEM
	#KK        = q/(1.+q) * nn*sep * sin(Inc)/sqrt(1.-e*e)   #/q for sec
	vRsec       = vmean + KK*( np.cos(ww + fot) + ecc*np.cos(ww) ) ##div by q to get seconday vel
	#vRsec       =  vmean + KK*( cosw*np.cos(fot) - np.sqrt(1.-cosw*cosw)*np.sin(fot) + ecc*cosw ) ## cosw*np.cos(fot) +- sqrt(1. sample two diff halves of w space
	vRprm        = qq*vRsec 
			
	

	
	### NOW COMPUTE REL BEAMING FORM RAD VEL
	r_sec      = semimaj * (1. - ecc*ecc) / (1. + ecc* np.cos( fot + ww ))
	vorbS      = 1./(1.+qq) * np.sqrt(G*Mbin *( 2./r_sec - 1./semimaj ))  #1./(1.+q) * sqrt(G*M/a)
	vorbP	   = qq*vorbS
	GamS       = 1./np.sqrt(1. - vorbS*vorbS/c/c)
	GamP       = 1./np.sqrt(1. - vorbP*vorbP/c/c)
	
	DopS      = (1./(GamS * (1. - vRsec/c)) )**(3.-alp)
	DopP      = (1./(GamP * (1. + vRprm/c)) )**(3.-alp)


	###TAKE INTO ACCOUNT TIME DELAY
	t_rtd_sec = t - r_sec/c/(1.+qq)*(1.-np.cos( fot + np.pi/2.)) ##is this correct??
	t_rtd_prm = t - qq*r_sec/c*(1.-np.cos( fot + np.pi/2. ))



	DopLum	  = DopS*fs*MagPSsec_ecc(plens, t_rtd_sec/day2sec) + DopP*(1.-fs)*MagPSprm_ecc(plens, t_rtd_prm/day2sec)

	#DopLum	  = DopS*fs*MagPSsec_ecc(plens, t/day2sec) + DopP*(1.-fs)*MagPSprm_ecc(plens, t/day2sec)
	## if fs=0.5 and q=1, this goes to zero in lim of no lensing and small v, in practice, when v not very small, get small double period modulation
	
	#mags      = -5./2. * np.log10(DopLum)  ## mag - mag0= -2.5 * log10(F(t)/F_0)  =  -2.5 * log10(DopLum) 
	

	# ## put limits on ecc
	# if (ecc<0. or ecc>=1.0):
	# 	DopLum = 99999999999.0*(t+1.)/(t+1.)
	# 	mags = 99999999999.0*(t+1.)/(t+1.)
	# if (KK/c>=1. or KK/c <0.):
	# 	DopLum = 99999999999.0*(t+1.)/(t+1.)
	# 	mags = 99999999999.0*(t+1.)/(t+1.)
	# if (cosw>1. or cosw < -1):
	# 	DopLum = -99999999999.0*(t+1.)/(t+1.)
	# 	mags = -99999999999.0*(t+1.)/(t+1.)
	# if (fs<0. or fs > 1):
	# 	DopLum = -99999999999.0*(t+1.)/(t+1.)
	# 	mags = -99999999999.0*(t+1.)/(t+1.)



	return DopLum #return fluxes here because data is countspersec ###mags



