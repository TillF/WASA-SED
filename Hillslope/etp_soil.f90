SUBROUTINE etp_soil(i_subbas3,vegi,defi,facwi,act,actveg,acts,  &
        height_act,lai_act,alb_act,rsfinal)

!Jose: an upper theshold for soil surface resistance "rss" was introduced. For very high values of "rss" a divi! jose miguel: for very low soil water content, the rss tends to infinity. This is not accepted by the computations, so the maximum of surface resistance of bare soil is limited to 1.0e6 (arbitrary huge number beyond maximum value obtained by Domingo et al. 1999) in order to prevent errors. See Domingo et al./Agricultural and Forest Meteorology 95 (1999) 76-77 
!2012-12-05
 
!Till: computationally irrelevant: minor changes to improve compiler compatibility
!2011-04-29

!Till: excluded domean, donight as logicals to module common_h
! excluded daily_delta_temp (daily temperature variation),hours_of_daylight
! 2008-02-05

!Till: added safety checks to prevent NaN values when LAI or plant height=0
! renamed muni to i_subbas3
!2006-08-16
!

! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:46:00

use climo_h
use common_h
use hymo_h
use params_h
use time_h

!** EVAPOTRANSPIRATION

!   Update 08/2000
!   Two-layer and two-source model gemaess Brenner & Incoln (1997).
!   modification: soil heat flux is taken to be different for
!   soil below shrubs and for bare soil (see also Domingo et al. (1999)).


! Ausgehend davon, das Temperatur, Globalstrahlung und Luftfeuchte als taegl.,
! Munizipien (Subbasins) abhaengige  Werte gegeben sind, also temp(id,imun), rad(id,imun),
! rhum(id,imun), wind(d,imun)

IMPLICIT NONE


INTEGER, INTENT(IN)                  :: i_subbas3
INTEGER, INTENT(IN)                  :: vegi
REAL, INTENT(IN)                         :: defi
REAL, INTENT(IN)                         :: facwi
REAL, INTENT(OUT)                        :: act
REAL, INTENT(OUT)                        :: actveg
REAL, INTENT(OUT)                        :: acts
REAL, INTENT(IN)                         :: height_act(nveg)
REAL, INTENT(IN)                         :: lai_act(nveg)
REAL, INTENT(IN)                         :: alb_act(nveg)
REAL, INTENT(OUT)                        :: rsfinal

REAL :: seconds_of_daylight
REAL :: tempval,f2
REAL :: dumv,dums,dumbs


! subscripts of current subbasin, landscape unit and vegetation unit


! real evapotranspiration of different components
! (vegetation, soil under canopy, bare soil,
!  total, total day, total night (mm))
REAL :: actbs, actd,actn

! vegetation characteristics of actual day





! soil moisture deficit in the less easily
! plant available part of nFK (Hough, 1997: 230)

! relative soil moisture content in the freely
! plant available part of nFK

! net radiation and extraterestrial radiation (W/m^2)
REAL :: rnetto, rnettos
! available energy for latent and sensible heat flux (W/m^2)
! according to (S&W (1985), p.843)
REAL :: anetto, anettos, anettop, anettobs
! soil heat flux (W/m²) (fraction of net radiation at soil surface)
! (Gstream is assumed to be euqal below shrub and in open soil)
! (Gstream=Gstreambs) (not clear in Brenner & Incoll)
REAL :: gstream, gstreambs
! fraction of net radiation to be transformed to ground heat flux
REAL :: gfracd,gfracn
! transformation parameter of ETP from latent heat flux (W/m**2)
! into water equivalent (mm)
! (day and night parameter, depending on daylight length)
REAL :: transd,transn
! extinction coefficient (for net radiation)
REAL :: ext
! canopy extinction coefficient for PAR
REAL :: extpar
! parameter for canopy PAR extinction calculation
REAL :: bq
! vapour pressure e (hPa), saturation vapour pressure es (hPa)
! and slope of the saturation vapour pressure curve s (hPa/K)
REAL :: es,s
! daily mean vapour pressure (hPa)
REAL :: emean
! water vapour pressure deficit at mean surface flow height zm
REAL :: d0
! day and night temperature
REAL :: tempd,tempn
! albedo
REAL :: alpha
! areal fraction of vegetation cover
REAL :: fcov
! Bedeckungsgrad, Bewoelkugsfaktor
REAL :: nn,f
! psychrometric constant  (hPa/K)
REAL :: gamma
! aerodynamic resistances (s/m)
REAL :: rap, ras, rabs, raa
! canopy /soil resistances (s/m)
REAL :: rs, rsp, rss, rsbs
! height of measurement of wind speed zr, roughness lengths (m)
REAL :: zr, z0s, z0
! mean canpoy air stream height (fixed=0.76*canpoy height)
REAL :: zm
! variable for z0 variation
REAL :: x
! displacement height
REAL :: dp
! eddy decay coefficient
REAL :: n
! Karman's constant
REAL :: k
! eddy diffusivity at canopy top
REAL :: kh
! friction velocity
REAL :: ustern
! air density (kg/m^3)
REAL :: roh
! heat capacity of air (J/(kg´*K)  (=/ parameter cp below !)
REAL :: heatc
! parameters Rx
REAL :: grs,grp,grbs,gra
! parameters Cx
REAL :: cs,cp,cbs
! intermediate evap parameters of different sources
REAL :: pmp,pms,pmbs



! length of day (seconds)
seconds_of_daylight=hours_of_daylight*60.*60.
! conversion factor for latent heat flux W/m2 -> mm
transd=seconds_of_daylight/2.45E6
transn=(86400.-seconds_of_daylight)/2.45E6


! check wind
! program is not defined for windspeed=0.
!average wind velocity (m/s) is set constant
wind(d,i_subbas3)=1.0

! measurement height of wind speed
zr=6.0
! if canopy height > measurement height
! wind speed in canopy height is assumed to be that of measurement height

  zr=max(zr,height_act(vegi))

! constant parameters for day and night
!      fcov=1.-0.7**lai_act(vegi)
fcov=1.
! !! fcov effects LAI (Bezug: total or canopy area !?)

k=0.41
z0s=0.01
zm=0.76*height_act(vegi)
if (lai_act(vegi)==0.) then !Till: without this problems occur with certain compiler settings
	x=0.
	dp=0.
else
	x=0.07*lai_act(vegi)
	dp=1.1*height_act(vegi)*LOG(1.+x**0.25)
end if

IF (x < 0.2) THEN
  z0=z0s+0.3*height_act(vegi)*x**0.5
ELSE
  z0=0.3*(height_act(vegi)*(1.-dp/height_act(vegi)))
END IF
ustern=k*wind(d,i_subbas3)/LOG((zr-dp)/z0)
kh=k*ustern*(height_act(vegi)-dp)

alpha=alb_act(vegi)
gamma=0.67
roh=1.18
heatc=1005.
n=2.5

! calculation of cloudiness needs only relative radiation
! (-> no raditation correction for daylight length required)
nn = rad(d,i_subbas3)/radex(d)/0.55-0.18/0.55
nn=MAX(0.0,nn)
nn=MIN(1.0,nn)
f  = 0.1+0.9*nn
!     f  = 0.2+0.8*nN
ext=0.5
extpar=0.5
bq=100.


! aerodynamic resistances
!      ras=height_act(vegi)*EXP(n)/(n*Kh)*
!     .    (EXP(-1.*n*z0/height_act(vegi))-
!     .     EXP(-1.*n*zm/height_act(vegi)))

ras=height_act(vegi)*EXP(n)/(n*kh)* (EXP(-1.*n*z0s/height_act(vegi))-  &
    EXP(-1.*n*zm/height_act(vegi)))

raa=(1./(k*ustern))* LOG((zr-dp)/(height_act(vegi)-dp))+  &
    (height_act(vegi)/(n*kh))* (EXP(n*(1.-zm/height_act(vegi)))-1.)

if (lai_act(vegi)==0.) then	!Till: this is a clumsy workaround to prevent errors when the LAI is 0 ii
	rap=1000000.0
else
	rap=25./(2.*lai_act(vegi))
end if

rabs=((LOG(zm/z0s))**2)/ ((k**2)*wind(d,i_subbas3)*  &
    LOG(zm/z0)/LOG(zr/z0))
rabs=rabs+(ras-rabs)*fcov
!      write(*,*) 'ras,rap,rabs,raa',ras,rap,rabs,raa

! surface resistance of bare soil
! jose miguel: for very low soil water content, the rss tends to infinity. This is not accepted by the computations, so the maximum of surface resistance of bare soil is limited to 1.0e6 (arbitrary huge number beyond maximum value obtained by Domingo et al. 1999) in order to prevent errors. See Domingo et al./Agricultural and Forest Meteorology 95 (1999) 76-77 
if (facwi<1.0e6) then
   rss=facwi
else
   rss=1.0e6
endif
!      rss=1.0e6
rsbs=rss
!      write(*,*) 'rss',rss

! mean daily vapor pressure emean
es = 6.11*EXP(17.62*temp(d,i_subbas3)/(243.12+temp(d,i_subbas3)))
emean  = rhum(d,i_subbas3)*es/100.


! if not daily mean values are used
IF (.NOT. domean) THEN
! ............................................................
! daytime values
  
  tempd=temp(d,i_subbas3)+daily_delta_temp
  es = 6.11*EXP(17.62*tempd/(243.12+tempd))
  s  = es*4284./((243.12+tempd)**2)
  gfracd=0.2
  
!  stomata/canopy resistance
!  given in data file is minimum leaf resistance based on unit LAI
!  for no stress conditions and light saturation
!  (according to Körner (1994))
  rs=resist(vegi)
  
! stomata resistance modification by soil water stress (defi) and
! water vapour pressure deficit (f2) (Hanan & Prince, 1997)
! (based on a multiplicative Jarvis-type model)
  f2=1./(1.+0.03*(es-emean))
  rs=1./((1./rs)*defi*f2)
  
!  mean leaf resistance in canopy depends on extinction of
!  photosnthetically active radiation within canopy
!  as function of LAI (Brenner & Lincoln, 1997, p.190).
!  (here parameterized according to Saugier & Katerji (1991) for
!   total short-wave radiation)
  rs=1./(((1./rs)/extpar)* LOG((bq+extpar*rad(d,i_subbas3)*24./hours_of_daylight)/  &
      (bq+extpar*rad(d,i_subbas3)*24./hours_of_daylight* EXP(-1.*3.1416*MAX(0.01,lai_act(vegi))))))
  !Till: added MAX() function to prevent NaN when LAI is 0

! canopy conductance is mean leaf conductance (per unit LAI)
! multiplied by LAI
! integration is already done by above equation
!      rsp=1./((1./rs)*lai_act(vegi))
  rsp=rs
  rsfinal=rs
  
  
  
  grs =(s+gamma)*ras +gamma*rss
  grp =(s+gamma)*rap +gamma*rsp
  grbs=(s+gamma)*rabs+gamma*rsbs
  gra =(s+gamma)*raa
  
  tempval=grs*grp*grbs+ (1.-fcov)*grs*grp*gra+  &
      fcov*grbs*grs*gra+ fcov*grbs*grp*gra
  cs =grbs*grp*(grs +gra)/tempval
  cp =grbs*grs*(grp +gra)/tempval
  cbs=grs *grp*(grbs+gra)/tempval
  
! daytime time mean of short wave radiation = daily mean * 24/hours_of_daylight
! Rnetto=Rkurz+Rlang
  
  rnetto = -f*(0.52-0.065*SQRT(emean))* 5.67E-8*(tempd+273.2)**4 +&
	   (1.-alpha)*rad(d,i_subbas3)*24./hours_of_daylight
  rnettos =rnetto*EXP(-1.*ext*lai_act(vegi))
!      Gstream =gfracd*Rnettos
!      Anetto  =Rnetto-Gstream
!      Anettos =Rnettos-Gstream
!      Anettobs=Rnetto-Gstream
!      Anettop =Anetto-Anettos
  
  gstream  =gfracd*rnettos
  gstreambs=gfracd*rnetto
  anetto   =rnetto-(gstream*fcov+gstreambs*(1.-fcov))
  anettos  =rnettos-gstream
  anettobs =rnetto-gstreambs
  anettop  =anetto-anettos
  
  
  pmp=(s*anetto+((roh*heatc*(es-emean)-s*rap*anettos)/ (raa+rap)))/  &
      (s+gamma*(1.+rsp/(raa+rap)))*transd
  pms=(s*anetto+((roh*heatc*(es-emean)-s*ras*anettop)/ (raa+ras)))/  &
      (s+gamma*(1.+rss/(raa+ras)))*transd
  pmbs=(s*anettobs+roh*heatc*(es-emean)/(raa+rabs))/  &
      (s+gamma*(1.+rsbs/(raa+rabs)))*transd
  
  actd=fcov*(cp*pmp+cs*pms)+(1.-fcov)*cbs*pmbs
  
! water vapour saturation deficit at exchange height
  d0=(es-emean)+(s*anetto-(s+gamma)*(actd/transd))* raa/(roh*heatc)
  
  
! evaporation of different sources (mm)
  actveg=fcov*transd* (s*anettop+roh*heatc*d0/rap)/  &
      (s+gamma*(1.+rsp/rap))
  acts  =fcov*transd* (s*anettos+roh*heatc*d0/ras)/  &
      (s+gamma*(1.+rss/ras))
  actbs =(1.-fcov)*transd* (s*anettobs+roh*heatc*d0/rabs)/  &
      (s+gamma*(1.+rsbs/rabs))
  
  act=actd
!      write(*,*) 'ETP,day, cov',actd,fcov
!      write(*,*) 'ETP,day,check',actveg+acts+actbs
!      write(*,*) 'ETP,day,veg,s,bs',actveg,acts,actbs
  
  
! ............................................................
  IF (donight) THEN
! night values
    
! canopy resistance as function of nFK deficit (Hough, 1997)
    rsp=2500.
    
    tempn=temp(d,i_subbas3)-daily_delta_temp
    es = 6.11*EXP(17.62*tempn/(243.12+tempn))
    s  = es*4284./((243.12+tempn)**2)
    gfracn=0.7
    
    grs =(s+gamma)*ras +gamma*rss
    grp =(s+gamma)*rap +gamma*rsp
    grbs=(s+gamma)*rabs+gamma*rsbs
    gra =(s+gamma)*raa
    
    tempval=grs*grp*grbs+ (1.-fcov)*grs*grp*gra+  &
        fcov*grbs*grs*gra+ fcov*grbs*grp*gra
    cs =grbs*grp*(grs +gra)/tempval
    cp =grbs*grs*(grp +gra)/tempval
    cbs=grs *grp*(grbs+gra)/tempval
    
! night short wave radiation = 0.
! Rnetto=Rlang
    
    rnetto = -f*(0.52-0.065*SQRT(emean))* 5.67E-8*(tempn+273.2)**4
    rnettos =rnetto*EXP(-1.*ext*lai_act(vegi))
!      Gstream =gfracn*Rnettos
!      Anetto  =Rnetto-Gstream
!      Anettos =Rnettos-Gstream
!      Anettobs=Rnetto-Gstream
!      Anettop =Anetto-Anettos
    
    gstream  =gfracn*rnettos
    gstreambs=gfracn*rnetto
    anetto   =rnetto-(gstream*fcov+gstreambs*(1.-fcov))
    anettos  =rnettos-gstream
    anettobs =rnetto-gstreambs
    anettop  =anetto-anettos
    
    
    pmp=(s*anetto+((roh*heatc*(es-emean)-s*rap*anettos)/ (raa+rap)))/  &
        (s+gamma*(1.+rsp/(raa+rap)))*transn
    pms=(s*anetto+((roh*heatc*(es-emean)-s*ras*anettop)/ (raa+ras)))/  &
        (s+gamma*(1.+rss/(raa+ras)))*transn
    pmbs=(s*anettobs+roh*heatc*(es-emean)/(raa+rabs))/  &
        (s+gamma*(1.+rsbs/(raa+rabs)))*transn
    
    actn=fcov*(cp*pmp+cs*pms)+(1.-fcov)*cbs*pmbs
    
! water vapour saturation deficit at exchange height
    d0=(es-emean)+(s*anetto-(s+gamma)*(actn/transn))* raa/(roh*heatc)
!      write(*,*) 'def',es-emean,d0
    
! evaporation of different sources (mm)
    dumv=actveg
    actveg=actveg+ fcov*transn*  &
        (s*anettop+roh*heatc*d0/rap)/ (s+gamma*(1.+rsp/rap))
    dums=acts
    acts  =acts+ fcov*transn*  &
        (s*anettos+roh*heatc*d0/ras)/ (s+gamma*(1.+rss/ras))
    dumbs=actbs
    actbs =actbs+ (1.-fcov)*transn*  &
        (s*anettobs+roh*heatc*d0/rabs)/ (s+gamma*(1.+rsbs/rabs))
    
    act=act+actn
!      write(*,*) 'ETP, night',actn
!      write(*,*) 'ETP total',act
!      write(*,*) 'check:',actveg+acts+actbs
!      write(*,*) 'ETP,night,veg,s,bs',actveg-dumv,acts-dums,actbs-dumbs
! End of night calculation
  END IF
  
  
! aggregate soil evaporation from bare soil and soil under canopy
  acts=acts+actbs
  
! .........................................................................
! End of day/night separate calculations, do daily mean calculations
ELSE
  
! length of day (seconds)
  seconds_of_daylight=hours_of_daylight*60.*60.
! conversion factor for latent heat flux W/m2 -> mm
  transd=seconds_of_daylight/2.45E6
  
  es = 6.11*EXP(17.62*temp(d,i_subbas3)/(243.12+temp(d,i_subbas3)))
  s  = es*4284./((243.12+temp(d,i_subbas3))**2)
  
! soil heat flux assumed to be negligable
  gfracd=0.
  
!  stomata/canopy resistance
!  given in data file is minimum leaf resistance based on unit LAI
!  for no stress conditions and light saturation
!  (according to Körner (1994)
  rs=resist(vegi)
  
! stomata resistance modification by soil water stress (defi) and
! water vapour pressure deficit (f2) (Hanan & Prince, 1997)
! (based on a multiplicative Jarvis-type model)
  f2=1./(1.+0.03*(es-emean))
  rs=1./((1./rs)*defi*f2)
  
!  mean leaf resistance in canopy depends on extinction of
!  photosnthetically active radiation within canopy
!  as function of LAI (Brenner & Lincoln, 1997, p.190).
!  (here parameterized according to Saugier & Katerji (1991) for
!   total short-wave radiation)
  rs=1./(((1./rs)/extpar)* LOG((bq+extpar*rad(d,i_subbas3)*24./hours_of_daylight)/  &
      (bq+extpar*rad(d,i_subbas3)*24./hours_of_daylight* EXP(-1.*3.1416*lai_act(vegi)))))
  
! canopy conductance is mean leaf conductance (per unit LAI)
! multiplied by LAI
! integration is already done by above equation
!      rsp=1./((1./rs)*lai_act(vegi))
  rsp=rs
  rsfinal=rs
  
  grs =(s+gamma)*ras +gamma*rss
  grp =(s+gamma)*rap +gamma*rsp
  grbs=(s+gamma)*rabs+gamma*rsbs
  gra =(s+gamma)*raa
  
  tempval=grs*grp*grbs+ (1.-fcov)*grs*grp*gra+  &
      fcov*grbs*grs*gra+ fcov*grbs*grp*gra
  cs =grbs*grp*(grs +gra)/tempval
  cp =grbs*grs*(grp +gra)/tempval
  cbs=grs *grp*(grbs+gra)/tempval
  
! daytime time mean of short wave radiation = daily mean * 24/hours_of_daylight
! Rnetto=Rkurz+Rlang
  
  rnetto = -f*(0.52-0.065*SQRT(emean))* 5.67E-8*(temp(d,i_subbas3)+273.2)**4 +&  
      (1.-alpha)*rad(d,i_subbas3)*24./hours_of_daylight
  rnettos =rnetto*EXP(-1.*ext*lai_act(vegi))
!      Gstream  =gfracd*Rnettos
!      Anetto   =Rnetto-Gstream
!      Anettos  =Rnettos-Gstream
!      Anettobs =Rnetto-Gstream
!      Anettop  =Anetto-Anettos
  
  gstream  =gfracd*rnettos
  gstreambs=gfracd*rnetto
  anetto   =rnetto-(gstream*fcov+gstreambs*(1.-fcov))
  anettos  =rnettos-gstream
  anettobs =rnetto-gstreambs
  anettop  =anetto-anettos
  
  
  pmp=(s*anetto+((roh*heatc*(es-emean)-s*rap*anettos)/ (raa+rap)))/  &
      (s+gamma*(1.+rsp/(raa+rap)))*transd
  pms=(s*anetto+((roh*heatc*(es-emean)-s*ras*anettop)/ (raa+ras)))/  &
      (s+gamma*(1.+rss/(raa+ras)))*transd
  pmbs=(s*anettobs+roh*heatc*(es-emean)/(raa+rabs))/  &
      (s+gamma*(1.+rsbs/(raa+rabs)))*transd
  
  actd=fcov*(cp*pmp+cs*pms)+(1.-fcov)*cbs*pmbs
  
! water vapour saturation deficit at exchange height
  d0=(es-emean)+(s*anetto-(s+gamma)*(actd/transd))* raa/(roh*heatc)
  
  
! evaporation of different sources (mm)
  actveg=fcov*transd* (s*anettop+roh*heatc*d0/rap)/  &
      (s+gamma*(1.+rsp/rap))
  acts  =fcov*transd* (s*anettos+roh*heatc*d0/ras)/  &
      (s+gamma*(1.+rss/ras))
  actbs =(1.-fcov)*transd* (s*anettobs+roh*heatc*d0/rabs)/  &
      (s+gamma*(1.+rsbs/rabs))
  
  act=actd
!      write(*,*) 'ETP,daymean, cov',actd,fcov
!      write(*,*) 'ETP,daymean,check',actveg+acts+actbs
  
! aggregate soil evaporation from bare soil and soil under canopy
  acts=acts+actbs
  
END IF


RETURN
END SUBROUTINE etp_soil
