SUBROUTINE etp_max(isubbas,vegi,etpmax,  &
        height_act,lai_act,alb_act,day,hh)

! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:45:55

use climo_h
use common_h
use params_h
use hymo_h
use time_h

!** POTENTIAL EVAPOTRANSPIRATION

!   10/2001
!   Penman-Monteith model, modified for day/night separate calculations
!   as in WASIM-ETH and in the actual evaporation routine in WASA

! Ausgehend davon, das Temperatur, Globalstrahlung und Luftfeuchte als taegl.,
! Munizipien/Subbasins abhaengige  Werte gegeben sind, also temp(id,imun), rad(id,imun),
! rhum(id,imun), wind(d,imun)
IMPLICIT NONE


INTEGER, INTENT(IN)                  :: isubbas
INTEGER, INTENT(IN)                  :: vegi
REAL, INTENT(OUT)                        :: etpmax
REAL, INTENT(IN)                         :: height_act(nveg)
REAL, INTENT(IN)                         :: lai_act(nveg)
REAL, INTENT(IN)                         :: alb_act(nveg)
INTEGER, INTENT(IN)                  :: day, hh


REAL,allocatable,save :: etpmax_stored(:,:) !Till: this serves to store values that have been calculated for a given subbasin and day. Later calls can be accelerated this way


REAL :: seconds_of_daylight
REAL :: tempval
REAL :: dumv,dums,dumbs


! subscripts of current subbasin, landscape unit and soil-vegetation unit


! real evapotranspiration of different components
! (vegetation, soil under canopy, bare soil,
!  total, total day, total night (mm))
REAL :: actveg,acts,actbs,actd,actn

! vegetation characteristics of actual day



! soil moisture deficit in the less easily
! plant available part of nFK (Hough, 1997: 230)
! REAL :: defi
! relative soil moisture content in the freely
! plant available part of nFK
! REAL :: facwi
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


if (.NOT. allocated(etpmax_stored)) then	!first call
	allocate(etpmax_stored(subasin,nveg))
	store_day=-1
	store_timestep=-1
end if

if (store_day/=day .OR. store_timestep/=hh) then	!new timestep begun - reset storage array
	etpmax_stored=-1.
	store_day=day
	store_timestep=hh
elseif (etpmax_stored(isubbas,vegi)/=-1.) then	!new timestep begun - reset storage array
	etpmax=etpmax_stored(isubbas,vegi)	!use saved value instead of recalculation
	return
end if





! length of day (seconds)
seconds_of_daylight=hours_of_daylight*60.*60.
! conversion factor for latent heat flux W/m2 -> mm
transd=seconds_of_daylight/2.45E6
transn=(86400.-seconds_of_daylight)/2.45E6



! check wind
! program is not defined for windspeed=0.
!wind velocity (m/s) is set constant for all sub-basins
!wind(d,isubbas)=1.0

! measurement height of wind speed
zr=6.0
! if canopy height > measurement height
! wind speed in canopy height is assumed to be that of measurement height
IF (zr < height_act(vegi)) THEN
  zr=height_act(vegi)
END IF


! constant parameters for day and night
!      fcov=1.-0.7**lai_act(vegi)
fcov=1.
! !! fcov effects LAI (Bezug: total or canopy area !?)

x=0.07*lai_act(vegi)
k=0.41
z0s=0.01
zm=0.76*height_act(vegi)
dp=1.1*height_act(vegi)*LOG(1.+x**0.25)
IF (x < 0.2) THEN
  z0=z0s+0.3*height_act(vegi)*x**0.5
ELSE
  z0=0.3*(height_act(vegi)*(1.-dp/height_act(vegi)))
END IF
ustern=k*wind(d,isubbas)/LOG((zr-dp)/z0)
kh=k*ustern*(height_act(vegi)-dp)

alpha=alb_act(vegi)
gamma=0.67
roh=1.18
heatc=1005.
n=2.5

! calculation of cloudiness needs only relative radiation
! (-> no raditation correction for daylight length required)
nn = rad(d,isubbas)/radex(d)/0.55-0.18/0.55
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

if (lai_act(vegi)/=0) then
    rap=25./(2.*lai_act(vegi))
else
    rap=1e10   !Till: no leaves, huge resistance
end if


rabs=((LOG(zm/z0s))**2)/ ((k**2)*wind(d,isubbas)*  &
    LOG(zm/z0)/LOG(zr/z0))
rabs=rabs+(ras-rabs)*fcov
!      write(*,*) 'ras,rap,rabs,raa',ras,rap,rabs,raa

! mimimum surface resistance of bare soil
rss=0.0
rsbs=rss
!      write(*,*) 'rss',rss

! mean daily vapor pressure emean
es = 6.11*EXP(17.62*temp(d,isubbas)/(243.12+temp(d,isubbas)))
emean  = rhum(d,isubbas)*es/100.

! minimum canopy resistance
    rsp=0.


! if not daily mean values are used
IF (domean .eqv. .FALSE.) THEN
! ............................................................

! daytime values

  tempd=temp(d,isubbas)+daily_delta_temp
  es = 6.11*EXP(17.62*tempd/(243.12+tempd))
  s  = es*4284./((243.12+tempd)**2)
  gfracd=0.2

! minimum surface resistance
  rs=0.

  grs =(s+gamma)*ras +gamma*rss
  grp =(s+gamma)*rap +gamma*rsp
  grbs=(s+gamma)*rabs+gamma*rsbs
  gra =(s+gamma)*raa

  tempval=grs*grp*grbs+ (1.-fcov)*grs*grp*gra+&
	    fcov*grbs*grs*gra+ fcov*grbs*grp*gra
  cs =grbs*grp*(grs +gra)/tempval
  cp =grbs*grs*(grp +gra)/tempval
  cbs=grs *grp*(grbs+gra)/tempval

! daytime time mean of short wave radiation = daily mean * 24/hours_of_daylight
! Rnetto=Rkurz+Rlang

  rnetto = -f*(0.52-0.065*SQRT(emean))* 5.67E-8*(tempd+273.2)**4&
	+(1.-alpha)*rad(d,isubbas)*24./hours_of_daylight
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

  !etpmax=actd
  etpmax=MAX(0.,actd)	!Till: prevented actual transpiration during day from getting negative

!      write(*,*) 'ETP,day, cov',actd,fcov
!      write(*,*) 'ETP,day,check',actveg+acts+actbs
!      write(*,*) 'ETP,day,veg,s,bs',actveg,acts,actbs


! ............................................................
  IF (donight) THEN
! night values


    tempn=temp(d,isubbas)-daily_delta_temp
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

    etpmax=etpmax+MAX(0.,actn)	!Till: negative values are discarded

!      write(*,*) 'ETP, night',actn
!      write(*,*) 'ETP total',act
!      write(*,*) 'check:',actveg+acts+actbs
!      write(*,*) 'ETP,night,veg,s,bs',actveg-dumv,acts-dums,actbs-dumbs
! End of night calculation
  END IF

! .........................................................................
! End of day/night separate calculations, do daily mean calculations
ELSE

! length of day (seconds)
  seconds_of_daylight=hours_of_daylight*60.*60.
! conversion factor for latent heat flux W/m2 -> mm
  transd=seconds_of_daylight/2.45E6

  es = 6.11*EXP(17.62*temp(d,isubbas)/(243.12+temp(d,isubbas)))
  s  = es*4284./((243.12+temp(d,isubbas))**2)

! soil heat flux assumed to be negligable
  gfracd=0.

! minimum surface resistiance
  rs=0.

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

  rnetto = -f*(0.52-0.065*SQRT(emean))* 5.67E-8*(temp(d,isubbas)+273.2)**4 +&
       (1.-alpha)*rad(d,isubbas)*24./hours_of_daylight
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

   !etpmax=actd
  etpmax=MAX(0.,actd)	!Till: prevented actual transpiration during day from getting negative

!      write(*,*) 'ETP,daymean, cov',actd,fcov
!      write(*,*) 'ETP,daymean,check',actveg+acts+actbs

END IF

etpmax_stored(isubbas,vegi)=etpmax	!save value for possible reuse

RETURN
END SUBROUTINE etp_max
