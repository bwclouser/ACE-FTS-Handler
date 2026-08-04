;+
; CLASS_NAME:
;       ACETest
;
; PURPOSE:
;       A test-only subclass of ACEData. Subclassing is what gives the tests
;       access to the protected instance data (self.hmz, self.raw, self.dmp),
;       so synthetic data can be injected without needing real ACE-FTS NetCDF
;       files on disk.
;
;       Everything injected here is deterministic and generated from a closed
;       form, which is what lets ace_run_tests state an exact expected value
;       for each statistic instead of comparing against a recorded baseline.
;       That matters because these tests exist mainly to guard refactors of
;       the binning core: an independently computed expectation catches a
;       change in results, where a recorded baseline would only catch a change
;       from whatever the code happened to do when the baseline was taken.
;
; SUPERCLASSES:
;       ACEData
;
; METHODS:
;       ACETest::injectHmz      -- synthetic harmonized data
;       ACETest::injectRaw      -- synthetic raw + DMP data, for harmonize
;
; MODIFICATION HISTORY:
;       Written by: Ben Clouser, 8/3/25
;-

FUNCTION ACETest::Init

  COMPILE_OPT IDL2

  RETURN,self.ACEData::Init()

END

;+
; =============================================================
;
; METHODNAME:
;       ACETest::injectHmz
;
; PURPOSE:
;       Populates self.hmz directly with a synthetic harmonized data set,
;       bypassing loadFiles and harmonize.
;
;       The generating formulae are:
;
;         lat[p]           = 5 + p                 (5 .. 5+npass-1)
;         lon[p]           = 100 + p               (100 .. 100+npass-1)
;         mixrat[a,p,0]    = (a + 1 + p/10) * 1e-6
;         mixrat[a,p,m]    = (m+1) * mixrat[a,p,0] for m > 0
;         temp[a,p]        = 200 + a/2 + p/10
;         pres[a,p]        = exp(-a/7)
;
;       So ratio[a,p,m-1] is exactly m+1 for every point, which makes the ratio
;       statistics of a multi-molecule set trivially checkable.
;
;       Latitude and longitude both advance with the pass index, so a lat
;       cut and a lon cut select the same passes and the expected membership
;       of a block stays easy to reason about.
;
; KEYWORD PARAMETERS:
;
;       npass     number of occultations to generate. Default 40.
;
;       nmols     how many molecules to generate. Any number >= 1. Default 1.
;
;       badErr    [level, firstPass, lastPass]. Sets mixrat0_err above any
;                 sane eLim for that inclusive pass range at that level, so
;                 the error flag drops exactly those points.
;
;       fillVal   [level, firstPass, lastPass]. Sets mixrat0 to the -999
;                 fill value over that inclusive pass range at that level.
;
;       global    Spread latitude over [-85,85], longitude over [-180,180] and
;                 the dates over 2004-2023, instead of the one-degree-per-pass
;                 march and the single June 2010 date. Intended for timing
;                 work, where a realistic pass count is needed and the exact
;                 membership of a block does not matter.
;-

PRO ACETest::injectHmz,npass=npass,nmols=nmols,badErr=badErr,fillVal=fillVal, $
                       global=global

  COMPILE_OPT IDL2

  IF N_ELEMENTS(npass) EQ 0 THEN npass=40
  IF N_ELEMENTS(nmols) EQ 0 THEN nmols=1

  nlev=150
  alt=FINDGEN(nlev)+0.5

  p=FINDGEN(npass)
  a=FINDGEN(nlev)

  IF KEYWORD_SET(global) THEN BEGIN

    lat=-85.0+170.0*p/(npass-1)
    ; A short cycle in longitude so successive passes are not colinear in
    ; lat/lon, which would put every occultation in one diagonal of bins.
    lon=-180.0+360.0*(p MOD 97)/96.0

    jd0=JULDAY(1,1,2004,0,0,0)
    jd1=JULDAY(12,31,2023,0,0,0)
    CALDAT,jd0+(jd1-jd0)*p/(npass-1),month,day,year,hour,mn,sc

  ENDIF ELSE BEGIN

    lat=5.0+p
    lon=100.0+p

    ; June 2010 for every pass, so the global time flag includes all of them.
    year=INTARR(npass)+2010
    month=INTARR(npass)+6
    day=FIX(p) MOD 28 + 1
    hour=FLTARR(npass)+12.0

  ENDELSE

  ; mixrat[a,p,0] = (a + 1 + p/10) * 1e-6, built by outer sum.
  mixrat0=((REBIN(a,nlev,npass)+1.0)+REBIN(TRANSPOSE(p),nlev,npass)/10.0)*1.0e-6

  temp=200.0+REBIN(a,nlev,npass)/2.0+REBIN(TRANSPOSE(p),nlev,npass)/10.0
  pres=EXP(-REBIN(a,nlev,npass)/7.0)

  ;----- Molecule m is (m+1) times the primary, so ratio m-1 is exactly m+1. --;
  mixrat=FLTARR(nlev,npass,nmols)
  mixrat_err=FLTARR(nlev,npass,nmols)+1.0e-9
  names=STRARR(nmols)
  molnames=['H2O','HDO','H2O_18','H2O_17','CH4','N2O','CO']

  FOR m=0,nmols-1 DO BEGIN
    mixrat[*,*,m]=(m+1.0)*mixrat0
    names[m]=(m LT N_ELEMENTS(molnames)) ? molnames[m] : 'MOL'+STRTRIM(m,2)
  ENDFOR

  ;----- badErr and fillVal act on the primary molecule. -------;
  IF N_ELEMENTS(badErr) EQ 3 THEN mixrat_err[badErr[0],badErr[1]:badErr[2],0]=1.0
  IF N_ELEMENTS(fillVal) EQ 3 THEN mixrat[fillVal[0],fillVal[1]:fillVal[2],0]=-999.0

  h={nmols:nmols,names:names,npass:npass,alt:alt, $
     srss:INTARR(npass),orbit:LINDGEN(npass)+1000L, $
     year:year,month:month,day:day,hour:hour, $
     juldays:JULDAY(month,day,year,hour,0d0,0d0),lat:lat,lon:lon, $
     mixrat:mixrat,mixrat_err:mixrat_err,temp:temp,pres:pres, $
     pt:temp*(1d0/pres)^.286d0,t_dmp:pres*!VALUES.F_NAN}

  IF nmols GT 1 THEN BEGIN
    ratio=FLTARR(nlev,npass,nmols-1)
    FOR m=1,nmols-1 DO ratio[*,*,m-1]=mixrat[*,*,m]/mixrat[*,*,0]
    h=CREATE_STRUCT(h,'ratio',ratio)
  ENDIF

  (*self.hmz)=h

END

;+
; =============================================================
;
; METHODNAME:
;       ACETest::injectRaw
;
; PURPOSE:
;       Populates self.raw with two synthetic molecules and self.dmp with a
;       DMP that deliberately covers only the first nDMP orbits, so that
;       harmonize has to deal with occultations that have no DMP match.
;
; KEYWORD PARAMETERS:
;
;       npass     occultations for the FIRST molecule. Default 6.
;       nDMP      how many of those orbits the DMP covers. Default 3.
;       nmols     how many molecules to generate. Default 2.
;       drop      number of trailing occultations to withhold from every
;                 molecule after the first, so harmonize has to intersect
;                 sets of different sizes. Default 0.
;       disjoint  shift the orbit numbers of every molecule after the first so
;                 that they share none with the primary.
;-

PRO ACETest::injectRaw,npass=npass,nDMP=nDMP,nmols=nmols,drop=drop,disjoint=disjoint

  COMPILE_OPT IDL2

  IF N_ELEMENTS(npass) EQ 0 THEN npass=6
  IF N_ELEMENTS(nDMP) EQ 0 THEN nDMP=3
  IF N_ELEMENTS(nmols) EQ 0 THEN nmols=2
  IF N_ELEMENTS(drop) EQ 0 THEN drop=0

  nlev=150
  L=LIST()
  molnames=['H2O','HDO','H2O_18','H2O_17','CH4','N2O','CO']

  FOR im=0,nmols-1 DO BEGIN
    nm=(im LT N_ELEMENTS(molnames)) ? molnames[im] : 'MOL'+STRTRIM(im,2)
    np=(im EQ 0) ? npass : npass-drop
    ;----- Reverse the orbit order for molecule 1 so that harmonize cannot rely
    ; on the molecules being stored in the same sequence. -------;
    orb=LINDGEN(np)+100L
    IF im EQ 1 THEN orb=REVERSE(orb)
    IF KEYWORD_SET(disjoint) && im GT 0 THEN orb+=1000L

    L.add,{name:nm,npass:np,alt:FINDGEN(nlev)+0.5, $
           srss:INTARR(np),orbit:orb, $
           year:INTARR(np)+2010,month:INTARR(np)+6, $
           day:INDGEN(np)+1,hour:FLTARR(np)+12.0, $
           lat:FLTARR(np)+15.0,lon:FLTARR(np)+105.0, $
           mixrat:FLTARR(nlev,np)+(im+1.0)*5.0e-6, $
           mixrat_err:FLTARR(nlev,np)+1.0e-9, $
           temp:FLTARR(nlev,np)+220.0,pres:FLTARR(nlev,np)+0.1}
  ENDFOR

  (*self.raw)=L

  ;----- nDMP=0 means "no DMP at all", so leave the pointer unassigned rather
  ; than trying to build zero-length arrays. -------;
  IF nDMP LE 0 THEN RETURN

  ; Theta and temperature are constant so the tests only have to check which
  ; columns are finite, not what value they hold.
  (*self.dmp)={srss:INTARR(nDMP),orbit:LINDGEN(nDMP)+100L, $
               theta:FLTARR(75,nDMP)+400.0,temp:FLTARR(75,nDMP)+215.0}

END

PRO ACETest__define

  COMPILE_OPT IDL2

  struct={ACETest,INHERITS ACEData}

END
