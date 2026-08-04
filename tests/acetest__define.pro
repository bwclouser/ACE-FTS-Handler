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
;         mixrat0[a,p]     = (a + 1 + p/10) * 1e-6
;         mixrat1[a,p]     = 2 * mixrat0[a,p]      (so every ratio is exactly 2)
;         temp[a,p]        = 200 + a/2 + p/10
;         pres[a,p]        = exp(-a/7)
;
;       Latitude and longitude both advance with the pass index, so a lat
;       cut and a lon cut select the same passes and the expected membership
;       of a block stays easy to reason about.
;
; KEYWORD PARAMETERS:
;
;       npass     number of occultations to generate. Default 40.
;
;       nmols     1 or 2. With 2, mixrat1 and ratio are filled in. Default 1.
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

  ; mixrat0[a,p] = (a + 1 + p/10) * 1e-6, built by outer sum.
  mixrat0=((REBIN(a,nlev,npass)+1.0)+REBIN(TRANSPOSE(p),nlev,npass)/10.0)*1.0e-6
  mixrat0_err=FLTARR(nlev,npass)+1.0e-9

  temp=200.0+REBIN(a,nlev,npass)/2.0+REBIN(TRANSPOSE(p),nlev,npass)/10.0
  pres=EXP(-REBIN(a,nlev,npass)/7.0)

  IF N_ELEMENTS(badErr) EQ 3 THEN mixrat0_err[badErr[0],badErr[1]:badErr[2]]=1.0
  IF N_ELEMENTS(fillVal) EQ 3 THEN mixrat0[fillVal[0],fillVal[1]:fillVal[2]]=-999.0

  h={nmols:nmols,name0:'H2O',npass:npass,alt:alt, $
     srss:INTARR(npass),orbit:LINDGEN(npass)+1000L, $
     year:year,month:month,day:day,hour:hour, $
     juldays:JULDAY(month,day,year,hour,0d0,0d0),lat:lat,lon:lon, $
     mixrat0:mixrat0,mixrat0_err:mixrat0_err,temp:temp,pres:pres, $
     pt:temp*(1d0/pres)^.286d0,t_dmp:pres*!VALUES.F_NAN}

  IF nmols GT 1 THEN BEGIN
    mixrat1=2.0*mixrat0
    h=CREATE_STRUCT(h,'name1','HDO','mixrat1',mixrat1, $
                    'mixrat1_err',FLTARR(nlev,npass)+1.0e-9, $
                    'ratio',mixrat1/mixrat0)
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
;       npass     occultations per molecule. Default 6.
;       nDMP      how many of those orbits the DMP covers. Default 3.
;-

PRO ACETest::injectRaw,npass=npass,nDMP=nDMP

  COMPILE_OPT IDL2

  IF N_ELEMENTS(npass) EQ 0 THEN npass=6
  IF N_ELEMENTS(nDMP) EQ 0 THEN nDMP=3

  nlev=150
  L=LIST()

  FOREACH nm,['H2O','HDO'] DO BEGIN
    L.add,{name:nm,npass:npass,alt:FINDGEN(nlev)+0.5, $
           srss:INTARR(npass),orbit:LINDGEN(npass)+100L, $
           year:INTARR(npass)+2010,month:INTARR(npass)+6, $
           day:INDGEN(npass)+1,hour:FLTARR(npass)+12.0, $
           lat:FLTARR(npass)+15.0,lon:FLTARR(npass)+105.0, $
           mixrat:FLTARR(nlev,npass)+5.0e-6,mixrat_err:FLTARR(nlev,npass)+1.0e-9, $
           temp:FLTARR(nlev,npass)+220.0,pres:FLTARR(nlev,npass)+0.1}
  ENDFOREACH

  (*self.raw)=L

  ; Theta and temperature are constant so the tests only have to check which
  ; columns are finite, not what value they hold.
  (*self.dmp)={srss:INTARR(nDMP),orbit:LINDGEN(nDMP)+100L, $
               theta:FLTARR(75,nDMP)+400.0,temp:FLTARR(75,nDMP)+215.0}

END

PRO ACETest__define

  COMPILE_OPT IDL2

  struct={ACETest,INHERITS ACEData}

END
