;+
; NAME:
;       ace_run_tests
;
; PURPOSE:
;       Regression tests for the ACEData class, run against synthetic data
;       injected through the ACETest subclass.
;
;       Every expected value here is derived from the generating formulae
;       documented in ACETest::injectHmz, computed independently of the
;       binning code under test. So these are not a recorded baseline: they
;       assert what the answer should be, which means they stay meaningful
;       across refactors of calcBlocks, setOpMask and reduceBlock.
;
; CALLING SEQUENCE:
;
;       ace_run_tests
;
;   or, to make the process exit status reflect the result:
;
;       idl -e 'ace_run_tests,/batch'
;
; KEYWORD PARAMETERS:
;
;       batch     EXIT with status 1 if any test failed, 0 otherwise.
;
; NOTES:
;       nTimes is 1 throughout. The nTimes > 1 path is knowingly broken --
;       setOpMask masks on the global time flag rather than the per-bin
;       limits, so every time bin comes out identical -- and there is no
;       point asserting against behaviour that is already known to be wrong.
;       Add coverage here when that is fixed.
;
; MODIFICATION HISTORY:
;       Written by: Ben Clouser, 8/3/25
;-

PRO ace_test_init

  COMPILE_OPT IDL2
  COMMON ace_test_tally,nRun,nFail

  nRun=0L
  nFail=0L

END


PRO ace_test_report,name,ok,actual,expected

  COMPILE_OPT IDL2
  COMMON ace_test_tally,nRun,nFail

  nRun+=1
  IF ok THEN BEGIN
    PRINT,'  pass  '+name
  ENDIF ELSE BEGIN
    nFail+=1
    PRINT,'  FAIL  '+name
    IF N_PARAMS() GE 3 THEN PRINT,'          got      ',actual
    IF N_PARAMS() GE 4 THEN PRINT,'          expected ',expected
  ENDELSE

END


PRO ace_assert_eq,name,actual,expected

  COMPILE_OPT IDL2

  ok=(N_ELEMENTS(actual) EQ N_ELEMENTS(expected))
  IF ok THEN ok=(TOTAL(actual NE expected) EQ 0)
  ace_test_report,name,ok,actual,expected

END


PRO ace_assert_close,name,actual,expected,tol

  COMPILE_OPT IDL2

  IF N_ELEMENTS(tol) EQ 0 THEN tol=1.0d-5

  ; Relative where there is something to be relative to, absolute at zero.
  scale=ABS(expected) > 1.0d-30
  ok=(N_ELEMENTS(actual) EQ N_ELEMENTS(expected))
  IF ok THEN ok=(TOTAL(FINITE(actual) NE 1) EQ 0)
  IF ok THEN ok=(MAX(ABS(actual-expected)/scale) LE tol)
  ace_test_report,name,ok,actual,expected

END


PRO ace_assert_true,name,cond

  COMPILE_OPT IDL2

  ace_test_report,name,KEYWORD_SET(cond)

END


;----- Applies the standard synthetic grid: 2 latitude bands x 1 longitude
; x 2 altitude levels x 1 time bin, with every occultation inside the global
; spatial and temporal limits. -------;
PRO ace_test_setgrid,d,nLats=nLats,latlims=latlims,altlims=altlims,nAlts=nAlts

  COMPILE_OPT IDL2

  IF N_ELEMENTS(nLats) EQ 0 THEN nLats=2
  IF N_ELEMENTS(latlims) EQ 0 THEN latlims=[4.5,44.5]
  IF N_ELEMENTS(altlims) EQ 0 THEN altlims=[20,22]
  IF N_ELEMENTS(nAlts) EQ 0 THEN nAlts=2

  d.setParam,'latlims',latlims
  d.setParam,'lonlims',[99.5,139.5]
  d.setParam,'altlims',altlims
  d.setParam,'nLats',nLats
  d.setParam,'nLons',1
  d.setParam,'nAlts',nAlts
  d.setParam,'yearlims',[2010,2010]
  d.setParam,'monthlims',[5,7]
  d.setParam,'daylims',[1,31]
  d.setParam,'hourlims',[1d-4,24d0-1d-4]
  d.setParam,'nTimes',1
  d.setParam,'elim',1d-3
  d.setParam,'cuttype','pointwise'

END


PRO ace_run_tests,batch=batch

  COMPILE_OPT IDL2
  COMMON ace_test_tally,nRun,nFail

  ace_test_init

  ; ================================================================
  PRINT,'Block statistics on a populated 2 x 1 x 2 grid'
  ; ================================================================
  ; Each block holds 20 consecutive passes at a single altitude level.
  ; mixrat0[a,p] = (a + 1 + p/10) * 1e-6, so a block at level a covering
  ; passes p0..p0+19 has mean (a + 1 + mean(p0..p0+19)/10) * 1e-6.

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=40
  ace_test_setgrid,d
  d.calcBlocks
  b=d.getBlocks()

  ace_assert_eq,'nPts is 20 in every block',REFORM(b.nPts),[20,20,20,20]

  lo=FINDGEN(20)          ; passes 0..19  -> latitude band 0
  hi=FINDGEN(20)+20.0     ; passes 20..39 -> latitude band 1

  ace_assert_close,'mean0, lat band 0, level 20',b[0,0,0,0].mean0,MEAN((21.0+lo/10.0)*1.0e-6)
  ace_assert_close,'mean0, lat band 1, level 20',b[1,0,0,0].mean0,MEAN((21.0+hi/10.0)*1.0e-6)
  ace_assert_close,'mean0, lat band 0, level 21',b[0,0,1,0].mean0,MEAN((22.0+lo/10.0)*1.0e-6)
  ace_assert_close,'mean0, lat band 1, level 21',b[1,0,1,0].mean0,MEAN((22.0+hi/10.0)*1.0e-6)

  ace_assert_close,'sd0, lat band 0, level 20',b[0,0,0,0].sd0,STDDEV((21.0+lo/10.0)*1.0e-6)

  ; Median from the boxplot data, against the same 20 values.
  ace_assert_close,'median (bw0[2]), lat band 0, level 20',b[0,0,0,0].bw0[2], $
                   MEDIAN((21.0+lo/10.0)*1.0e-6,/EVEN),1.0d-4

  ace_assert_close,'meant, lat band 0, level 20',b[0,0,0,0].meant,MEAN(200.0+20.0/2.0+lo/10.0)
  ace_assert_close,'meanLat, lat band 0',b[0,0,0,0].meanLat,MEAN(5.0+lo)
  ace_assert_close,'meanLon, lat band 1',b[1,0,0,0].meanLon,MEAN(100.0+hi)

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'Empty blocks are NaN with nPts 0'
  ; ================================================================
  ; Latitude bands [-90,-45], [-45,0] and [45,90] contain no occultations;
  ; only [0,45] does.

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=40
  ace_test_setgrid,d,nLats=4,latlims=[-90,90],altlims=[20,21],nAlts=1
  d.calcBlocks
  b=d.getBlocks()

  ace_assert_eq,'nPts is 0 in the three empty bands',REFORM(b[[0,1,3],0,0,0].nPts),[0,0,0]
  ace_assert_eq,'nPts is 40 in the populated band',b[2,0,0,0].nPts,40
  ace_assert_true,'mean0 is NaN in the empty bands', $
                  TOTAL(FINITE(REFORM(b[[0,1,3],0,0,0].mean0))) EQ 0
  ace_assert_true,'sd0 is NaN in the empty bands', $
                  TOTAL(FINITE(REFORM(b[[0,1,3],0,0,0].sd0))) EQ 0
  ace_assert_true,'mean0 is finite in the populated band',FINITE(b[2,0,0,0].mean0)
  ace_assert_true,'orbIDs is empty in an empty band', $
                  N_ELEMENTS(*b[0,0,0,0].orbIDs) EQ 0

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'Error limit and fill value both drop points'
  ; ================================================================

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=40,badErr=[20,0,4],fillVal=[21,0,2]
  ace_test_setgrid,d
  d.calcBlocks
  b=d.getBlocks()

  ace_assert_eq,'5 points cut by eLim at level 20',b[0,0,0,0].nPts,15
  ace_assert_eq,'3 points cut by -999 fill at level 21',b[0,0,1,0].nPts,17
  ace_assert_eq,'the other latitude band is untouched',REFORM(b[1,0,*,0].nPts),[20,20]
  ace_assert_close,'mean0 reflects only the surviving points',b[0,0,0,0].mean0, $
                   MEAN((21.0+(FINDGEN(15)+5.0)/10.0)*1.0e-6)

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'Two-molecule statistics'
  ; ================================================================
  ; mixrat[*,*,1] is exactly 2 * mixrat[*,*,0], so every ratio is exactly 2
  ; and mean1 must be twice mean0. This is what catches reduceBlock reading
  ; the wrong molecule into the second slot.

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=40,nmols=2
  ace_test_setgrid,d
  d.calcBlocks
  b=d.getBlocks()

  ace_assert_close,'meanRat is exactly 2',REFORM(b.meanRat),DBLARR(4)+2.0
  ace_assert_close,'sdRat is 0',REFORM(b.sdRat),DBLARR(4),1.0d-6
  ace_assert_close,'mean1 is twice mean0',REFORM(b.mean1),2.0*REFORM(b.mean0)
  ace_assert_true,'mean1 differs from mean0',TOTAL(REFORM(b.mean1) EQ REFORM(b.mean0)) EQ 0

  ; The generalized arrays must agree with the legacy scalar aliases.
  ace_assert_close,'mean[0] matches the mean0 alias',REFORM(b.mean[0]),REFORM(b.mean0)
  ace_assert_close,'mean[1] matches the mean1 alias',REFORM(b.mean[1]),REFORM(b.mean1)
  ace_assert_close,'sd[1] matches the sd1 alias',REFORM(b.sd[1]),REFORM(b.sd1)
  ace_assert_close,'meanRatio[0] matches the meanRat alias', $
                   REFORM(b.meanRatio[0]),REFORM(b.meanRat)
  ace_assert_close,'bwRatio[2,0] matches the bwRat[2] alias', $
                   REFORM(b.bwRatio[2,0]),REFORM(b.bwRat[2])

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'Arbitrary molecule count'
  ; ================================================================
  ; Molecule m is (m+1) times the primary, so ratio m-1 is exactly m+1
  ; everywhere and every statistic has a known closed form.

  FOREACH nm,[1,2,3,5] DO BEGIN

    d=OBJ_NEW('ACETest')
    d.injectHmz,npass=40,nmols=nm
    ace_test_setgrid,d
    d.calcBlocks
    b=d.getBlocks()

    tag='nmols='+STRTRIM(nm,2)+': '
    ace_assert_eq,tag+'nmols recorded on the block',b[0].nmols,nm
    ace_assert_eq,tag+'mean has one entry per molecule', $
                  N_ELEMENTS(b[0].mean),nm
    ; IDL drops a trailing dimension of size 1, so bw reports [5] when nmols=1.
    ace_assert_eq,tag+'bw is 5 x nmols',SIZE(b[0].bw,/DIMENSIONS), $
                  (nm EQ 1) ? [5] : [5,nm]
    ace_assert_eq,tag+'names has one entry per molecule', $
                  N_ELEMENTS(b[0].names),nm

    lo=FINDGEN(20)
    base=MEAN((21.0+lo/10.0)*1.0e-6)

    FOR m=0,nm-1 DO $
      ace_assert_close,tag+'mean['+STRTRIM(m,2)+'] is '+STRTRIM(m+1,2)+'x the primary', $
                       b[0,0,0,0].mean[m],(m+1)*base

    IF nm GT 1 THEN BEGIN
      ace_assert_eq,tag+'meanRatio has nmols-1 entries', $
                    N_ELEMENTS(b[0].meanRatio),nm-1
      FOR m=0,nm-2 DO $
        ace_assert_close,tag+'meanRatio['+STRTRIM(m,2)+'] is exactly '+STRTRIM(m+2,2), $
                         REFORM(b.meanRatio[m]),DBLARR(4)+(m+2.0)
      ; The ratio is constant, so its spread should be zero to within float
      ; rounding. It is exactly zero only when m+2 is a power of two, since
      ; then (m+2)*x/x is exact; 3*x/x and 5*x/x are not.
      FOR m=0,nm-2 DO $
        ace_assert_true,tag+'sdRatio['+STRTRIM(m,2)+'] is negligible', $
                        MAX(ABS(REFORM(b.sdRatio[m]))) LT 1.0e-6*(m+2.0)
      ; Legacy aliases must still track the first molecule and first ratio.
      ace_assert_close,tag+'mean0 alias tracks mean[0]',REFORM(b.mean0),REFORM(b.mean[0])
      ace_assert_close,tag+'meanRat alias tracks meanRatio[0]', $
                       REFORM(b.meanRat),REFORM(b.meanRatio[0])
    ENDIF ELSE BEGIN
      ace_assert_true,tag+'no ratio fields for a single molecule', $
                      TOTAL(TAG_NAMES(b[0]) EQ 'MEANRATIO') EQ 0
    ENDELSE

    ; Empty bins stay NaN across every molecule and ratio.
    ace_test_setgrid,d,nLats=4,latlims=[-90,90],altlims=[20,21],nAlts=1
    d.calcBlocks
    b2=d.getBlocks()
    ace_assert_true,tag+'empty bins are NaN for every molecule', $
                    TOTAL(FINITE(b2[0,0,0,0].mean)) EQ 0
    IF nm GT 1 THEN ace_assert_true,tag+'empty bins are NaN for every ratio', $
                    TOTAL(FINITE(b2[0,0,0,0].meanRatio)) EQ 0

    OBJ_DESTROY,d

  ENDFOREACH

  ; ================================================================
  PRINT,'harmonize across more than two molecules'
  ; ================================================================
  ; Molecule 1 is stored in reverse orbit order and molecules after the first
  ; withhold their last 2 occultations, so the intersection has to cope with
  ; different sizes and different orderings.

  d=OBJ_NEW('ACETest')
  d.injectRaw,npass=8,nmols=4,drop=2,nDMP=3
  d.harmonize,[0,1,2,3]
  h=d.getHmzData()

  ace_assert_eq,'harmonize keeps the 6 common occultations',h.npass,6
  ace_assert_eq,'harmonize records 4 molecules',h.nmols,4
  ace_assert_eq,'mixrat is nlev x npass x nmols',SIZE(h.mixrat,/DIMENSIONS),[150,6,4]
  ace_assert_eq,'ratio is nlev x npass x nmols-1',SIZE(h.ratio,/DIMENSIONS),[150,6,3]
  ace_assert_eq,'orbits come out in primary order',h.orbit,LINDGEN(6)+100L
  ace_assert_true,'names are carried through', $
                  ARRAY_EQUAL(h.names,['H2O','HDO','H2O_18','H2O_17'])

  ; injectRaw makes molecule m equal to (m+1)*5e-6, so despite molecule 1 being
  ; reversed, a correct match gives ratio m-1 = m+2 at every point.
  FOR m=0,2 DO $
    ace_assert_close,'ratio['+STRTRIM(m,2)+'] is exactly '+STRTRIM(m+2,2), $
                     REFORM(h.ratio[0,*,m]),FLTARR(6)+(m+2.0)

  OBJ_DESTROY,d

  ; A molecule sharing no orbits with the primary must be rejected.
  d=OBJ_NEW('ACETest')
  d.injectRaw,npass=4,nmols=2,/disjoint
  d.harmonize,[0,1]
  ace_assert_true,'harmonize rejects a molecule with no common orbits', $
                  N_ELEMENTS(d.getHmzData()) EQ 0
  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'harmonize with a partial DMP'
  ; ================================================================
  ; The DMP covers 3 of 6 orbits. The other 3 must come back NaN rather
  ; than picking up the last DMP column, and levels above the 75-level DMP
  ; grid must be NaN rather than 0.

  d=OBJ_NEW('ACETest')
  d.injectRaw,npass=6,nDMP=3
  d.harmonize,[0,1],/addDMP
  h=d.getHmzData()

  ace_assert_eq,'all 6 orbits harmonized',h.npass,6
  ace_assert_close,'theta is 400 for the 3 matched orbits',REFORM(h.pt[0,0:2]),FLTARR(3)+400.0
  ace_assert_true,'theta is NaN for the 3 unmatched orbits', $
                  TOTAL(FINITE(REFORM(h.pt[0,3:5]))) EQ 0
  ace_assert_close,'DMP temperature is 215 for the matched orbits', $
                   REFORM(h.t_dmp[0,0:2]),FLTARR(3)+215.0
  ace_assert_true,'theta is NaN above the DMP grid', $
                  TOTAL(FINITE(REFORM(h.pt[100,*]))) EQ 0

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'getOcc'
  ; ================================================================

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=40
  o=d.getOcc([1000L,1005L],/SEARCH)
  ace_assert_eq,'/SEARCH finds 2 known orbits',o.npass,2
  ace_assert_close,'/SEARCH returns the right latitudes',o.lat,[5.0,10.0]

  ; Trailing dimension of 1 is dropped again: one molecule reports [150,2].
  ace_assert_eq,'getOcc mirrors the harmonized mixrat layout', $
                SIZE(o.mixrat,/DIMENSIONS),[150,2]

  o2=d.getOcc([999999L],/SEARCH)
  ace_assert_true,'/SEARCH on an absent orbit returns !NULL',N_ELEMENTS(o2) EQ 0

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'Parameter handling'
  ; ================================================================

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=40

  d.setParam,'nLats',5
  ace_assert_close,'setParam sets an exact name',d.getParam('nLats'),5.0
  d.setParam,'NLATS',6
  ace_assert_close,'setParam is case insensitive',d.getParam('nLats'),6.0

  ; 'nLat' is a substring of NLATS. The old substring match wrote to it.
  d.setParam,'nLat',99
  ace_assert_close,'setParam ignores a partial name',d.getParam('nLats'),6.0

  ; eLim used to be compared against the string 'elims', so the error stage
  ; was never marked dirty when the limit changed.
  d.setParam,'elim',1d-4
  ace_assert_close,'changing eLim marks the error stage dirty',d.getParam('eMod'),1.0

  SETENV,'ACE_DATA_ROOT=/tmp/aceroot'
  d.setParam,'dmp_path','v5'
  ace_assert_true,'the v5 alias resolves under ACE_DATA_ROOT', $
                  d.getParam('dmp_path') EQ $
                  '/tmp/aceroot/ACE-DMPSv5p2/ACEFTS_L2_v5p2_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
  d.setParam,'dmp_path','v4p1p2'
  ace_assert_true,'the v4p1p2 alias resolves under ACE_DATA_ROOT', $
                  d.getParam('dmp_path') EQ $
                  '/tmp/aceroot/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
  d.setParam,'dmp_path','/a/literal/path.nc4'
  ace_assert_true,'a literal DMP path passes through unchanged', $
                  d.getParam('dmp_path') EQ '/a/literal/path.nc4'
  SETENV,'ACE_DATA_ROOT='

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'Block array structure'
  ; ================================================================

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=40
  ace_test_setgrid,d
  d.calcBlocks
  b=d.getBlocks()
  tags=TAG_NAMES(b[0])

  ; areBlocksDefined and resetBlocks used to build different tag lists.
  ace_assert_true,'block struct carries the potential temperature fields', $
                  TOTAL(tags EQ 'MEANPT') EQ 1 AND TOTAL(tags EQ 'BWPT') EQ 1
  ace_assert_true,'block struct carries the DMP temperature fields', $
                  TOTAL(tags EQ 'MEANT_DMP') EQ 1 AND TOTAL(tags EQ 'BWT_DMP') EQ 1
  ace_assert_true,'block struct carries jullims', TOTAL(tags EQ 'JULLIMS') EQ 1

  ; nPts was a 16-bit int, which overflows at 32767 points in a bin.
  ace_assert_eq,'nPts is a 32-bit long',SIZE(b[0].nPts,/TYPE),3

  ; A live pointer in the REPLICATE template gave every block the same heap
  ; variable, so all of them reported the same orbit list.
  ace_assert_true,'blocks hold distinct orbID lists', $
                  TOTAL((*b[0,0,0,0].orbIDs) NE (*b[1,0,0,0].orbIDs)) GT 0

  OBJ_DESTROY,d

  ; ================================================================
  PRINT,'Prerequisite guards'
  ; ================================================================

  d=OBJ_NEW('ACETest')
  d.calcBlocks
  ace_assert_true,'calcBlocks without harmonized data returns quietly', $
                  N_ELEMENTS(d.getBlocks()) EQ 0
  OBJ_DESTROY,d

  d=OBJ_NEW('ACETest')
  d.harmonize,[0,1]
  ace_assert_true,'harmonize without raw data returns quietly', $
                  N_ELEMENTS(d.getHmzData()) EQ 0
  OBJ_DESTROY,d

  d=OBJ_NEW('ACETest')
  d.injectRaw,npass=6
  d.harmonize,[0,5]
  ace_assert_true,'harmonize rejects out-of-range elements', $
                  N_ELEMENTS(d.getHmzData()) EQ 0
  OBJ_DESTROY,d

  ; ================================================================
  PRINT,''
  PRINT,STRING(nRun-nFail,FORMAT='(I0)')+' / '+STRING(nRun,FORMAT='(I0)')+' passed'
  IF nFail GT 0 THEN PRINT,STRING(nFail,FORMAT='(I0)')+' FAILED'

  IF KEYWORD_SET(batch) THEN EXIT,STATUS=(nFail GT 0)

END
