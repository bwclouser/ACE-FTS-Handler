;+
; NAME:
;       ace_benchmark
;
; PURPOSE:
;       Times ACEData::calcBlocks on synthetic data shaped like the allMerid
;       preset: 17 latitude bands x 1 longitude x 35 altitude levels x 1 time
;       bin, which is 595 blocks.
;
;       The checksums are printed so a change in timing can be shown not to
;       have changed the results. They are sums over all blocks, so they are
;       sensitive to any block whose statistics move.
;
; CALLING SEQUENCE:
;
;       ace_benchmark
;       ace_benchmark,npass=20000
;
; KEYWORD PARAMETERS:
;
;       npass     occultations to generate. Default 4000. ACE-FTS itself has
;                 tens of thousands per molecule, so the default is modest.
;
; MODIFICATION HISTORY:
;       Written by: Ben Clouser, 8/4/25
;-

PRO ace_benchmark,npass=npass

  COMPILE_OPT IDL2

  IF N_ELEMENTS(npass) EQ 0 THEN npass=4000

  d=OBJ_NEW('ACETest')
  d.injectHmz,npass=npass,nmols=2,/global

  d.setParam,'latlims',[-85,85]
  d.setParam,'lonlims',[-180,180]
  d.setParam,'altlims',[5,40]
  d.setParam,'nLats',17
  d.setParam,'nLons',1
  d.setParam,'nAlts',35
  d.setParam,'yearlims',[2004,2023]
  d.setParam,'monthlims',[1,12]
  d.setParam,'daylims',[1,31]
  d.setParam,'hourlims',[1d-4,24d0-1d-4]
  d.setParam,'nTimes',1
  d.setParam,'elim',1d-3
  d.setParam,'cuttype','pointwise'

  t0=SYSTIME(1)
  d.calcBlocks
  t1=SYSTIME(1)

  b=d.getBlocks()

  PRINT,'occultations     : ',STRTRIM(npass,2)
  PRINT,'blocks           : ',STRTRIM(N_ELEMENTS(b),2)
  PRINT,'points binned    : ',STRTRIM(LONG(TOTAL(b.nPts)),2)
  PRINT,'checksum mean0   : ',TOTAL(b.mean0,/NAN)
  PRINT,'checksum meanRat : ',TOTAL(b.meanRat,/NAN)
  PRINT,'checksum meanLat : ',TOTAL(b.meanLat,/NAN)
  PRINT,'calcBlocks       : ',STRING(t1-t0,FORMAT='(F8.2)')+' s'

  OBJ_DESTROY,d

END
