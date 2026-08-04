;+
; CLASS_NAME:
;       ACEData
;
; PURPOSE:
;       An ACEData object provides a simple interface with which to 
;       load, process, and display data from the ACE-FTS experiment.
;
; CATEGORY:
;       Data handling and plotting.
;
; SUPERCLASSES:
;       This class inherits from no other classes.
;
; SUBCLASSES:
;       This class has no subclasses.
;
; CREATION:
;       See ACEData::Init
;
; METHODS:
;       Intrinsic Methods
;       This class has the following methods:
;       
;       ACEData::Cleanup
;       ACEData::freeRaw
;       ACEData::loadFiles
;       ACEData::loadDMP
;       ACEData::harmonize
;       ACEData::setParam
;       ACEData::getParam
;       ACEData::getOcc
;       ACEData::getRawData
;       ACEData::getHmzData
;       ACEData::getBlocks
;       ACEData::areFlagsDefined
;       ACEData::setTimeFlag
;       ACEData::setSpaceFlag
;       ACEData::setErrorFlag
;       ACEData::areBlocksDefined
;       ACEData::blockTemplate
;       ACEData::resetBlocks
;       ACEData::emptyBlock
;       ACEData::blockPasses
;       ACEData::reduceCache
;       ACEData::reduceBlock
;       ACEData::calcBlocks
;       ACEData::makeContour
;       ACEData::makeMeridianPlot
;       ACEData::makeTimeSeries
;       ACEData::makeHeightTimePlot
;       ACEData::makeMeanSD
;       ACEData::makeSpotMap
;
;       Deprecated:
;       ACEData::isOpMaskDefined
;       ACEData::setOpMask       -- superseded by ACEData::blockPasses
;
; DATA LAYOUT:
;
;       Any number of molecules can be harmonized together. In the harmonized
;       structure returned by ACEData::getHmzData, and in the per-occultation
;       structure returned by ACEData::getOcc, the molecules stack along the
;       third dimension:
;
;           mixrat[nLevels,nPass,nMols]
;           mixrat_err[nLevels,nPass,nMols]
;           ratio[nLevels,nPass,nMols-1]    each molecule over the primary
;           names[nMols]
;
;       In the reduced data returned by ACEData::getBlocks, the statistics are
;       indexed the same way:
;
;           mean[nMols]      sd[nMols]      bw[5,nMols]
;           meanRatio[nMols-1]  sdRatio[...]  bwRatio[5,...]
;
;       Each block also carries the scalar names that predate this: mean0, sd0,
;       bw0, and for two or more molecules mean1, sd1, bw1, meanRat, sdRat and
;       bwRat. Those track the first molecule and the first ratio, and are kept
;       so that existing scripts and saved block arrays keep working.
;
; USAGE:
;   These programs are designed for use with ACE-FTS v4.1 retrievals. There is no guarantee they
;   will work properly with other versions of the data. Ideally, the user will only need to interact
;   with the following methods to generate plots and analyze data:
;   
;       ACEData::loadFiles
;       ACEData::harmonize
;       ACEData::setParam
;       ACEData::calcBlocks
;       ACEData::makeContour
;       
;   The following sample call will generate contour plots of CO2 at 16.5 km, then destroy the object.
;
;       myData=obj_new('ACEData')
;       myData.loadFiles,['ACEFTS_L2_v4p1_CO2.nc']
;       myData.harmonize,[0]
;       myData.setParam,'ctStandard'
;       myData.calcBlocks
;       myData.makeContour,ratPlot,mol0Plot,16.5
;       obj_destroy,myData
;       
; DEPENDENCIES:
; 
;       The IDL NCDF routines are used to load the ACE data from file.
;    
;       
; MODIFICATION HISTORY:
;       Written by: Ben Clouser, 03/3/21
;   
;       v0.1.0  03/04/21
;       v0.2.0
;         Added PT
;         Added ACEData::makeMeridianPlot
;         Added ACEData::loadDMP
;         Fixed Error in orbIDs pointers
;       v0.2.1
;         Added support for v4.1/4.2 data
;       v0.2.2  07/01/22-07/11/22
;         Added Julian days to harmonize
;         Added averaging for continuous time blocks
;         Fixed makeMeridianPlot to display both mol and ratio plots
;         Set up makeHeightTimePlot
;       v0.2.3
;         Fixed makeContour to use interpolate instead of griddata. Needs testing on non-linux platforms.
;  
; TODO:
;               
;         
;-

;+
; =============================================================
;
; METHODNAME:
;       ACEData::Init
;
; PURPOSE:
;       The ACEData::Init function method initializes the
;       ACE Data object.
;
;       NOTE: Init methods are special lifecycle methods, and as such
;       cannot be called outside the context of object creation.  This
;       means that in most cases, you cannot call the Init method
;       directly.  There is one exception to this rule: If you write
;       your own subclass of this class, you can call the Init method
;       from within the Init method of the subclass.
;
; CALLING SEQUENCE:
;       myData = OBJ_NEW('ACEData')
;
;       or
;
;       Result = myData.[ACEData::]Init()
;
; OPTIONAL INPUTS:  NONE
;
;
; KEYWORD PARAMETERS: NONE
;
; OUTPUTS:
;       1: successful, 0: unsuccessful.
;
; EXAMPLE:
;       myData = OBJ_NEW('ACEData')
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-


FUNCTION ACEData::Init

  COMPILE_OPT IDL2

  self.raw=PTR_NEW(/ALLOCATE)
  self.hmz=PTR_NEW(/ALLOCATE)
  self.dmp=PTR_NEW(/ALLOCATE)
  RETURN,1
  
END

;+
; =============================================================
;
; METHODNAME:
;       IDLgrLegend::Cleanup
;
; PURPOSE:
;       The ACEData::Cleanup procedure method preforms all cleanup
;       on the object.
;
;       NOTE: Cleanup methods are special lifecycle methods, and as such
;       cannot be called outside the context of object destruction.  This
;       means that in most cases, you cannot call the Cleanup method
;       directly.  There is one exception to this rule: If you write
;       your own subclass of this class, you can call the Cleanup method
;       from within the Cleanup method of the subclass.
;
; CALLING SEQUENCE:
;       OBJ_DESTROY, myData
;
;       or
;
;       myData.[ACEData::]Cleanup
;
; INPUTS:
;       There are no inputs for this method.
;
; KEYWORD PARAMETERS:
;       There are no keywords for this method.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO ACEData::Cleanup

  COMPILE_OPT IDL2

  IF PTR_VALID(self.raw) THEN PTR_FREE,self.raw
  IF PTR_VALID(self.hmz) THEN PTR_FREE,self.hmz
  IF PTR_VALID(self.flags) THEN PTR_FREE,self.flags
  IF PTR_VALID(self.opMask) THEN PTR_FREE,self.opMask
  IF PTR_VALID(self.blocks) THEN BEGIN
    ;----- N_ELEMENTS rather than subscripting [0], which errors when the heap
    ; variable was allocated but never assigned. -------;
    IF N_ELEMENTS(*self.blocks) GT 0 THEN PTR_FREE,(*self.blocks).orbIDs
    PTR_FREE,self.blocks
  ENDIF
  IF PTR_VALID(self.dmp) THEN PTR_FREE,self.dmp

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::freeRaw
;
; PURPOSE:
;       The ACEData::freeRaw method procedure
;       frees the pointer to the raw data loaded into the object
;       by the ACEData::loadFiles procedure. The raw data can take up significant
;       space in memory, especially if several instances of the ACEData object
;       are active
;
; CALLING SEQUENCE:
;       myData.[IDLgrLegend::]freeRaw
;
; INPUTS:
;       NONE
;
; EXAMPLE:
;       myData.freeRaw
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO ACEData::freeRaw

  COMPILE_OPT IDL2

  IF self.raw ne !null then begin
    PTR_FREE,self.raw
    self.raw=PTR_NEW()
    result=1
  ENDIF ELSE BEGIN
    result=0
  ENDELSE

  ;RETURN,result

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::loadFiles
;
; PURPOSE:
;       The ACEData::loadFiles procedure method loads the raw data from one
;       or more ACE-FTS v4.1 data file into the data object.
;
;       NOTE: The data is loaded essentially 1 for 1 from the data file into
;       the data object, so the data object can easily become quite large. Consider
;       using the freeRaw procedure to free the pointers to the raw data once it has
;       been processed to conserve memory.
;       
;       The data is loaded into a list of structures, with one element of the list
;       for each filename provided by the user.
;
; CALLING SEQUENCE:
;
;       myData.loadFiles,fileNames
;       
; INPUTS:
;   
;       fileNames is an array of strings containing the paths to the files to be loaded
;       into memory.
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: NONE
;
; EXAMPLE:
; 
;       myData.loadFiles,['ACEFTS_L2_v4p1_H2O.nc','ACEFTS_L2_v4p1_H2O_162.nc']
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Added support for v4p1p2, 12/8/21
;-

PRO ACEData::loadFiles,fnames,version=version,primeMerid=primeMerid

  COMPILE_OPT IDL2

  ;----- The VMR for the molecule sits at a fixed variable index in the ACE-FTS
  ; L2 files, with its uncertainty at the following index. That index moved
  ; between data versions, which is what the version keyword selects. -------;
  IF KEYWORD_SET(version) THEN BEGIN
    CASE version OF
      'p1p2': molnameid=11
      ELSE: molnameid=10
    ENDCASE
  ENDIF ELSE BEGIN
    molnameid=10
  ENDELSE

  n_files=N_ELEMENTS(fnames)
  mollist=LIST()

  FOR i=0,n_files-1 DO BEGIN

    IF ~FILE_TEST(fnames[i],/READ) THEN BEGIN
      PRINT,'ERROR: CANNOT READ FILE: '+fnames[i]
      CONTINUE
    ENDIF

    id=NCDF_OPEN(fnames[i])

    ;----- Release the file handle even if a variable turns out to be missing,
    ; then let the error propagate. -------;
    CATCH,ncerr
    IF ncerr NE 0 THEN BEGIN
      CATCH,/CANCEL
      NCDF_CLOSE,id
      MESSAGE,/REISSUE_LAST
    ENDIF

    mol=NCDF_VARINQ(id,molnameid)

    alt=self.ncGet(id,'altitude')
    srss=self.ncGet(id,'sunset_sunrise')
    year=self.ncGet(id,'year')
    month=self.ncGet(id,'month')
    day=self.ncGet(id,'day')
    hour=self.ncGet(id,'hour')
    lat=self.ncGet(id,'latitude')
    lon=self.ncGet(id,'longitude')
    orbit=self.ncGet(id,'orbit')
    temp=self.ncGet(id,'temperature')
    tempfit=self.ncGet(id,'temperature_fit')
    pres=self.ncGet(id,'pressure')

    NCDF_VARGET,id,molnameid,mixrat
    NCDF_VARGET,id,molnameid+1,mixrat_err

    CATCH,/CANCEL
    NCDF_CLOSE,id

    IF KEYWORD_SET(primeMerid) THEN BEGIN
      jneg=WHERE(lon LT 0,nneg)
      IF nneg GT 0 THEN lon[jneg]+=360
    ENDIF

    mollist.add,{name:mol.name,npass:N_ELEMENTS(lon),alt:alt,srss:srss, $
                 orbit:orbit,year:year,month:month,day:day,hour:hour, $
                 lat:lat,lon:lon,mixrat:mixrat,mixrat_err:mixrat_err, $
                 temp:temp,tempfit:tempfit,pres:pres}

  ENDFOR

  IF mollist.count() EQ 0 THEN BEGIN
    PRINT,'ERROR: NO FILES LOADED'
    RETURN
  ENDIF

  (*self.raw)=mollist

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::ncGet
;
; PURPOSE:
;       The ACEData::ncGet function method is a private method and is not
;       intended to be called directly. It reads a NetCDF variable by name and
;       reports the name in the error message if the variable is absent, which
;       NCDF_VARGET on an id of -1 does not do.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 8/4/25
;-

FUNCTION ACEData::ncGet,id,vname

  COMPILE_OPT IDL2

  vid=NCDF_VARID(id,vname)
  IF vid EQ -1 THEN MESSAGE,'Variable not found in file: '+vname
  NCDF_VARGET,id,vid,data
  RETURN,data

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::loadDMP
;
; PURPOSE:
;       The ACEData::loadFiles procedure method loads the JPL Derived Meteorological Product file
;
; CALLING SEQUENCE:
;
;       myData.loadDMP
;
; INPUTS: NONE  (DMP path is a parameter which can be edited using ACEData::setParam)
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: NONE
;
; EXAMPLE:
;
;       myData.loadDMP
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 4/15/21
;-

PRO ACEData::loadDMP

  COMPILE_OPT IDL2

  IF ~FILE_TEST(self.params.dmp_path,/READ) THEN BEGIN
    PRINT,'ERROR: CANNOT READ DMP FILE: '+self.params.dmp_path
    PRINT,'       Set it with setParam,''dmp_path'',<path or version alias>.'
    RETURN
  ENDIF

  id=NCDF_OPEN(self.params.dmp_path)

  CATCH,ncerr
  IF ncerr NE 0 THEN BEGIN
    CATCH,/CANCEL
    NCDF_CLOSE,id
    MESSAGE,/REISSUE_LAST
  ENDIF

  srss=self.ncGet(id,'Sunset Sunrise')
  orbit=self.ncGet(id,'Orbit')
  alt=self.ncGet(id,'Altitude')

  pres=self.ncGet(id,'Pressure')
  uwind=self.ncGet(id,'U')
  vwind=self.ncGet(id,'V')
  temp=self.ncGet(id,'Temperature')
  theta=self.ncGet(id,'Theta')
  gph=self.ncGet(id,'GPH')
  pv=self.ncGet(id,'PV')
  spv=self.ncGet(id,'sPV')
  rv=self.ncGet(id,'RV')
  dtdz=self.ncGet(id,'dT_dZ')
  ss=self.ncGet(id,'Static_Stability')
  eql=self.ncGet(id,'EqL')
  hpvg=self.ncGet(id,'HorsPVgrad')
  apvg=self.ncGet(id,'avgsPVgrad')
  tg=self.ncGet(id,'Tgrad')
  msf=self.ncGet(id,'MSF')
  o3=self.ncGet(id,'AssimO3')

  CATCH,/CANCEL
  NCDF_CLOSE,id

  ;----- The level count comes from the file rather than being assumed to be
  ; 75, so a DMP on a different vertical grid does not silently reshape into
  ; nonsense. -------;
  n_orbs=N_ELEMENTS(orbit)
  n_alts=N_ELEMENTS(alt)

  IF n_alts*n_orbs NE N_ELEMENTS(pres) THEN BEGIN
    PRINT,'ERROR: DMP DIMENSIONS INCONSISTENT. '+ $
          STRTRIM(n_alts,2)+' altitudes x '+STRTRIM(n_orbs,2)+' orbits '+ $
          'does not match '+STRTRIM(N_ELEMENTS(pres),2)+' pressure values.'
    RETURN
  ENDIF

  str={srss:srss,orbit:orbit,alt:alt, $
       pres:REFORM(pres,n_alts,n_orbs),U:REFORM(uwind,n_alts,n_orbs), $
       V:REFORM(vwind,n_alts,n_orbs),temp:REFORM(temp,n_alts,n_orbs), $
       theta:REFORM(theta,n_alts,n_orbs),GPH:REFORM(gph,n_alts,n_orbs), $
       pv:REFORM(pv,n_alts,n_orbs),spv:REFORM(spv,n_alts,n_orbs), $
       rv:REFORM(rv,n_alts,n_orbs),dtdz:REFORM(dtdz,n_alts,n_orbs), $
       ss:REFORM(ss,n_alts,n_orbs),eql:REFORM(eql,n_alts,n_orbs), $
       hpvg:REFORM(hpvg,n_alts,n_orbs),apvg:REFORM(apvg,n_alts,n_orbs), $
       tg:REFORM(tg,n_alts,n_orbs),msf:REFORM(msf,n_alts,n_orbs), $
       o3:REFORM(o3,n_alts,n_orbs)}

  (*self.DMP)=str

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::harmonize
;
; PURPOSE:
;       The ACEData::harmonize procedure method operates on the raw data loaded
;       into the ACE Data object by the ACEData::loadFiles procedure method. This method
;       finds retrievals common to ALL of the requested molecules and stores them in a
;       new data structure. Any number of molecules is supported. The first molecule in
;       the elements array is the 'primary' molecule, and ACEData::harmonize calculates
;       the ratio of each subsequent molecule to it, so N molecules give N-1 ratios.
;
;       Occultations are matched on orbit number together with the sunrise/sunset
;       flag. An occultation is kept only if every requested molecule has a retrieval
;       for it, so adding a molecule can only narrow the result. The harmonized arrays
;       come out in the primary molecule's occultation order.
;
;       See the DATA LAYOUT section of the class header for how the molecules and
;       ratios are indexed.
;
;
; CALLING SEQUENCE:
;
;       myData.harmonize,elements
;
; INPUTS:
;
;       Elements is an array of integers that determiines which of the molecules provided to
;       ACEData::loadFiles is supposed to be the primary molecule. The integers in the array
;       reference the array of filename strings provided to ACEData::loadFiles. For example,
;       assume the following calls are made:
;       
;       myData=OBJ_NEW('ACEData')
;       myData.loadFiles,['ACEFTS_L2_v4p1_H2O.nc','ACEFTS_L2_v4p1_H2O_162.nc'] 
;       
;       If the user subsequently calls
;       
;       myData.harmonize,[0,1]
;       
;       then H2O will be the primary moecule. However if the user calls
;       
;       myData.harmonize,[1,0]
;       
;       then H2O_162 will be the primary molecule.
;
;       The array can name any number of molecules. Given four files loaded in
;       order, 
;       
;       myData.harmonize,[0,1,2,3]
;       
;       keeps the occultations common to all four, with molecule 0 as primary and
;       three ratios against it. The order of the array sets the order of the
;       mixrat, names, mean and bw dimensions, so harmonize,[2,0,1] makes the
;       third file loaded the primary.
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: 
;
;       addDMP incorporates selected DMP fields into the harmonized data structure.
;
; EXAMPLE:
;
;       myData.harmonize,[0,1]
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   BWC, 6/8/24: 
;-

PRO ACEData::harmonize,elements,addDMP=addDMP

  COMPILE_OPT IDL2

  IF ~PTR_VALID(self.raw) || N_ELEMENTS(*self.raw) EQ 0 THEN BEGIN
    PRINT,'ERROR: NO RAW DATA. Call loadFiles first.'
    RETURN
  ENDIF

  nRaw=(*self.raw).count()
  IF MIN(elements) LT 0 OR MAX(elements) GE nRaw THEN BEGIN
    PRINT,'ERROR: ELEMENTS OUT OF RANGE. '+STRTRIM(nRaw,2)+' molecule(s) loaded.'
    RETURN
  ENDIF

  n_mols=N_ELEMENTS(elements)
  nrat=n_mols-1

  mol0=(*self.raw)[elements[0]]
  nlev=N_ELEMENTS(mol0.alt)

  ;----- An occultation is identified by its orbit number together with the
  ; sunrise/sunset flag. Matching is done on that key by sorting each molecule
  ; once and using VALUE_LOCATE. The previous implementation scanned one
  ; molecule in full for every occultation of the other and grew its index
  ; arrays by reallocation, which is quadratic, and it only ever looked at
  ; elements[0] and elements[1] however many were passed. -------;
  keys=LIST()
  srts=LIST()
  FOR m=0,n_mols-1 DO BEGIN
    mol=(*self.raw)[elements[m]]
    kk=LONG64(mol.orbit)*2LL+LONG64(mol.srss)
    keys.add,kk
    srts.add,SORT(kk)
  ENDFOR

  ;----- Start from the primary molecule and intersect with each of the others
  ; in turn. Working down from the primary keeps the harmonized arrays in the
  ; primary molecule's occultation order. -------;
  ckeys=keys[0]
  keep0=LINDGEN(N_ELEMENTS(ckeys))

  FOR m=1,n_mols-1 DO BEGIN
    ks=(keys[m])[srts[m]]
    loc=VALUE_LOCATE(ks,ckeys)
    g=WHERE(loc GE 0,ng)
    IF ng GT 0 THEN hit=WHERE(ks[loc[g]] EQ ckeys[g],nhit) ELSE nhit=0
    IF nhit EQ 0 THEN BEGIN
      PRINT,'ERROR: NO OCCULTATIONS COMMON TO ALL '+STRTRIM(n_mols,2)+' MOLECULES'
      RETURN
    ENDIF
    sel=g[hit]
    ckeys=ckeys[sel]
    keep0=keep0[sel]
  ENDFOR

  npass=N_ELEMENTS(ckeys)

  ;----- Index into each molecule's own arrays for the surviving occultations.
  ; Every key is now known to be present in every molecule, so VALUE_LOCATE
  ; lands on an exact match. -------;
  ind=LONARR(npass,n_mols)
  ind[*,0]=keep0
  FOR m=1,n_mols-1 DO BEGIN
    ks=(keys[m])[srts[m]]
    ind[*,m]=(srts[m])[VALUE_LOCATE(ks,ckeys)]
  ENDFOR

  ;----- Mixing ratios for every molecule stack along the third dimension, and
  ; each subsequent molecule is ratioed against the primary. -------;
  names=STRARR(n_mols)
  mixrat=FLTARR(nlev,npass,n_mols)
  mixrat_err=FLTARR(nlev,npass,n_mols)

  FOR m=0,n_mols-1 DO BEGIN
    mol=(*self.raw)[elements[m]]
    names[m]=mol.name
    mixrat[*,*,m]=mol.mixrat[*,ind[*,m]]
    mixrat_err[*,*,m]=mol.mixrat_err[*,ind[*,m]]
  ENDFOR

  IF nrat GT 0 THEN BEGIN
    ratio=FLTARR(nlev,npass,nrat)
    FOR m=1,n_mols-1 DO ratio[*,*,m-1]=mixrat[*,*,m]/mixrat[*,*,0]
  ENDIF

  juldays=JULDAY(mol0.month[keep0],mol0.day[keep0],mol0.year[keep0],mol0.hour[keep0],0d0,0d0)

  ;----- Potential temperature and temperature come from the DMP when one is
  ; loaded, and from the retrieval itself otherwise. -------;
  isDMP=0B
  IF KEYWORD_SET(addDMP) THEN BEGIN
    IF ~PTR_VALID(self.dmp) THEN BEGIN
      PRINT,'ERROR: POINTER TO DMP NOT VALID'
    ENDIF ELSE IF N_ELEMENTS(*self.dmp) EQ 0 THEN BEGIN
      PRINT,'ERROR: DMP NOT LOADED'
    ENDIF ELSE BEGIN
      isDMP=1B
    ENDELSE
  ENDIF

  IF isDMP THEN BEGIN

    self.params.isDMP=1B

    nDMPlev=(SIZE((*self.dmp).theta,/DIMENSIONS))[0]
    IF nDMPlev GT nlev THEN BEGIN
      PRINT,'ERROR: DMP GRID HAS MORE LEVELS ('+STRTRIM(nDMPlev,2)+ $
            ') THAN THE RETRIEVAL GRID ('+STRTRIM(nlev,2)+')'
      RETURN
    ENDIF

    ;----- A matching DMP entry is not guaranteed to exist. Occultations without
    ; one stay NaN rather than silently picking up the last DMP column, and so
    ; do levels above the DMP's own vertical grid. -------;
    pt=FLTARR(nlev,npass)+!VALUES.F_NAN
    t_dmp=FLTARR(nlev,npass)+!VALUES.F_NAN

    dk=LONG64((*self.dmp).orbit)*2LL+LONG64((*self.dmp).srss)
    ds=SORT(dk)
    dks=dk[ds]

    matched=BYTARR(npass)
    dind=LONARR(npass)

    loc=VALUE_LOCATE(dks,ckeys)
    g=WHERE(loc GE 0,ng)
    IF ng GT 0 THEN BEGIN
      hit=WHERE(dks[loc[g]] EQ ckeys[g],nhit)
      IF nhit GT 0 THEN BEGIN
        matched[g[hit]]=1B
        dind[g[hit]]=ds[loc[g[hit]]]
      ENDIF
    ENDIF

    jDMP=WHERE(matched,nDMP,NCOMPLEMENT=nNoDMP)
    IF nDMP GT 0 THEN BEGIN
      pt[0:nDMPlev-1,jDMP]=(*self.dmp).theta[*,dind[jDMP]]
      t_dmp[0:nDMPlev-1,jDMP]=(*self.dmp).temp[*,dind[jDMP]]
    ENDIF
    IF nNoDMP GT 0 THEN PRINT,'WARNING: NO DMP MATCH FOR '+STRTRIM(nNoDMP,2)+ $
                              ' OF '+STRTRIM(npass,2)+' OCCULTATIONS'

  ENDIF ELSE BEGIN

    pt=mol0.temp[*,keep0]*(1d0/mol0.pres[*,keep0])^.286d0
    t_dmp=mol0.pres[*,keep0]*!VALUES.F_NAN

  ENDELSE

  molout={nmols:n_mols,names:names,npass:npass,alt:mol0.alt, $
          srss:mol0.srss[keep0],orbit:mol0.orbit[keep0], $
          year:mol0.year[keep0],month:mol0.month[keep0], $
          day:mol0.day[keep0],hour:mol0.hour[keep0],juldays:juldays, $
          lat:mol0.lat[keep0],lon:mol0.lon[keep0], $
          mixrat:mixrat,mixrat_err:mixrat_err, $
          temp:mol0.temp[*,keep0],pres:mol0.pres[*,keep0], $
          pt:pt,t_dmp:t_dmp}

  IF nrat GT 0 THEN molout=CREATE_STRUCT(molout,'ratio',ratio)

  (*self.hmz)=molout

end

;+
; =============================================================
;
; METHODNAME:
;       ACEData::setParam
;
; PURPOSE:
;       The ACEData::setParam procedure method sets the parameters
;       used to determine the global flags, as well as the number of spatial
;       and temporal bins, and the error criteria.
;
;
; CALLING SEQUENCE:
;
;       myData.setParam,pName,pVal
;
; INPUTS:
;
;       pName is a string containing the name of the parameter to be modified.
;       
;       pVal is the value that the parameter should be set to.
;       
;       CURRENTLY SUPPORTED PARAMETERS:
;       
;       latLims
;       lonLims
;       altLims
;       nLats
;       nLons
;       nAlts
;       yearLims
;       monthLims
;       dayLims
;       hourLims
;       nTimes
;       eLims             --            the error limit beyond which data is cut
;       cutType           --            the manner in which data is cut, currently only 'pointwise' is supported
;       
;       ctStandard        --            Instead of setting an individual parameter, this input will load a suite of parameters considered
;                                       typical for making contour plots of JJA data. In this usage, pVal does not need to be set.
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: NONE
;
; EXAMPLE:
;
;       myData.setParam,'latLims',[-70,70]
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO ACEData::setParam,pName,pVal

  COMPILE_OPT IDL2

  ;----- Exact match on the tag name. The old test was a substring match, so
  ; setParam,'nLat',... silently wrote to NLATS, and any short name could land
  ; on the wrong field. -------;
  names=TAG_NAMES(self.params)
  upName=STRUPCASE(pName)
  isParam=WHERE(names EQ upName,nParam)

  if nParam gt 0 then begin
    IF upName EQ 'DMP_PATH' THEN pVal=self.resolveDMPPath(pVal)

    self.params.(isParam[0])=pVal

    ;----- Mark the dependent stages dirty so calcBlocks knows to recompute. --;
    IF TOTAL(upName EQ ['LATLIMS','LONLIMS','ALTLIMS','NLATS','NLONS','NALTS']) GT 0 THEN self.params.sMod=1
    IF TOTAL(upName EQ ['YEARLIMS','MONTHLIMS','DAYLIMS','HOURLIMS','JULLIMS','NTIMES']) GT 0 THEN self.params.tMod=1
    IF TOTAL(upName EQ ['ELIM','CUTTYPE']) GT 0 THEN self.params.eMod=1
    self.params.bMod=1
  endif else begin
    CASE pName OF
      'ctStandard': BEGIN
        self.params.latlims=[-70,70]
        self.params.lonlims=[-180,180]
        self.params.altlims=[16,17]
        self.params.nLats=10
        self.params.nLons=15
        self.params.nAlts=1
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[6,8]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
        'htStandard': BEGIN
          self.params.latlims=[-60,60]
          self.params.lonlims=[-180,180]
          self.params.altlims=[12,40]
          self.params.nLats=1
          self.params.nLons=1
          self.params.nAlts=28
          self.params.yearlims=[2004,2022]
          self.params.monthlims=[1,12]
          self.params.daylims=[1,31]
          self.params.hourlims=[1d-4,24d0-1d-4]
          self.params.nTimes=18
          self.params.eLim=1d-3
          self.params.cuttype='pointwise'
          self.params.sMod=1
          self.params.tMod=1
          self.params.eMod=1
          self.params.bMod=1
          ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
        END
      'monsProf': BEGIN
        self.params.latlims=[5,35]
        self.params.lonlims=[60,120]
        self.params.altlims=[8,30]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=22
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[6,8]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'OsanProf': BEGIN
        self.params.latlims=[10,45]
        self.params.lonlims=[100,160]
        self.params.altlims=[8,30]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=22
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[6,8]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'NAMProf': BEGIN
        self.params.latlims=[10,50]
        self.params.lonlims=[-130,-70]
        self.params.altlims=[8,30]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=22
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[6,8]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'TRANSProf': BEGIN
        self.params.latlims=[40,70]
        self.params.lonlims=[-230,-120]
        self.params.altlims=[8,30]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=22
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[7,9]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'CRProf': BEGIN
        self.params.latlims=[-10,20]
        self.params.lonlims=[-100,-60]
        self.params.altlims=[8,30]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=22
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[7,9]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END


      'NHTime': BEGIN
        self.params.latlims=[0,90]
        self.params.lonlims=[-180,180]
        self.params.altlims=[20,25]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=1
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[1,12]
        self.params.daylims=[0d0,0d0]
        self.params.hourlims=[0d0,0d0]
        self.params.nTimes=18.*12.
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'SHTime': BEGIN
        self.params.latlims=[-90,0]
        self.params.lonlims=[-180,180]
        self.params.altlims=[30,31]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=1
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[1,12]
        self.params.daylims=[0d0,0d0]
        self.params.hourlims=[0d0,0d0]
        self.params.nTimes=18.*24.
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'htHunga': BEGIN
        self.params.latlims=[-45,0]
        self.params.lonlims=[-180,180]
        self.params.altlims=[12,40]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=28
        self.params.yearlims=[2022,2023]
        self.params.monthlims=[1,12]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
        self.params.solidBlock=1
      END
      'SHPole': BEGIN
        self.params.latlims=[-85,-60]
        self.params.lonlims=[-180,180]
        self.params.altlims=[12,40]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=28
        self.params.yearlims=[2004,2022]
        self.params.monthlims=[7,9]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'SHPole2020': BEGIN
        self.params.latlims=[-85,-60]
        self.params.lonlims=[-180,180]
        self.params.altlims=[12,40]
        self.params.nLats=1
        self.params.nLons=1
        self.params.nAlts=28
        self.params.yearlims=[2020,2020]
        self.params.monthlims=[7,9]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
      'allMerid': BEGIN
        self.params.latlims=[-85,85]
        self.params.lonlims=[-180,180]
        self.params.altlims=[5,40]
        self.params.nLats=17
        self.params.nLons=1
        self.params.nAlts=35
        self.params.yearlims=[2004,2024]
        self.params.monthlims=[1,12]
        self.params.daylims=[1,31]
        self.params.hourlims=[1d-4,24d0-1d-4]
        self.params.nTimes=1
        self.params.eLim=1d-3
        self.params.cuttype='pointwise'
        self.params.sMod=1
        self.params.tMod=1
        self.params.eMod=1
        self.params.bMod=1
        self.params.solidBlock=1
        ;self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
      END
     ELSE: print,'Parameter name not found'
     ENDCASE
  endelse

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::getParam
;
; PURPOSE:
;       The ACEData::getParam function method returns the value of the
;       selected parameter.
;
;
; CALLING SEQUENCE:
;
;       result = myData.getParam(pName)
;
; INPUTS:
;
;       pName is a string containing the name of the parameter to be modified.
;
;       CURRENTLY SUPPORTED PARAMETERS:
;
;       latLims
;       lonLims
;       altLims
;       nLats
;       nLons
;       nAlts
;       yearLims
;       monthLims
;       dayLims
;       hourLims
;       nTimes
;       eLims             --            the error limit beyond which data is cut
;       cutType           --            the manner in which data is cut, currently only 'pointwise' is supported
;
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: NONE
; 
; OUTPUTS:
; 
;         The output is the value of the parameter selected with pName. The type of data is 
;         that associated with pName.
;
; EXAMPLE:
;
;       latLims=myData.getParam('latLims')
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Added print all, 7/12/21
;-

FUNCTION ACEData::getParam,pName

  COMPILE_OPT IDL2

IF N_PARAMS() LT 1 THEN BEGIN
  
  PRINT,'Latitude Range',self.params.latLims
  PRINT,'Latitude Boxes',self.params.nLats
  PRINT,'Longitude Range',self.params.lonLims
  PRINT,'Longitude Boxes',self.params.nLons
  PRINT,'Altitude Range',self.params.altLims
  PRINT,'Altitude Boxes',self.params.nAlts
  PRINT,'Year Range',self.params.yearLims
  PRINT,'Month Range',self.params.monthLims
  PRINT,'Day Range',self.params.dayLims
  PRINT,'Hour Range',self.params.hourLims
  PRINT,'Julian Range',self.params.julLims
  PRINT,'Time Boxes',self.params.nTimes
  
  RETURN,1
  
ENDIF ELSE BEGIN

  names=TAG_NAMES(self.params)
  isParam=WHERE(names EQ STRUPCASE(pName),nParam)
  if nParam gt 0 then begin
    result=self.params.(isParam[0])
  endif else begin
    result=!values.f_nan
    print,'Parameter name not found'
  endelse

  RETURN,result

ENDELSE

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::resolveDMPPath
;
; PURPOSE:
;       The ACEData::resolveDMPPath function method is a private method and is
;       not intended to be called directly. It expands the short DMP version
;       aliases into full paths.
;
;       The aliases resolve against the ACE_DATA_ROOT environment variable:
;
;           setenv,'ACE_DATA_ROOT=/media/benjamin/Data/ACE'
;           myData.setParam,'dmp_path','v5'
;
;       If ACE_DATA_ROOT is not set, the original hardcoded root is used so
;       that existing scripts keep working unchanged. Anything that is not a
;       recognised alias is passed straight through as a literal path.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 8/4/25
;-

FUNCTION ACEData::resolveDMPPath,pVal

  COMPILE_OPT IDL2

  root=GETENV('ACE_DATA_ROOT')
  IF root EQ '' THEN root='/media/benjamin/Data/ACE'

  CASE STRLOWCASE(pVal) OF
    'v4p1':   subdir='ACE-DMPSv4p120210415111208'
    'v4p1p2': subdir='ACE-DMPSv4p1p2'
    'v5':     subdir='ACE-DMPSv5p2'
    ELSE:     RETURN,pVal
  ENDCASE

  IF STRLOWCASE(pVal) EQ 'v5' THEN BEGIN
    fname='ACEFTS_L2_v5p2_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
  ENDIF ELSE BEGIN
    fname='ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
  ENDELSE

  RETURN,FILEPATH(fname,ROOT_DIR=root,SUBDIRECTORY=subdir)

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::getOcc
;
; PURPOSE:
;       The ACEData::getOcc function method returns the retrievals for
;       one or more orbit IDs supplied by the user. Orbits not found in the
;       harmonized data are indicated to the user by text.
;
;
; CALLING SEQUENCE:
;
;       result = myData.getOcc(orbIDs)
;
; INPUTS:
;
;       orbIDs is an array of integers representing orbit IDs.
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: NONE
;
; OUTPUTS:
;
;         The output is a single STRUCTURE containing ARRAYS of data for the orbits found. Note that
;         the number of elements in the arrays may be smaller than the number of elements in orbIDs
;         since orbits are not guaranteed to exist in the harmonized data. (Additionaly, there are several
;         tags in the structure such as altitude that do not vary with the number of orbits.) 
;
; EXAMPLE:
;
;       orbits = myData.getOcc(orbIDs)
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

FUNCTION ACEData::getOcc,orbids,SEARCH=search

  COMPILE_OPT IDL2

  nOccs=n_elements(orbids)

IF KEYWORD_SET(search) THEN BEGIN
  inds=[]
  missing=[]

  ;----- Accumulate only the orbits that were actually found. The missing IDs are
  ; tracked separately because one orbit ID can match more than one index, so the
  ; positions in inds do not correspond to positions in orbids. -------;
  for i=0,nOccs-1 do begin
    x=where((*self.hmz).orbit eq orbids[i],nx)
    if nx gt 0 then inds=[inds,x] else missing=[missing,orbids[i]]
  endfor

  if n_elements(missing) gt 0 then print,'Orbits ',missing,' not found.'

  if n_elements(inds) eq 0 then begin
    print,'No requested orbits found in the harmonized data.'
    RETURN,!NULL
  endif

  nOccs=n_elements(inds)
ENDIF ELSE BEGIN
  inds=orbids
ENDELSE

;----- Mirrors the harmonized layout: mixing ratios stack along the third
; dimension and ratios are held separately, so one build covers any number of
; molecules. This replaced a CASE on nmols that had branches for 1 and 2 only
; and no ELSE. -------;
out={nmols:(*self.hmz).nmols,names:(*self.hmz).names,npass:nOccs, $
     alt:(*self.hmz).alt,srss:(*self.hmz).srss[inds], $
     orbit:(*self.hmz).orbit[inds],year:(*self.hmz).year[inds], $
     month:(*self.hmz).month[inds],day:(*self.hmz).day[inds], $
     hour:(*self.hmz).hour[inds],lat:(*self.hmz).lat[inds], $
     lon:(*self.hmz).lon[inds], $
     mixrat:(*self.hmz).mixrat[*,inds,*], $
     mixrat_err:(*self.hmz).mixrat_err[*,inds,*], $
     temp:(*self.hmz).temp[*,inds],pres:(*self.hmz).pres[*,inds], $
     pt:(*self.hmz).pt[*,inds],t_dmp:(*self.hmz).t_dmp[*,inds]}

IF (*self.hmz).nmols GT 1 THEN out=CREATE_STRUCT(out,'ratio',(*self.hmz).ratio[*,inds,*])

RETURN,out

END

FUNCTION ACEData::getRawData

  COMPILE_OPT IDL2

  IF PTR_VALID(self.raw) && N_ELEMENTS(*self.raw) GT 0 THEN BEGIN
    RETURN,(*self.raw)
  ENDIF ELSE BEGIN
    PRINT,'Raw data not defined.'
    RETURN,!NULL
  ENDELSE

END

FUNCTION ACEData::getHmzData

  COMPILE_OPT IDL2

  IF PTR_VALID(self.hmz) && N_ELEMENTS(*self.hmz) GT 0 THEN BEGIN
    RETURN,(*self.hmz)
  ENDIF ELSE BEGIN
    PRINT,'Harmonized data not defined.'
    RETURN,!NULL
  ENDELSE


END

FUNCTION ACEData::getBlocks

  COMPILE_OPT IDL2

  IF PTR_VALID(self.blocks) && N_ELEMENTS(*self.blocks) GT 0 THEN BEGIN
    RETURN,(*self.blocks)
  ENDIF ELSE BEGIN
    PRINT,'Blocks not defined.'
    RETURN,!NULL
  ENDELSE

END


;+
; =============================================================
;
; METHODNAME:
;       ACEData::areFlagsDefined
;
; PURPOSE:
;       The ACEData::areFlagsDefined procedure method is a private
;       method and is not intended to be called directly. It declares the global flag arrays
;       based on the size of the given data set.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-


PRO ACEData::areFlagsDefined

  COMPILE_OPT IDL2

  npass=(*self.hmz).npass
  nlev=N_ELEMENTS((*self.hmz).alt)

  if self.flags eq !null then begin

    self.flags=PTR_NEW(/ALLOCATE)
    ;----- Anonymous, not named. A named structure definition is cached for the
    ; whole IDL session, so {flags,...} sized by npass raised a conflicting
    ; definition error as soon as a second object with a different pass count
    ; was created. -------;
    B={tflag:boolarr(npass),sflag:boolarr(npass),eflag:boolarr(npass),aflag:boolarr(nlev,npass)}
    (*self.flags)=B
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::setTimeFlag
;
; PURPOSE:
;       The ACEData::setTimeFlag procedure method is a private
;       method and is not intended to be called directly. It checks the user defined limits on
;       the time range and bins, and calculates the appropriate global temporal flag.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Added support for julian days, 7/12/22
;-

pro ACEData::setTimeFlag,sunrise=sunrise,sunset=sunset

  COMPILE_OPT IDL2
 
  mth_lab=['J','F','M','A','M','J','J','A','S','O','N','D']
  mth_str=''

  self.areFlagsDefined
  
  IF self.params.solidBlock EQ 0 THEN BEGIN
    
    yearInds=indgen(self.params.yearlims[1]-self.params.yearlims[0]+1)+self.params.yearlims[0]
    yearIncl=intarr((*self.hmz).npass)
    for ii=0,n_elements(yearInds)-1 do begin
      jyr=where((*self.hmz).year eq yearInds[ii],njyr)
      if njyr gt 0 then yearIncl[jyr]=1
    endfor
    
    if self.params.monthlims[1] ne 0 then begin 
      monthIncl=INTARR((*self.hmz).npass)
      if self.params.monthlims[0] ge self.params.monthlims[1] THEN BEGIN 
      ;monthInds=[indgen(12.-self.params.monthlims[0]+1)+self.params.monthlims[0],indgen(self.params.monthlims[1])+1] 
      ;ENDIF ELSE BEGIN
      ;  monthInds=indgen(self.params.monthlims[1]-self.params.monthlims[0]+1)+self.params.monthlims[0]
      ;ENDELSE
      ;monthIncl=intarr((*self.hmz).npass)
      ;for ii=0,n_elements(monthInds)-1 do begin
      ;  jmth=where((*self.hmz).month eq monthInds[ii])
      ;  monthIncl[jmth]=1
      ;  mth_str+=mth_lab[monthInds[ii]-1]
      ;endfor
      jmth=WHERE((*self.hmz).month GE self.params.monthlims[0] OR (*self.hmz).month LE self.params.monthlims[1],njmth)
      mth_str=mth_lab[[indgen(12.-self.params.monthlims[0]+1)+self.params.monthlims[0],indgen(self.params.monthlims[1])+1]-1]
      ENDIF ELSE BEGIN
      jmth=WHERE((*self.hmz).month GE self.params.monthlims[0] AND (*self.hmz).month LE self.params.monthlims[1],njmth)
      mth_str=mth_lab[indgen(self.params.monthlims[1]-self.params.monthlims[0]+1)+self.params.monthlims[0]-1]
      ENDELSE
      if njmth gt 0 then monthIncl[jmth]=1
    endif else begin
      monthIncl=intarr((*self.hmz).npass)+1
    endelse
    
    ;stop
    
    if self.params.daylims[1] ne 0 then begin 
      if self.params.daylims[0] ge self.params.daylims[1] then begin
        ;add code here that can handle rolling over months... maybe this will have to be lumped together with monthlims?
      endif else begin
        dayInds=indgen(self.params.daylims[1]-self.params.daylims[0]+1)+self.params.daylims[0]
      endelse
      dayIncl=intarr((*self.hmz).npass)
      for ii=0,n_elements(dayInds)-1 do begin
        jdy=where((*self.hmz).day eq dayInds[ii],njdy)
        if njdy gt 0 then dayIncl[jdy]=1
        ;stop
      endfor
    endif else begin
      dayIncl=intarr((*self.hmz).npass)+1
    endelse
    
    if self.params.hourlims[1] ne 0 then begin
      if self.params.hourlims[0] le self.params.hourlims[1] then begin
        hourIncl=(((*self.hmz).hour ge self.params.hourlims[0]) AND ((*self.hmz).hour le self.params.hourlims[1]))
      endif else begin
        hourIncl=(((*self.hmz).hour le self.params.hourlims[1]) OR ((*self.hmz).hour ge self.params.hourlims[0]))
      endelse
    endif else begin
      hourIncl=intarr((*self.hmz).npass)+1
    endelse
  
    if keyword_set(sunrise) OR keyword_set(sunset) then begin
      if keyword_set(sunrise) then begin
        isSunrise=((*self.hmz).srss eq 0)
        (*self.flags).tflag=(yearIncl AND monthIncl AND dayIncl AND hourIncl AND isSunrise)
      endif
    
      if keyword_set(sunset) then begin
        isSunset=((*self.hmz).srss eq 1)
        (*self.flags).tflag=(yearIncl AND monthIncl AND dayIncl AND hourIncl AND isSunset)
      endif
    endif else begin
      (*self.flags).tflag=(yearIncl AND monthIncl AND dayIncl AND hourIncl)
    endelse

    mstr=''
    FOREACH mth,mth_str DO BEGIN
      mstr+=mth
    ENDFOREACH
    self.params.mth_str=mstr
  
  ENDIF ELSE BEGIN
    
    dayLast=julday(self.params.monthlims[1],self.params.daylims[1],self.params.yearlims[1],self.params.hourlims[1],0,0)
    dayFirst=julday(self.params.monthlims[0],self.params.daylims[0],self.params.yearlims[0],self.params.hourlims[0],0,0)

    ;----- WHERE returns indices, not a mask, so build the flag by subscripting
    ; a zeroed array rather than assigning the index list to tflag. -------;
    julIncl=WHERE((*self.hmz).juldays GE self.params.jullims[0] AND (*self.hmz).juldays LE self.params.jullims[1],njulIncl)
    (*self.flags).tflag[*]=0
    IF njulIncl GT 0 THEN (*self.flags).tflag[julIncl]=1
    
  ENDELSE

  self.params.tMod=0

end

;+
; =============================================================
;
; METHODNAME:
;       ACEData::setSpaceFlag
;
; PURPOSE:
;       The ACEData::setSpaceFlag procedure method is a private
;       method and is not intended to be called directly. It checks the user defined limits on
;       lat, lon, and alt ranges and bins, and calculates the appropriate global spatial flag.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO ACEData::setSpaceFlag

  COMPILE_OPT IDL2

  self.areFlagsDefined

  kran=WHERE((*self.hmz).lat GE self.params.latlims[0] AND (*self.hmz).lat LE self.params.latlims[1] AND (*self.hmz).lon ge self.params.lonlims[0] AND (*self.hmz).lon LE self.params.lonlims[1],nkran,complement=jran,ncomplement=njran)
  IF nkran GT 0 THEN (*self.flags).sflag[kran]=1
  IF njran GT 0 THEN (*self.flags).sflag[jran]=0
  
  self.params.sMod=0
  
END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::setErrorFlag
;
; PURPOSE:
;       The ACEData::setErrorFlag procedure method is a private
;       method and is not intended to be called directly. It checks the user defined error settings,
;        and calculates the appropriate global error flags.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Added cutByTemp and scanwise options, 5/29/21
;-

PRO ACEData::setErrorFlag,cutByTemp=cutByTemp,cutByRatio=cutByRatio

  COMPILE_OPT IDL2
  
  nlev=N_ELEMENTS((*self.hmz).alt)

  (*self.flags).eflag=1
  IF KEYWORD_SET(cutByTemp) THEN BEGIN
    k=WHERE((*self.hmz).temp LE cutByTemp[0],nk)
    IF nk GT 0 THEN BEGIN
      z=k MOD nlev
      v=WHERE(z GE cutByTemp[1] AND z LE cutByTemp[2],nv)
      IF nv GT 0 THEN BEGIN
        kOrb=k[v]/nlev
        (*self.flags).eflag[kOrb[uniq(kOrb)]]=0
      ENDIF
    ENDIF
    ;stop
  ENDIF

  (*self.flags).aflag[*,*]=1
  npass=(*self.hmz).npass
  nbad1=0
  nbad2=0

  ;----- Accumulate over all molecules rather than branching on a count. The
  ; fill-value and pointwise-limit tests are ORed across molecules, so a point
  ; is dropped if any molecule is bad there, and the scanwise test is ANDed, so
  ; a scan is dropped only if every molecule exceeds the limit. Both match what
  ; the old one- and two-molecule branches did. -------;
  nmols=(*self.hmz).nmols

  isBad=BYTARR(nlev,npass)
  isOverPoint=BYTARR(nlev,npass)
  isOverScan=BYTARR(nlev,npass)+1B

  FOR m=0,nmols-1 DO BEGIN
    mr=(*self.hmz).mixrat[*,*,m]
    me=(*self.hmz).mixrat_err[*,*,m]
    isBad OR= (mr EQ -999. OR me EQ -999. OR me EQ -888.)
    isOverPoint OR= (me GE self.params.eLim)
    isOverScan AND= (me GE self.params.eLim)
  ENDFOR

  bad1=WHERE(isBad,nbad1)

  IF self.params.cuttype.contains('pointwise',/FOLD_CASE) THEN BEGIN
    bad2=WHERE(isOverPoint,nbad2)
  ENDIF
  IF self.params.cuttype.contains('scanwise',/FOLD_CASE) THEN BEGIN
    jscan=WHERE(isOverScan,njscan)
    IF njscan GT 0 THEN (*self.flags).eflag[jscan/nlev]=0
  ENDIF

  IF nbad1 GT 0 THEN (*self.flags).aflag[bad1]=0
  IF self.params.cuttype.contains('pointwise',/FOLD_CASE) THEN BEGIN
    IF nbad2 GT 0 THEN (*self.flags).aflag[bad2]=0
  ENDIF
  ;stop

  
  self.params.eMod=0

END

PRO ACEData::errorByExtrema

  COMPILE_OPT IDL2

  MESSAGE,'ACEData::errorByExtrema is not implemented.'

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::isOpMaskDefined
;
; PURPOSE:
;       The ACEData::isOpMaskDefined procedure method is a private
;       method and is not intended to be called directly. It declares the operational mask
;       (flags used to indicate to other methods where operations on the data should take place)
;       based on the size and parameters of the given data set.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Added support for Julian Day Calls, 06/30/22
;-

PRO ACEData::isOpMaskDefined

  COMPILE_OPT IDL2

  if self.opMask eq !null then begin
    self.opMask=PTR_NEW(/ALLOCATE)
    C=boolarr(N_ELEMENTS((*self.hmz).alt),(*self.hmz).npass)
    (*self.opMask)=C
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::setOpMask
;
; PURPOSE:
;       The ACEData::setOpMask procedure method is a private
;       method and is not intended to be called directly. It sets
;       the flags for the operation about to take place. The inputs
;       define the spatial and temporal limits of the mask.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Added Support for Julian Day Calls, 06/30/22
;-

PRO ACEData::setOpMask,blats,blons,balts,jdays=jdays,byears=byears,bmonths=bmonths,bdays=bdays,bhours=bhours,sunrise=sunrise,sunset=sunset

  COMPILE_OPT IDL2

  ;----- No longer used. Block selection moved to ACEData::blockPasses plus the
  ; altitude level ranges precomputed in calcBlocks, which do the same job
  ; without materialising a full nlev x npass mask per block. The body below is
  ; left for reference. -------;
  MESSAGE,'ACEData::setOpMask is no longer used. See ACEData::blockPasses.'

  self.isOpMaskDefined
  
;  yearInds=indgen(byears[1]-byears[0]+1)+byears[0]
;  yearIncl=intarr((*self.hmz).npass)
;  for ii=0,n_elements(yearInds)-1 do begin
;    jyr=where((*self.hmz).year eq yearInds[ii])
;    yearIncl[jyr]=1
;  endfor
;
;  if bmonths[1] ge bmonths[0] then begin
;    monthInds=indgen(bmonths[1]-bmonths[0]+1)+bmonths[0]
;    monthIncl=intarr((*self.hmz).npass)
;    for ii=0,n_elements(monthInds)-1 do begin
;      jmth=where((*self.hmz).month eq monthInds[ii])
;      monthIncl[jmth]=1
;      mth_str+=mth_lab[monthInds[ii]-1]
;    endfor
;  endif else begin
;    monthInds=[indgen(12.-bmonths[0]+1)+bmonths[0],indgen(bmonths[1])+1]
;    monthIncl=intarr((*self.hmz).npass)
;    for ii=0,n_elements(monthInds)-1 do begin
;      if monthInds[ii] le 12 then begin
;       jmth=where((*self.hmz).month eq monthInds[ii] AND (*self.hmz).year eq bmonths[0])
;      endif else begin     
;       jmth=where((*self.hmz).month eq monthInds[ii] AND (*self.hmz).year eq bmonths[1])
;      endelse
;      monthIncl[jmth]=1
;    endfor
;  endelse
;
;    if bdays[0] ge bdays[1] then begin
;      
;    endif else begin
;      dayInds=indgen(self.params.daylims[1]-self.params.daylims[0]+1)+self.params.daylims[0]
;    endelse
;    dayIncl=intarr((*self.hmz).npass)
;    for ii=0,n_elements(dayInds)-1 do begin
;      jdy=where((*self.hmz).day eq dayInds[ii])
;      dayIncl[jdy]=1
;    endfor
;
;  if self.params.hourlims[1] ne 0 then begin
;    if self.params.hourlims[0] le self.params.hourlims[1] then begin
;      hourIncl=(((*self.hmz).hour ge self.params.hourlims[0]) AND ((*self.hmz).hour le self.params.hourlims[1]))
;    endif else begin
;      hourIncl=(((*self.hmz).hour le self.params.hourlims[1]) OR ((*self.hmz).hour ge self.params.hourlims[0]))
;    endelse
;  endif else begin
;    hourIncl=intarr((*self.hmz).npass)+1
;  endelse
  IF KEYWORD_SET(jdays) THEN BEGIN
    k=WHERE((*self.hmz).juldays GE jdays[0] AND (*self.hmz).juldays LE jdays[1] AND ((*self.hmz).lat GE blats[0]) AND ((*self.hmz).lat LE blats[1]) AND ((*self.hmz).lon GE blons[0]) AND ((*self.hmz).lon LE blons[1]),nk,complement=j,ncomplement=nj)
  ENDIF ELSE BEGIN
    ;k=WHERE(((*self.hmz).year GE byears[0]) AND ((*self.hmz).year LE byears[1]) AND ((*self.hmz).month GE bmonths[0]) AND ((*self.hmz).month LE bmonths[1]) AND ((*self.hmz).day GE bdays[0]) AND ((*self.hmz).day LE bdays[1]) AND ((*self.hmz).hour GE bhours[0]) AND ((*self.hmz).hour LE bhours[1]) AND ((*self.hmz).lat GE blats[0]) AND ((*self.hmz).lat LE blats[1]) AND ((*self.hmz).lon GE blons[0]) AND ((*self.hmz).lon LE blons[1]),complement=j)
    k=WHERE((*self.flags).tflag AND ((*self.hmz).lat GE blats[0]) AND ((*self.hmz).lat LE blats[1]) AND ((*self.hmz).lon GE blons[0]) AND ((*self.hmz).lon LE blons[1]),nk,complement=j,ncomplement=nj)

  ENDELSE
  IF nk GT 0 THEN (*self.opMask)[*,k]=1
  IF nj GT 0 THEN (*self.opMask)[*,j]=0
  ;deal with altitudes here

  ;----- Level count and grid extent come from the retrieval grid rather than
  ; assuming 150 levels spanning 0-150 km. -------;
  alt=(*self.hmz).alt
  nlev=N_ELEMENTS(alt)

  if balts[0] gt MIN(alt) or balts[1] lt MAX(alt) then begin

    altind0=intarr(nlev)
    altind0[*]=1

    jlo=WHERE(alt le balts[0],njlo)
    if njlo gt 0 then altind0[jlo]=0
    jhi=WHERE(alt ge balts[1],njhi)
    if njhi gt 0 then altind0[jhi]=0

    altind=rebin(altind0,nlev,(*self.hmz).npass)
    
    (*self.opMask)=((*self.opMask) AND altind)
    
    
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::areBlocksDefined
;
; PURPOSE:
;       The ACEData::areBlocksDefined procedure method is a private
;       method and is not intended to be called directly. It checks
;       if the blocks that hold the reduced data have been defined,
;       and if not, declares them.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO ACEData::areBlocksDefined

  COMPILE_OPT IDL2

  IF ~PTR_VALID(self.blocks) THEN self.resetBlocks

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::blockTemplate
;
; PURPOSE:
;       The ACEData::blockTemplate function method is a private method and is
;       not intended to be called directly. It returns the structure that one
;       element of the reduced data array is built from.
;
;       This is the single definition of that structure. areBlocksDefined and
;       resetBlocks previously carried their own copies with different tag
;       lists, so which fields a block had depended on which path had
;       allocated the array.
;
;       orbIDs is deliberately a null pointer here. The caller allocates one
;       heap variable per block after REPLICATE; putting a live pointer in the
;       template would give every block the same one.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 8/4/25
;-

FUNCTION ACEData::blockTemplate

  COMPILE_OPT IDL2

  nmols=(*self.hmz).nmols
  nrat=nmols-1

  ;----- mean, sd and bw carry one entry per molecule; meanRatio, sdRatio and
  ; bwRatio carry one per ratio against the primary. -------;
  D={nmols:nmols,names:(*self.hmz).names, $
     yearlims:[0d0,0d0],monthlims:[0d0,0d0],daylims:[0d0,0d0], $
     hourlims:[0d0,0d0],jullims:[0d0,0d0], $
     latlims:[0d0,0d0],lonlims:[0d0,0d0],altlims:[0d0,0d0], $
     mean:FLTARR(nmols),sd:FLTARR(nmols),bw:FLTARR(5,nmols), $
     meant:0.,sdt:0.,bwt:FLTARR(5), $
     meanp:0.,sdp:0.,bwp:FLTARR(5), $
     meanpt:0.,sdpt:0.,bwpt:FLTARR(5), $
     meant_dmp:0.,sdt_dmp:0.,bwt_dmp:FLTARR(5), $
     meanLat:0.,sdLat:0.,meanLon:0.,sdLon:0., $
     nPts:0L,orbIDs:PTR_NEW()}

  IF nrat GT 0 THEN BEGIN
    D=CREATE_STRUCT(D,'meanRatio',FLTARR(nrat),'sdRatio',FLTARR(nrat), $
                      'bwRatio',FLTARR(5,nrat))
  ENDIF

  ;----- Legacy scalar aliases. These are what existing driver scripts and
  ; interactive sessions read, so they are maintained alongside the arrays for
  ; the first one or two molecules and the first ratio. They are plain scalars,
  ; not one-element arrays, because (*self.blocks).mean0 on a struct array has
  ; to keep returning one value per block. -------;
  D=CREATE_STRUCT(D,'name0',(*self.hmz).names[0],'mean0',0.,'sd0',0.,'bw0',FLTARR(5))

  IF nmols GT 1 THEN BEGIN
    D=CREATE_STRUCT(D,'name1',(*self.hmz).names[1], $
                      'mean1',0.,'sd1',0.,'bw1',FLTARR(5), $
                      'meanRat',0.,'sdRat',0.,'bwRat',FLTARR(5))
  ENDIF

  RETURN,D

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::resetBlocks
;
; PURPOSE:
;       The ACEData::resetBlocks procedure method is a private
;       method and is not intended to be called directly. If the parameters
;       that determine the dimensions of the reduced data are changed (i.e., the user
;       changes the number of latitude bins), then this procedure frees the pointers
;       to the previous reduced data and redeclares the blocks.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO ACEData::resetBlocks

  COMPILE_OPT IDL2

  ;----- Release the per-block heap variables before dropping the array that
  ; holds the pointers to them. -------;
  IF PTR_VALID(self.blocks) THEN BEGIN
    IF N_ELEMENTS(*self.blocks) GT 0 THEN PTR_FREE,(*self.blocks).orbIDs
    PTR_FREE,self.blocks
  ENDIF

  self.blocks=PTR_NEW(/ALLOCATE)
  (*self.blocks)=REPLICATE(self.blockTemplate(),self.params.nLats, $
                           self.params.nLons,self.params.nAlts,self.params.nTimes)

  nblocks=N_ELEMENTS(*self.blocks)
  FOR ib=0,nblocks-1 DO (*self.blocks)[ib].orbIDs=PTR_NEW(/ALLOCATE)

  self.params.bMod=0

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::emptyBlock
;
; PURPOSE:
;       The ACEData::emptyBlock procedure method is a private
;       method and is not intended to be called directly. It marks a single
;       reduced-data block as containing no data, by setting every statistic
;       to NaN and the point count to zero. This is what ACEData::reduceBlock
;       calls when no retrievals fall inside the block, so that empty bins are
;       distinguishable from bins that genuinely averaged to zero.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 8/3/25
;-

PRO ACEData::emptyBlock,i,j,k,l

  COMPILE_OPT IDL2

  nan=!VALUES.F_NAN
  nanbw=FLTARR(5)+nan
  nmols=(*self.hmz).nmols
  nrat=nmols-1

  (*self.blocks)[i,j,k,l].mean=FLTARR(nmols)+nan
  (*self.blocks)[i,j,k,l].sd=FLTARR(nmols)+nan
  (*self.blocks)[i,j,k,l].bw=FLTARR(5,nmols)+nan

  IF nrat GT 0 THEN BEGIN
    (*self.blocks)[i,j,k,l].meanRatio=FLTARR(nrat)+nan
    (*self.blocks)[i,j,k,l].sdRatio=FLTARR(nrat)+nan
    (*self.blocks)[i,j,k,l].bwRatio=FLTARR(5,nrat)+nan
  ENDIF

  (*self.blocks)[i,j,k,l].mean0=nan
  (*self.blocks)[i,j,k,l].sd0=nan
  (*self.blocks)[i,j,k,l].bw0=nanbw
  (*self.blocks)[i,j,k,l].meant=nan
  (*self.blocks)[i,j,k,l].sdt=nan
  (*self.blocks)[i,j,k,l].bwt=nanbw
  (*self.blocks)[i,j,k,l].meanp=nan
  (*self.blocks)[i,j,k,l].sdp=nan
  (*self.blocks)[i,j,k,l].bwp=nanbw
  (*self.blocks)[i,j,k,l].meanpt=nan
  (*self.blocks)[i,j,k,l].sdpt=nan
  (*self.blocks)[i,j,k,l].bwpt=nanbw
  (*self.blocks)[i,j,k,l].meant_dmp=nan
  (*self.blocks)[i,j,k,l].sdt_dmp=nan
  (*self.blocks)[i,j,k,l].bwt_dmp=nanbw
  (*self.blocks)[i,j,k,l].meanLat=nan
  (*self.blocks)[i,j,k,l].sdLat=nan
  (*self.blocks)[i,j,k,l].meanLon=nan
  (*self.blocks)[i,j,k,l].sdLon=nan
  (*self.blocks)[i,j,k,l].nPts=0

  IF nmols GT 1 THEN BEGIN
    (*self.blocks)[i,j,k,l].mean1=nan
    (*self.blocks)[i,j,k,l].sd1=nan
    (*self.blocks)[i,j,k,l].bw1=nanbw
    (*self.blocks)[i,j,k,l].meanRat=nan
    (*self.blocks)[i,j,k,l].sdRat=nan
    (*self.blocks)[i,j,k,l].bwRat=nanbw
  ENDIF

  IF PTR_VALID((*self.blocks)[i,j,k,l].orbIDs) THEN (*(*self.blocks)[i,j,k,l].orbIDs)=!NULL

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::reduceBlocks
;
; PURPOSE:
;       The ACEData::reduceBlocks procedure method is a private
;       method and is not intended to be called directly. This procedure
;       uses the global flags and the operational mask to determine which
;       data should be operated on, then calculates the means, standard deviatons,
;       boxplot, etc. for that block of data. The data is then stored in the
;       appropriate element of the reduced data array.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO ACEData::reduceBlock,i,j,k,l,c,passIdx,levLo,levHi

  COMPILE_OPT IDL2

  ;----- Nothing survived the global flags, so no block has any data. -------;
  IF c.nStart EQ 0 THEN BEGIN
    self.emptyBlock,i,j,k,l
    RETURN
  ENDIF

  ;----- No occultations in this lat/lon/time cell, or no levels in this
  ; altitude bin. -------;
  nSel=N_ELEMENTS(passIdx)
  IF nSel EQ 0 || levLo LT 0 THEN BEGIN
    self.emptyBlock,i,j,k,l
    RETURN
  ENDIF

  ;----- Work on the block's own sub-rectangle. This used to AND a full
  ; nlev x npass operational mask and run WHERE over the result, once per
  ; block, to find at most a few hundred points. -------;
  nSubLev=levHi-levLo+1

  subFlag=c.aflag[levLo:levHi,passIdx]
  kGood=where(subFlag eq 1,ndat)

  ;----- Note that the old ndat=n_elements(kGood) reported 1 for an empty bin,
  ; not 0. -------;
  IF ndat EQ 0 THEN BEGIN
    self.emptyBlock,i,j,k,l
    RETURN
  ENDIF

  ;----- kGood indexes the sub-rectangle, whose first dimension is level, so
  ; the quotient gives the position within passIdx. -------;
  finOrbs=kGood/nSubLev
  selOrbs=passIdx[finOrbs]

  nmols=(*self.hmz).nmols
  nrat=nmols-1

  ;----- One pass over the molecules, and one over the ratios, instead of a
  ; hardcoded primary plus optional second molecule. -------;
  bmean=FLTARR(nmols)
  bsd=FLTARR(nmols)+!VALUES.F_NAN
  bbw=FLTARR(5,nmols)

  FOR m=0,nmols-1 DO BEGIN
    mg=(c.mr[levLo:levHi,passIdx,m])[kGood]
    bmean[m]=mean(mg)
    IF ndat GE 2 THEN bsd[m]=stddev(mg)
    IF ndat GE 5 THEN bbw[*,m]=createboxplotdata(mg)
  ENDFOR

  (*self.blocks)[i,j,k,l].mean=bmean
  (*self.blocks)[i,j,k,l].sd=bsd
  IF ndat GE 5 THEN (*self.blocks)[i,j,k,l].bw=bbw

  IF nrat GT 0 THEN BEGIN
    rmean=FLTARR(nrat)
    rsd=FLTARR(nrat)+!VALUES.F_NAN
    rbw=FLTARR(5,nrat)
    FOR m=0,nrat-1 DO BEGIN
      rg=(c.rat[levLo:levHi,passIdx,m])[kGood]
      rmean[m]=mean(rg)
      IF ndat GE 2 THEN rsd[m]=stddev(rg)
      IF ndat GE 5 THEN rbw[*,m]=createboxplotdata(rg)
    ENDFOR
    (*self.blocks)[i,j,k,l].meanRatio=rmean
    (*self.blocks)[i,j,k,l].sdRatio=rsd
    IF ndat GE 5 THEN (*self.blocks)[i,j,k,l].bwRatio=rbw
  ENDIF

  ;----- Legacy scalar aliases, kept in step with the arrays above. -------;
  (*self.blocks)[i,j,k,l].mean0=bmean[0]
  (*self.blocks)[i,j,k,l].sd0=bsd[0]
  IF ndat GE 5 THEN (*self.blocks)[i,j,k,l].bw0=bbw[*,0]

  IF nmols GT 1 THEN BEGIN
    (*self.blocks)[i,j,k,l].mean1=bmean[1]
    (*self.blocks)[i,j,k,l].sd1=bsd[1]
    (*self.blocks)[i,j,k,l].meanRat=rmean[0]
    (*self.blocks)[i,j,k,l].sdRat=rsd[0]
    IF ndat GE 5 THEN BEGIN
      (*self.blocks)[i,j,k,l].bw1=bbw[*,1]
      (*self.blocks)[i,j,k,l].bwRat=rbw[*,0]
    ENDIF
  ENDIF

  tempg=(c.temp[levLo:levHi,passIdx])[kGood]
  presg=(c.pres[levLo:levHi,passIdx])[kGood]
  ptg=(c.pt[levLo:levHi,passIdx])[kGood]
  tdmpg=(c.tdmp[levLo:levHi,passIdx])[kGood]
  (*self.blocks)[i,j,k,l].meant=mean(tempg)
  (*self.blocks)[i,j,k,l].meanp=mean(presg)
  (*self.blocks)[i,j,k,l].meanpt=mean(ptg)
  (*self.blocks)[i,j,k,l].meant_dmp=mean(tdmpg)

  ;----- STDDEV requires at least two points. -------;
  IF ndat GE 2 THEN BEGIN
    (*self.blocks)[i,j,k,l].sdt=stddev(tempg)
    (*self.blocks)[i,j,k,l].sdp=stddev(presg)
    (*self.blocks)[i,j,k,l].sdpt=stddev(ptg)
    (*self.blocks)[i,j,k,l].sdt_dmp=stddev(tdmpg)
  ENDIF ELSE BEGIN
    (*self.blocks)[i,j,k,l].sdt=!VALUES.F_NAN
    (*self.blocks)[i,j,k,l].sdp=!VALUES.F_NAN
    (*self.blocks)[i,j,k,l].sdpt=!VALUES.F_NAN
    (*self.blocks)[i,j,k,l].sdt_dmp=!VALUES.F_NAN
  ENDELSE

  (*self.blocks)[i,j,k,l].meanLat=mean(c.lat[selOrbs])
  (*self.blocks)[i,j,k,l].meanLon=mean(c.lon[selOrbs])
  (*self.blocks)[i,j,k,l].nPts=ndat
  ;stop
  if ndat ge 5 then begin
    (*self.blocks)[i,j,k,l].bwt=createboxplotdata(tempg)
    (*self.blocks)[i,j,k,l].bwp=createboxplotdata(presg)
    (*self.blocks)[i,j,k,l].bwpt=createboxplotdata(ptg)
    IF self.params.isDMP THEN (*self.blocks)[i,j,k,l].bwt_dmp=createboxplotdata(tdmpg)
  endif
  ;stop
  ;(*(*self.blocks)[i,j,k,l].orbIDs)=(*self.hmz).orbit[startOrbs[finOrbs]]
  (*(*self.blocks)[i,j,k,l].orbIDs)=c.startOrbs[selOrbs]

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::blockPasses
;
; PURPOSE:
;       The ACEData::blockPasses function method is a private method and is not
;       intended to be called directly. It returns the positions, in the space
;       of the ACEData::reduceCache arrays, of the occultations that fall inside
;       one lat/lon/time cell.
;
;       Together with reduceCache this replaces the operational mask. The old
;       arrangement kept a full nlev x npass byte mask, rewrote all of it for
;       every block in setOpMask, then ANDed a full-size slice of it against the
;       altitude flags and ran WHERE over the result inside reduceBlock, all to
;       locate at most a few hundred points.
;
;       Block membership factorises instead: latitude, longitude and time select
;       occultations, which is one dimension, and altitude selects a contiguous
;       run of levels, because the retrieval grid is monotonic. So a block is a
;       sub-rectangle, and reduceBlock only has to touch that.
;
;       Latitude and longitude do not depend on the altitude bin, so calcBlocks
;       calls this once per lat/lon/time cell and reuses the result down the
;       whole altitude column.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 8/4/25
;-

FUNCTION ACEData::blockPasses,c,blats,blons,jdays=jdays

  COMPILE_OPT IDL2

  IF c.nStart EQ 0 THEN RETURN,!NULL

  ;----- Bounds are inclusive at both ends, which is what the operational mask
  ; did, so an occultation sitting exactly on an interior bin edge belongs to
  ; both neighbouring bins. That is preserved deliberately. HISTOGRAM-based
  ; binning would be marginally faster but assigns each occultation to exactly
  ; one bin, which would quietly move numbers in existing analyses. -------;
  IF N_ELEMENTS(jdays) EQ 2 THEN BEGIN
    ksel=WHERE(c.jul GE jdays[0] AND c.jul LE jdays[1] AND $
               c.lat GE blats[0] AND c.lat LE blats[1] AND $
               c.lon GE blons[0] AND c.lon LE blons[1],nsel)
  ENDIF ELSE BEGIN
    ;----- The global time flag is already folded into the cache, which is what
    ; the old non-solidBlock branch of setOpMask re-applied here. -------;
    ksel=WHERE(c.lat GE blats[0] AND c.lat LE blats[1] AND $
               c.lon GE blons[0] AND c.lon LE blons[1],nsel)
  ENDELSE

  IF nsel EQ 0 THEN RETURN,!NULL

  RETURN,ksel

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::reduceCache
;
; PURPOSE:
;       The ACEData::reduceCache function method is a private method and is not
;       intended to be called directly. It returns the quantities reduceBlock
;       needs that do not vary from block to block, so calcBlocks can compute
;       them once per pass rather than once per bin.
;
;       That is the global flag intersection, and the slices of the harmonized
;       arrays taken at the occultations those flags admit. reduceBlock used to
;       rebuild all of them on every call, which meant copying five or seven
;       full-size arrays per bin.
;
;       See the header above ACEData::blockPasses for how a block is selected.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 8/4/25
;-

FUNCTION ACEData::reduceCache

  COMPILE_OPT IDL2

  goodOrbs=((*self.flags).sflag AND (*self.flags).tflag AND (*self.flags).eflag)
  startOrbs=WHERE(goodOrbs EQ 1,nStart)

  IF nStart EQ 0 THEN RETURN,{nStart:0L}

  c={nStart:nStart,startOrbs:startOrbs, $
     aflag:(*self.flags).aflag[*,startOrbs], $
     mr:(*self.hmz).mixrat[*,startOrbs,*], $
     temp:(*self.hmz).temp[*,startOrbs], $
     pres:(*self.hmz).pres[*,startOrbs], $
     pt:(*self.hmz).pt[*,startOrbs], $
     tdmp:(*self.hmz).t_dmp[*,startOrbs], $
     lat:(*self.hmz).lat[startOrbs], $
     lon:(*self.hmz).lon[startOrbs], $
     jul:(*self.hmz).juldays[startOrbs]}

  IF (*self.hmz).nmols GT 1 THEN BEGIN
    c=CREATE_STRUCT(c,'rat',(*self.hmz).ratio[*,startOrbs,*])
  ENDIF

  RETURN,c

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::calcBlocks
;
; PURPOSE:
; 
;       The ACEData::calcBlocks calculates the parameters of interest for each lat/lot/alt/time
;       bin as defined by the flags and operational mask. Currently this includes the mean, sd,
;       and boxplot of each molecule and ratio, as well as temperature and pressure.
;
;
; CALLING SEQUENCE:
;
;       myData.calcBlocks
;
; INPUTS: NONE
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: NONE
;
; EXAMPLE:
;
;       myData.calcBlocks
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Added support for continuous calculations from julian day, 7/1/22
;-

PRO ACEData::calcBlocks,useParam=useParam

  COMPILE_OPT IDL2

  IF ~PTR_VALID(self.hmz) || N_ELEMENTS(*self.hmz) EQ 0 THEN BEGIN
    PRINT,'ERROR: NO HARMONIZED DATA. Call loadFiles then harmonize first.'
    RETURN
  ENDIF

  IF self.params.nLats LT 1 OR self.params.nLons LT 1 OR $
     self.params.nAlts LT 1 OR self.params.nTimes LT 1 THEN BEGIN
    PRINT,'ERROR: nLats, nLons, nAlts and nTimes must all be at least 1.'
    RETURN
  ENDIF

  daysInMonth=[31,28,31,30,31,30,31,31,30,31,30,31]

  ;----- Recompute the flags when a parameter changed, and also when they have
  ; never been computed at all, which reduceCache below depends on. -------;
  flagsStale=(self.params.sMod EQ 1 OR self.params.tMod EQ 1 OR $
              self.params.eMod EQ 1 OR self.params.bMod EQ 1)
  IF ~PTR_VALID(self.flags) THEN BEGIN
    flagsStale=1B
  ENDIF ELSE IF N_ELEMENTS(*self.flags) EQ 0 THEN BEGIN
    flagsStale=1B
  ENDIF

  IF flagsStale THEN BEGIN
    print,'Parameters have Changed, resetting blocks'
    self.setSpaceFlag
    self.setTimeFlag
    self.setErrorFlag
    self.resetBlocks
  ENDIF

  ;----- Built once for the whole pass. See ACEData::reduceCache. -------;
  rc=self.reduceCache()

  latsran=(self.params.latlims[1]-self.params.latlims[0])/self.params.nlats*dindgen(self.params.nlats+1)+self.params.latlims[0]
  longsran=(self.params.lonlims[1]-self.params.lonlims[0])/self.params.nlons*dindgen(self.params.nlons+1)+self.params.lonlims[0]
  altsran=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts+1)+self.params.altlims[0]
  
  dayLast=julday(self.params.monthlims[1],self.params.daylims[1],self.params.yearlims[1],self.params.hourlims[1],0,0)
  dayFirst=julday(self.params.monthlims[0],self.params.daylims[0],self.params.yearlims[0],self.params.hourlims[0],0,0)
  daysran=(dayLast-dayFirst)/self.params.nTimes*dindgen(self.params.nTimes+1)+dayFirst

  ;----- Each altitude bin maps to a contiguous run of retrieval levels, which
  ; relies on the altitude grid being monotonic. Say so rather than mis-slicing
  ; silently if that ever stops holding. -------;
  alt=(*self.hmz).alt
  nlev=N_ELEMENTS(alt)
  altMin=MIN(alt)
  altMax=MAX(alt)
  IF nlev GT 1 THEN BEGIN
    IF TOTAL(alt[1:*] LE alt[0:nlev-2]) GT 0 THEN $
      PRINT,'WARNING: ALTITUDE GRID IS NOT MONOTONIC; BLOCK LEVEL RANGES MAY BE WRONG'
  ENDIF

  levLo=LONARR(self.params.nAlts)
  levHi=LONARR(self.params.nAlts)
  FOR ka=0,self.params.nAlts-1 DO BEGIN
    alims=altsran[ka:ka+1]
    ;----- Reproduces the old setOpMask altitude cut exactly, including that it
    ; applied no cut at all when the bin spanned the whole grid. -------;
    IF alims[0] GT altMin OR alims[1] LT altMax THEN BEGIN
      jlev=WHERE(alt GT alims[0] AND alt LT alims[1],nlevsel)
    ENDIF ELSE BEGIN
      jlev=LINDGEN(nlev)
      nlevsel=nlev
    ENDELSE
    IF nlevsel GT 0 THEN BEGIN
      levLo[ka]=jlev[0]
      levHi[ka]=jlev[nlevsel-1]
    ENDIF ELSE BEGIN
      levLo[ka]=-1
      levHi[ka]=-1
    ENDELSE
  ENDFOR

  ;stop
  IF self.params.solidBlock EQ 1 THEN BEGIN
    ;----- Altitude is the innermost loop so that one pass selection is shared
    ; by the whole altitude column above a lat/lon/time cell. -------;
    for l=0,self.params.nTimes-1 do begin
      jlims=daysran[l:l+1]
      for j=0,self.params.nLons-1 do begin
        lonlims=longsran[j:j+1]
        for i=0,self.params.nLats-1 do begin
          latlims=latsran[i:i+1]
          passIdx=self.blockPasses(rc,latlims,lonlims,jdays=jlims)
          for k=0,self.params.nAlts-1 do begin
            alims=altsran[k:k+1]
            self.reduceBlock,i,j,k,l,rc,passIdx,levLo[k],levHi[k]
            (*self.blocks)[i,j,k,l].yearlims=[-1,-1]
            (*self.blocks)[i,j,k,l].monthlims=[-1,-1]
            (*self.blocks)[i,j,k,l].daylims=[-1,-1]
            (*self.blocks)[i,j,k,l].hourlims=[-1,-1]
            (*self.blocks)[i,j,k,l].jullims=jlims
            (*self.blocks)[i,j,k,l].latlims=latlims
            (*self.blocks)[i,j,k,l].lonlims=lonlims
            (*self.blocks)[i,j,k,l].altlims=alims

          endfor
        endfor
      endfor
    endfor
    ;stop
  ENDIF ELSE BEGIN
    timesran=fltarr(n_elements(daysran),4)
    caldat,daysran,months,days,years,hours,mins,sec
    timesran[0,0]=years
    timesran[0,1]=months
    timesran[0,2]=days
    timesran[0,3]=hours
    ;stop
    ;----- The per-bin year/month/day/hour limits below are recorded on each
    ; block but do not reach the selection: as in the previous implementation,
    ; the global time flag is what restricts occultations in this branch, so
    ; every time bin sees the same set. That is the known nTimes > 1 defect and
    ; is unchanged here. -------;
    for l=0,self.params.nTimes-1 do begin
      ylims=[timesran[l,0],timesran[l+1,0]]
      mlims=[timesran[l,1],timesran[l+1,1]]
      dlims=[timesran[l,2],timesran[l+1,2]]
      hlims=[timesran[l,3],timesran[l+1,3]]
      for j=0,self.params.nLons-1 do begin
        lonlims=longsran[j:j+1]
        for i=0,self.params.nLats-1 do begin
          latlims=latsran[i:i+1]
          passIdx=self.blockPasses(rc,latlims,lonlims)
          for k=0,self.params.nAlts-1 do begin
            alims=altsran[k:k+1]
            self.reduceBlock,i,j,k,l,rc,passIdx,levLo[k],levHi[k]
            (*self.blocks)[i,j,k,l].yearlims=ylims
            (*self.blocks)[i,j,k,l].monthlims=mlims
            (*self.blocks)[i,j,k,l].daylims=dlims
            (*self.blocks)[i,j,k,l].hourlims=hlims
            (*self.blocks)[i,j,k,l].jullims=[-1,-1]
            (*self.blocks)[i,j,k,l].latlims=latlims
            (*self.blocks)[i,j,k,l].lonlims=lonlims
            (*self.blocks)[i,j,k,l].altlims=alims
            ;stop

          endfor
        endfor
      endfor
    endfor
  ENDELSE
END

FUNCTION ACEData::makeHist,nbins,range,oneBin=oneBin

  COMPILE_OPT IDL2

  ;----- Unfinished: the oneBin branch was never written, and the branch below
  ; stops inside its own loop. Left in place for reference rather than deleted.
  ; -------;
  MESSAGE,'ACEData::makeHist is not implemented.'

  IF keyword_set(oneBin) THEN BEGIN

  ENDIF ELSE BEGIN
    nBlocks=self.params.nLats*self.params.nLons*self.params.nAlts*self.params.nTimes
    FOR i=0,nBlocks-1 DO BEGIN
      norbs=(*self.blocks)[i].npts
      dex=(*(*self.blocks)[i].orbIDs)
      orbs=self.getOcc(dex)
      ;FOR is=0,norbs-1 DO BEGIN
        ;jk=where((*self.flags).aflag[*,dex[is]] eq 1)
        altInd=(*self.blocks)[i].altLims[0]
        hst0=histogram(reform((orbs.ratio[altInd,*]-1d0)*1d3,norbs),min=range[0],max=range[1],nbins=nbins,locations=xbin)
        p=plot(xbin,hst0)
      ;ENDFOR
      
    ENDFOR
  ENDELSE

END

PRO ACEData::setAxes,xaxis,yaxis

  COMPILE_OPT IDL2

  MESSAGE,'ACEData::setAxes is not implemented.'

END


;+
; =============================================================
;
; METHODNAME:
;       ACEData::makeContour
;
; PURPOSE:
;
;       The ACEData::makeContour procedure method makes a contour plot
;       of the reduced data at a specified altitude level.
;
;
; CALLING SEQUENCE:
;
;       myData.makeContour,ratPlot,mol0Plot,aLevels
;
; INPUTS: 
;
;       ratPlot is a variable that will contain the ratio plot generated by the method
;       
;       mol0Plot is a variable that will contain the plot generated by the method
;       for the primary molecule.
;       
;       aLevels are the altitude levels at which the contour plots will be generated. Currently
;       only one level at a time is supported.
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS:
; 
;       Invoking the useMedian parameter will result in the plots being generated using
;       the median values as calculated by calcBlocks. The default behavior is for the mean 
;       to be used.
;
; EXAMPLE:
;
;       myData.makeContour,rat,mol0,16.5
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;   Fixed makeContour to use interpolate instead of griddata. Needs testing on non-linux platforms. 7/13/22
;
;-

PRO ACEData::makeContour,ratPlot,mol0Plot,aLevels,useMedian=useMedian, $
                         molIndex=molIndex,ratIndex=ratIndex

  COMPILE_OPT IDL2

  IF ~PTR_VALID(self.blocks) || N_ELEMENTS(*self.blocks) EQ 0 THEN BEGIN
    PRINT,'ERROR: NO REDUCED DATA. Call calcBlocks first.'
    RETURN
  ENDIF

  ;----- Which molecule and which ratio to draw, for the case of more than two
  ; molecules. Both default to the first, so two-molecule calls are unaffected.
  ; -------;
  IF N_ELEMENTS(molIndex) EQ 0 THEN molIndex=0
  IF N_ELEMENTS(ratIndex) EQ 0 THEN ratIndex=0

  IF molIndex LT 0 OR molIndex GE (*self.hmz).nmols THEN BEGIN
    PRINT,'ERROR: molIndex OUT OF RANGE. '+STRTRIM((*self.hmz).nmols,2)+' molecule(s) harmonized.'
    RETURN
  ENDIF
  IF (*self.hmz).nmols GT 1 THEN BEGIN
    IF ratIndex LT 0 OR ratIndex GE (*self.hmz).nmols-1 THEN BEGIN
      PRINT,'ERROR: ratIndex OUT OF RANGE. '+STRTRIM((*self.hmz).nmols-1,2)+' ratio(s) available.'
      RETURN
    ENDIF
  ENDIF

  alt_map=aLevels             ;currently only one level at a time is supported

  ct=colortable(72,/reverse)
  ;ct=colortable(25)
  
  ;-------- This block defines a few variables for the plots. The case structure
  ;defines the color scale base on altitude, and could probably use further tweaks. --------;
  posit=[0.05,0.1,0.87,0.96]
  dim=[1280,580]
  case (*self.hmz).names[0] of
    'H2O': BEGIN
      nlevs=16d0
      case alt_map of
        14:mol0range=dindgen(nlevs+1)/nlevs*20d0+2.5d0
        15:mol0range=dindgen(nlevs+1)/nlevs*12d0+3.2d0
        16:mol0range=dindgen(nlevs+1)/nlevs*4.8d0+3.2d0
        17:mol0range=dindgen(nlevs+1)/nlevs*4.8d0+3.2d0
        18:mol0range=dindgen(nlevs+1)/nlevs*4.8d0+3.2d0
        19:mol0range=dindgen(nlevs+1)/nlevs*4.8d0+3.2d0
        20:mol0range=dindgen(nlevs+1)/nlevs*4.8d0+3.2d0
        else:mol0range=dindgen(nlevs+1)/nlevs*4.8d0+3.2d0
      endcase
    END
    'CO2':BEGIN
      nlevs=16d0
      mol0range=dindgen(nlevs+1)/nlevs*35d0+390d0
    END
    ELSE: BEGIN
      ;----- Any other molecule: span whatever range is actually present in the
      ; blocks. This CASE previously had no ELSE, so loading anything but H2O
      ; or CO2 crashed here. -------;
      nlevs=16d0
      blockMean=(*self.blocks).mean[molIndex]
      kfin=WHERE(FINITE(blockMean),nfin)
      IF nfin EQ 0 THEN MESSAGE,'No finite block means to contour.'
      lo=MIN(blockMean[kfin])*1d6
      hi=MAX(blockMean[kfin])*1d6
      IF hi LE lo THEN hi=lo+1d0
      mol0range=dindgen(nlevs+1)/nlevs*(hi-lo)+lo
    END

  endcase

  ;-------- This section plots the ratio and primary molecule data. If savPlotRat or savPlot0 are set, then the plots are saved
  ;directly to file.
  
  lats=(dindgen(self.params.nLats)+0.5d0)/self.params.nLats*(self.params.latLims[1]-self.params.latLims[0])+self.params.latLims[0]
  lons=(dindgen(self.params.nLons)+0.5d0)/self.params.nLons*(self.params.lonLims[1]-self.params.lonLims[0])+self.params.lonLims[0]
  
  if (*self.hmz).nmols gt 1 then begin
    if keyword_set(useMedian) then rats=(*self.blocks).bwRatio[2,ratIndex] else rats=(*self.blocks).meanRatio[ratIndex]
    if keyword_set(savPlotRat) then mapRat=map('Geographic',center_longitude=180,limit=[-75,0,75,360],dimensions=dim,position=posit,font_size=18,/buffer) else mapRat=map('Geographic',center_longitude=180,limit=[-75,0,75,360],dimensions=dim,position=posit,font_size=18)
    ;stop
    ratios=interpolate(transpose((rats-1d0)*1d3),dindgen(self.params.nLons),dindgen(self.params.nLats),/grid)
    
    ;cc=contour((rats-1d0)*1d3,(*self.blocks).meanLon,(*self.blocks).meanLat,c_value=([dindgen(17)]/16d0*.2d0-.76d0)*1d3,rgb_table=ct,/fill,/overplot)
    cc=contour(ratios,lons,lats,c_value=([dindgen(21)]/20d0*.30d0-.84d0)*1d3,rgb_table=ct,/fill,/overplot)
    ;mc=mapcontinents(/continents,limit=[-65,-180,65,180])
    mc=mapcontinents(/continents,limit=[-75,0,75,360])
    grid=mapRat.mapgrid
    grid.linestyle='dotted'
    grid.font_size=16
    grid.label_position=0
    ;----- Label the ratio with delta notation only for the H2O/HDO pair. -----;
    IF (*self.hmz).names[0] EQ 'H2O' AND (*self.hmz).names[ratIndex+1] EQ 'HDO' THEN BEGIN
      ratLabel='$\delta D$'
    ENDIF ELSE BEGIN
      ratLabel=(*self.hmz).names[ratIndex+1]+'/'+(*self.hmz).names[0]
    ENDELSE
    mapRat.title=self.params.mth_str+' '+ratLabel+', '+string(alt_map,format='(F4.1)')+' km'
      mapRat.font_size=18
    cb=colorbar(target=cc,title='$\delta D (‰)$',range=[-.840d0,-.540d0]*1d3,orientation=1,textpos=1,font_size=14)
      ;PRINT,'MADE CB'
    ratPlot=mapRat
    if keyword_set(savPlotRat) then begin
      ratPlot.save,savPlotRat,resolution=300
      ratPlot.close
    endif
    ;stop
  endif
  if keyword_set(savPlot0) then map0=map('Geographic',center_longitude=180,limit=[-75,0,75,360],dimensions=dim,position=posit,font_size=18,/buffer) else map0=map('Geographic',center_longitude=180,limit=[-75,0,75,360],dimensions=dim,position=posit,font_size=18)
  if keyword_set(useMedian) then mol0=(*self.blocks).bw[2,molIndex] else mol0=(*self.blocks).mean[molIndex]
  mol0s=interpolate(transpose(mol0*1d6),dindgen(self.params.nLons),dindgen(self.params.nLats),/grid)
  ;cc=contour(mol0*1d6,(*self.blocks).meanLon,(*self.blocks).meanLat,c_value=mol0range,rgb_table=ct,/fill,/overplot)
  cc=contour(mol0s,lons,lats,c_value=mol0range,rgb_table=ct,/fill,/overplot)
  ;mc=mapcontinents(/continents,limit=[-65,-180,65,180])
  mc=mapcontinents(/continents,limit=[-75,0,75,360])
  grid=map0.mapgrid
  grid.linestyle='dotted'
  grid.font_size=16
  grid.label_position=0
  map0.title=self.params.mth_str+' '+(*self.hmz).names[0]+', '+string(alt_map,format='(F4.1)')+' km'
    map0.font_size=18
  cb=colorbar(target=cc,title=(*self.hmz).names[0]+' (ppm)',range=[mol0range[0],mol0range[n_elements(mol0range)-1]],orientation=1,textpos=1,font_size=14)
  mol0Plot=map0

  if keyword_set(savPlot0) then begin
    mol0Plot.save,savPlot0,resolution=300
    mol0Plot.close
  endif

end


;+
; =============================================================
;
; METHODNAME:
;       ACEData::makeMeridianPlot
;
; PURPOSE:
;
;       The ACEData::makeMeridianPlot procedure method makes a meridional plot
;       of the reduced data between the specified latitudes.
;
;
; CALLING SEQUENCE:
;
;       myData.makeMeridianPlot,ratPlot,mol0Plot
;
; INPUTS:
;
;       ratPlot is a variable that will contain the ratio plot generated by the method
;
;       mol0Plot is a variable that will contain the plot generated by the method
;       for the primary molecule.
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS:
;
;       Invoking the useMedian parameter will result in the plots being generated using
;       the median values as calculated by calcBlocks. The default behavior is for the mean
;       to be used.
;
; EXAMPLE:
;
;       myData.makeMeridianPlot,rat,mol0
;       
; TODO: Add ability to handle blocks with longitudinal subsets.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/22/21
;-

PRO ACEData::makeMeridianPlot,ratPlot,mol0Plot,useMedian=useMedian,lats=lats,alts=alts,title=title,dimensions=dimensions,position=position,fs=fs, $
                              molIndex=molIndex,ratIndex=ratIndex

  COMPILE_OPT IDL2

  IF N_ELEMENTS(molIndex) EQ 0 THEN molIndex=0
  IF N_ELEMENTS(ratIndex) EQ 0 THEN ratIndex=0

  IF NOT KEYWORD_SET(title) THEN title=''
  IF NOT KEYWORD_SET(dimensions) THEN dims=[1024,1024] ELSE dims=dimensions
  IF NOT KEYWORD_SET(position) THEN pos=[0.1,0.1,0.85,0.9] ELSE pos=position
  IF NOT KEYWORD_SET(fs) THEN fs=20 ELSE fs=fs
  
  ct=colortable(72,/reverse)

  if keyword_set(useMedian) then rats=(*self.blocks).bwRatio[2,ratIndex] else rats=(*self.blocks).meanRatio[ratIndex]
  n_longs=self.params.nLons
  alts=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts)+self.params.altlims[0]+.5
  alts=transpose(rebin(alts,self.params.nAlts,self.params.nLats))
  ;lats=(self.params.latlims[1]-self.params.latlims[0])/self.params.nLats*dindgen(self.params.nLats)+self.params.latlims[0]
  cbrange=([dindgen(21)]/20d0*0.46d0-.76d0)*1d3

  for i=0,n_longs-1 do begin
    ratios=interpolate(reform((rats[*,i,*]-1d0)*1d3),dindgen(self.params.nLats),dindgen(self.params.nAlts),/grid)
    lats=reform((*self.blocks)[*,i,*].meanLat,self.params.nLats,self.params.nAlts)
    ratPlot=contour(reform((rats[*,i,*]-1d0)*1d3,self.params.nLats,self.params.nAlts),lats,alts,c_value=cbrange,xrange=latlims,yrange=[5,40],font_size=fs,axis_style=2,xthick=2,ythick=2,xtitle='Latitude',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=dims,position=pos,/fill,title=title)
    ;cb=colorbar(target=ratPlot,title='$\delta D (‰)$',range=[-.760,-.3]*1d3,orientation=1,textpos=1,font_size=14,position=[.88,.25,.91,.75])
    

  endfor

  if keyword_set(useMedian) then mol0=(*self.blocks).bw[2,molIndex] else mol0=(*self.blocks).mean[molIndex]
  n_longs=self.params.nLons
  alts=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts)+self.params.altlims[0]+.5
  alts=transpose(rebin(alts,self.params.nAlts,self.params.nLats))
  ;lats=(self.params.latlims[1]-self.params.latlims[0])/self.params.nLats*dindgen(self.params.nLats)+self.params.latlims[0]
  cbrange=[dindgen(21)]/20d0*12d0+2d0
  for i=0,n_longs-1 do begin
    mol0s=interpolate(reform(mol0[*,i,*]*5d10),dindgen(self.params.nLats),dindgen(self.params.nAlts),/grid)
    lats=reform((*self.blocks)[*,i,*].meanLat,self.params.nLats,self.params.nAlts)
    mol0Plot=contour(mol0s,lats,alts,c_value=cbrange,xrange=latlims,yrange=[5,40],font_size=fs,axis_style=2,xthick=2,ythick=2,xtitle='Latitude',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=dims,position=pos,/fill,title=title)
    ;cb=colorbar(target=mol0Plot,title='$H_2O (ppm)$',range=[3,15],orientation=1,textpos=1,font_size=14,position=[.88,.25,.91,.75])
    print,fs
  endfor
  ;stop
END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::makeTimeSeries
;
; PURPOSE:
;
;       The ACEData::makeTimeSeries method plots a lat/lon/height box over a specified set of times
;
;
; CALLING SEQUENCE:
;
;       myData.makeTimeSeries,ratPlot,mol0Plot
;
; INPUTS:
;
;       ratPlot is a variable that will contain the ratio plot generated by the method
;
;       mol0Plot is a variable that will contain the plot generated by the method
;       for the primary molecule.
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS:
;
;       Invoking the useMedian parameter will result in the plots being generated using
;       the median values as calculated by calcBlocks. The default behavior is for the mean
;       to be used.
;
; EXAMPLE:
;
;       myData.makeTimeSeries,rat,mol0
;
; 
; 
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 6/30/22
;-

PRO ACEData::makeTimeSeries,ratPlot,mol0Plot,useMedian=useMedian,ratIndex=ratIndex

  COMPILE_OPT IDL2

  IF N_ELEMENTS(ratIndex) EQ 0 THEN ratIndex=0

  ct=colortable(72,/reverse)

  if keyword_set(useMedian) then rats=(*self.blocks).bwRatio[2,ratIndex] else rats=(*self.blocks).meanRatio[ratIndex]
  n_longs=self.params.nLons
  alts=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts)+self.params.altlims[0]+.5
  alts=transpose(rebin(alts,self.params.nAlts,self.params.nLats))
  ;lats=(self.params.latlims[1]-self.params.latlims[0])/self.params.nLats*dindgen(self.params.nLats)+self.params.latlims[0]
  cbrange=([dindgen(22)]/21d0*0.42d0-.76d0)*1d3
  for i=0,n_longs-1 do begin
    lats=reform((*self.blocks)[*,i,*].meanLat,self.params.nLats,self.params.nAlts)
    cc=contour((reform(rats[*,i,*],self.params.nLats,self.params.nAlts)-1)*1d3,lats,alts,c_value=cbrange,xrange=latlims,yrange=[8,40],font_size=18,xtitle='Latitude',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=[1024,1024],position=[0.1,0.1,0.85,0.9],/fill)
    cb=colorbar(target=cc,title='$\delta D (‰)$',range=[-.760,-.495]*1d3,orientation=1,textpos=1,font_size=14)
    ;stop
  endfor

END


PRO ACEData::makeZonalPlot,mols,prof,xaxis,yaxis,yearlims=yearlims,monthlims=monthlims,hourlims=hourlims,latlims=latlims,lonlims=lonlims,n_lats=n_lats,n_longs=n_longs,ratio=ratio,ratlims=ratlims,raw=raw,acecut=acecut,fullcut=fullcut

  COMPILE_OPT IDL2

  MESSAGE,'ACEData::makeZonalPlot is not implemented.'

END

PRO ACEData::makeHeightTimePlot,mol0,rat0,alt_map=alt_map,yearlims=yearlims,monthlims=monthlims,hourlims=hourlims,latlims=latlims,lonlims=lonlims,n_lats=n_lats,n_longs=n_longs,ratio=ratio,ratlims=ratlims,raw=raw,acecut=acecut,fullcut=fullcut,xrange=xrange, $
                                  molIndex=molIndex,ratIndex=ratIndex

  COMPILE_OPT IDL2

  IF NOT KEYWORD_SET(xrange) THEN xrange=[]
  IF N_ELEMENTS(molIndex) EQ 0 THEN molIndex=0
  IF N_ELEMENTS(ratIndex) EQ 0 THEN ratIndex=0
  
  dayLast=julday(self.params.monthlims[1],self.params.daylims[1],self.params.yearlims[1],self.params.hourlims[1],0,0)
  dayFirst=julday(self.params.monthlims[0],self.params.daylims[0],self.params.yearlims[0],self.params.hourlims[0],0,0)
  daysran=(dayLast-dayFirst)/self.params.nTimes*dindgen(self.params.nTimes+1)+dayFirst
  jdays=((daysran+shift(daysran,-1))/2d0)[0:self.params.nTimes]
  times=DBLARR(self.params.nTimes)
  FOR i=0,self.params.nTimes-1 DO BEGIN
    tt=date_conv(jdays[i],'R')
    times[i]=FLOAT(STRMID(tt,7,4))+(FLOAT(STRMID(tt,11))-1.)/365.25
    ;stop
  ENDFOR
  ;stop
  ct=colortable(72,/reverse)
  
  ntimes=self.params.nTimes
  nalts=self.params.nAlts
  alts=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts)+self.params.altlims[0]+.5
  alts=transpose(rebin(alts,self.params.nAlts,self.params.nTimes))
  ;stop
  ;times=dindgen(self.params.nTimes)
  times=rebin(times,self.params.nTimes,self.params.nAlts)
  mol0=contour(transpose(reform((*self.blocks).bw[2,molIndex],nAlts,nTimes))*1d6,times,alts,c_value=[dindgen(22)]*.25+3d0,xrange=xrange,yrange=[10,40],font_size=18,axis_style=2,xthick=2,ythick=2,xtitle='Year',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=[1024,1024],position=[0.1,0.1,0.85,0.9],/fill)
  cb=colorbar(target=mol0,title='$H_2O (ppm)$',range=[3,8.25],orientation=1,textpos=1,font_size=14)

  rat0=contour((transpose(reform((*self.blocks).bwRatio[2,ratIndex],nAlts,nTimes))-1d0)*1d3,times,alts,c_value=([dindgen(22)]/21d0*0.63d0-.76d0)*1d3,xrange=xrange,yrange=[10,40],font_size=18,axis_style=2,xthick=2,ythick=2,xtitle='Year',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=[1024,1024],position=[0.1,0.1,0.85,0.9],/fill)
  cb=colorbar(target=rat0,title='$\delta D (‰)$',range=[-.760,-.130]*1d3,orientation=1,textpos=1,font_size=14)
  ;stop
END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::selParse
;
; PURPOSE:
;       The ACEData::selParse procedure method is a private method and is not
;       intended to be called directly. It maps a plotting selector name onto a
;       kind and an index, so that the plotting methods extend to any number of
;       harmonized molecules without a branch per molecule.
;
;           'mol0', 'mol1', 'mol2', ...   kind 'mol', index as written
;           'delta', 'delta0', 'delta1'   kind 'rat', index as written, 0 if bare
;           anything else                 kind is the name itself, index -1
;
;       So the selectors that existed before keep their meaning: 'mol0' and
;       'mol1' are the first two molecules and a bare 'delta' is the first
;       ratio against the primary.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 8/4/25
;-

PRO ACEData::selParse,sname,kind,index

  COMPILE_OPT IDL2

  index=-1L
  kind=sname

  mm=STREGEX(sname,'^mol([0-9]+)$',/SUBEXPR,/EXTRACT)
  IF mm[0] NE '' THEN BEGIN
    kind='mol'
    index=LONG(mm[1])
    RETURN
  ENDIF

  mm=STREGEX(sname,'^delta([0-9]*)$',/SUBEXPR,/EXTRACT)
  IF mm[0] NE '' THEN BEGIN
    kind='rat'
    IF mm[1] EQ '' THEN index=0L ELSE index=LONG(mm[1])
    RETURN
  ENDIF

END

pro ACEData::makeMeanSD,xname,yname,spag,savPlot=savPlot,spaghetti=spaghetti,colorby=colorby,outcols=outcols,title=title,nomean=nomean

  COMPILE_OPT IDL2

  IF ~PTR_VALID(self.blocks) || N_ELEMENTS(*self.blocks) EQ 0 THEN BEGIN
    PRINT,'ERROR: NO REDUCED DATA. Call calcBlocks first.'
    RETURN
  ENDIF

  ct=colortable(72,/reverse)
  IF NOT KEYWORD_SET(title) THEN title=''

  nmols=(*self.hmz).nmols
  names=(*self.hmz).names

  self.selParse,xname,xkind,xidx
  self.selParse,yname,ykind,yidx

  IF xkind EQ 'mol' && (xidx GE nmols) THEN $
    MESSAGE,'xname '+xname+' but only '+STRTRIM(nmols,2)+' molecule(s) harmonized.'
  IF ykind EQ 'mol' && (yidx GE nmols) THEN $
    MESSAGE,'yname '+yname+' but only '+STRTRIM(nmols,2)+' molecule(s) harmonized.'
  IF xkind EQ 'rat' && (xidx GE nmols-1) THEN $
    MESSAGE,'xname '+xname+' but only '+STRTRIM(nmols-1,2)+' ratio(s) available.'
  IF ykind EQ 'rat' && (yidx GE nmols-1) THEN $
    MESSAGE,'yname '+yname+' but only '+STRTRIM(nmols-1,2)+' ratio(s) available.'

      CASE ykind OF
        'alt':BEGIN
          ymean=self.params.altlims[0]+dindgen(self.params.nAlts)+0.5d0
          ytitle='Altitude (km)'
          yrange=self.params.altlims
        END
        'pt':BEGIN
          ymean=reform((*self.blocks).meanpt,self.params.nAlts)
          ytitle='$\Theta (K)$'
          yrange=[300,500]
        END
        'mol':BEGIN
          ymean=reform((*self.blocks).mean[yidx]*1d6,self.params.nAlts)
          ysd=reform((*self.blocks).sd[yidx]*1d6,self.params.nAlts)
          ybw=reform((*self.blocks).bw[*,yidx],5,self.params.nAlts)
          ytitle=names[yidx]+' (ppm)'
          yrange=[0,100]
        END
        'rat':BEGIN
          ymean=reform(((*self.blocks).meanRatio[yidx]-1d0)*1d3,self.params.nAlts)
          ysd=reform((*self.blocks).sdRatio[yidx]*1d3,self.params.nAlts)
          ybw=reform(((*self.blocks).bwRatio[*,yidx]-1d0)*1d3,5,self.params.nAlts)
          ytitle=names[yidx+1]+'/'+names[0]+' ($10^{-3}$)'
          yrange=[-2000,1000]
        END
        'temp':BEGIN
          ymean=reform((*self.blocks).meant,self.params.nAlts)
          ysd=reform((*self.blocks).sdt,self.params.nAlts)
          ybw=reform((*self.blocks).bwt,5,self.params.nAlts)
          ytitle='Temp. (K)'
          yrange=[180,270]
        END
        ELSE: MESSAGE,'Unrecognised yname: '+STRTRIM(yname,2)+ $
                      '. Expected alt, pt, temp, mol<n> or delta<n>.'
      ENDCASE

      CASE xkind OF
        'rat':BEGIN
          xmean=reform(((*self.blocks).meanRatio[xidx]-1d0)*1d3,self.params.nAlts)
          xsd=reform((*self.blocks).sdRatio[xidx]*1d3,self.params.nAlts)
          xbw=reform(((*self.blocks).bwRatio[*,xidx]-1d0)*1d3,5,self.params.nAlts)
          ;----- Keep delta notation for the H2O/HDO pair, otherwise name the
          ; two molecules involved. -------;
          IF names[0] EQ 'H2O' AND names[xidx+1] EQ 'HDO' THEN BEGIN
            xtitle='$\delta D (‰)$'
          ENDIF ELSE BEGIN
            xtitle=names[xidx+1]+'/'+names[0]+' ($10^{-3}$)'
          ENDELSE
          xrange=[-2000,1000]
        END
        'mol':BEGIN
          xmean=reform((*self.blocks).mean[xidx]*1d6,self.params.nAlts)
          xsd=reform((*self.blocks).sd[xidx]*1d6,self.params.nAlts)
          xbw=reform((*self.blocks).bw[*,xidx],5,self.params.nAlts)
          xtitle=names[xidx]+' (ppm)'
          xrange=[0,500]
        END
        'temp':BEGIN
          xmean=reform((*self.blocks).meant,self.params.nAlts)
          xsd=reform((*self.blocks).sdt,self.params.nAlts)
          xbw=reform((*self.blocks).bwt,5,self.params.nAlts)
          xtitle='Temp. (K)'
          xrange=[180,270]
        END
        'pres':BEGIN
          xmean=reform((*self.blocks).meanp,self.params.nAlts)
          xsd=reform((*self.blocks).sdp,self.params.nAlts)
          xbw=reform((*self.blocks).bwp,5,self.params.nAlts)
          xtitle='Pres. (hPa)'
          xrange=[0,1]
        END
        ELSE: MESSAGE,'Unrecognised xname: '+STRTRIM(xname,2)+ $
                      '. Expected temp, pres, mol<n> or delta<n>.'
      ENDCASE

      if keyword_set(savPlot) then spag=boxplot(ymean,xbw,xtitle=xtitle,ytitle=ytitle,font_size=16,xthick=2,ythick=2,/horizontal,xrange=xrange,yrange=yrange,dimensions=[640,720],title=title,/buffer) else spag=boxplot(ymean,xbw,xtitle=xtitle,ytitle=ytitle,font_size=16,xthick=2,ythick=2,/horizontal,xrange=xrange,yrange=yrange,dimensions=[640,720],title=title)
      map=map('Geographic',center_longitude=90,limit=[-75,-180,75,180],position=[0.63,0.82,.93,1.06],font_size=18,/current,axis_style=0,linestyle=6,label_show=0)
      mc=mapcontinents(/continents,limit=[-75,-180,75,180])
      coords=[[(*self.blocks)[0].lonlims[0],(*self.blocks)[0].latlims[1]],[(*self.blocks)[0].lonlims[1],(*self.blocks)[0].latlims[1]],[(*self.blocks)[0].lonlims[1],(*self.blocks)[0].latlims[0]],[(*self.blocks)[0].lonlims[0],(*self.blocks)[0].latlims[0]]]
      region=polygon(coords,target=map,/data,fill_color='green')
      spag.select
      IF NOT KEYWORD_SET(nomean) THEN BEGIN
        meanRat=plot(xmean,ymean,thick=2,'r',/overplot)
        sdRat=polygon([xmean-xsd,reverse(xmean+xsd)],[ymean,reverse(ymean)],/data,/fill_background, FILL_COLOR="light steel blue",PATTERN_ORIENTATION=45,pattern_spacing=4,/overplot)
      ENDIF
      if keyword_set(savPlot) then begin
        spag.save,savPlot+mth_str+'_Rat '+string(mean(mols.lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(mols.lon[jrat[insrat]]),format='(F7.2)')+'.png',resolution=300
      endif


  if keyword_set(spaghetti) then begin
    IF KEYWORD_SET(outcols) THEN outcol=LIST()
    maxPts=MAX((*self.blocks).npts,mx)
    norbs=(*self.blocks)[mx].npts
    dex=(*(*self.blocks)[mx].orbids)
    orbs=self.getOcc(dex)
    ;print,dex
    for is=0,norbs-1 do begin
      
      IF KEYWORD_SET(colorby) THEN BEGIN
        CASE colorby OF
          'time': BEGIN
            color=[FIX(1.*is/norbs*255S),0S,0S]
          END
          'latitude': BEGIN
            color=[(orbs.lat[is]-self.params.latlims[0])/(self.params.latlims[1]-self.params.latlims[0])*255S,0S,0S]
          END
            
        ENDCASE
        ct=BYTARR(256,3)
        ct[*,0]=BINDGEN(256)
      ENDIF ELSE BEGIN
        color=[0S,0S,0S]
      ENDELSE
      ;print,color
      jk=where((*self.flags).aflag[*,dex[is]] eq 1,njk)
      IF njk EQ 0 THEN CONTINUE
      jk1=where(ymean ge orbs.alt[jk[0]])
      
      ;----- Select this occultation's profile explicitly. This replaces an
      ; EXECUTE that assembled the plot call from cmdx/cmdy strings. Writing
      ; the subscripts out showed that the pt case was reading orbs.pt[jk],
      ; which indexes the flattened 2-D array and so always returned the first
      ; occultation regardless of is. -------;
      CASE ykind OF
        'alt':  ydat=orbs.alt[jk]
        'pt':   ydat=orbs.pt[jk,is]
        'mol':  ydat=orbs.mixrat[jk,is,yidx]*1d6
        'rat':  ydat=(orbs.ratio[jk,is,yidx]-1d0)*1d3
        'temp': ydat=orbs.temp[jk,is]
      ENDCASE

      CASE xkind OF
        'rat':  xdat=(orbs.ratio[jk,is,xidx]-1d0)*1d3
        'mol':  xdat=orbs.mixrat[jk,is,xidx]*1d6
        'temp': xdat=orbs.temp[jk,is]
        'pres': xdat=orbs.pres[jk,is]
      ENDCASE

      p=plot(xdat,ydat,transparency=85,color=color,/overplot)
      
      IF KEYWORD_SET(outcols) THEN BEGIN
        onecall=DBLARR(3,N_ELEMENTS(jk))
        onecall[0,*]=orbs.alt[jk]
        onecall[1,*]=orbs.mixrat[jk,is,0]
        onecall[2,*]=(orbs.ratio[jk,is,0]-1d0)*1d3
        ;stop
        outcol.add,onecall
      ENDIF

    endfor
    spag=boxplot(ymean,xbw,xtitle=xtitle,ytitle=ytitle,font_size=16,thick=2,xthick=2,ythick=2,/horizontal,xrange=xrange,yrange=yrange,dimensions=[640,720],title=title,color='m',/OVERPLOT)
    IF KEYWORD_SET(outcols) THEN outcols=outcol.toarray(dimension=2)
    IF KEYWORD_SET(colorby) THEN cb=colorbar(target=spag,range=[self.params.latlims[0],self.params.latlims[1]],position=[0.03,0.0,0.35,0.03],title=colorby,font_size=10,textpos=1,rgb_table=ct)
  endif

end


pro ACEData::makeMeanSDTP,xname,yname,savPlot=savPlot

  COMPILE_OPT IDL2

  ;----- Unfinished: this body refers to n_longs, lat, jrat, latsran, temp,
  ; pres and alt, none of which are defined in this scope or passed in. It has
  ; never been able to run. Left in place for reference rather than deleted.
  ; -------;
  MESSAGE,'ACEData::makeMeanSDTP is not implemented.'

  ct=colortable(72,/reverse)

  temp_bins=dblarr(n_longs,n_lats)
  pres_bins=dblarr(n_longs,n_lats)
  nums=intarr(n_longs,n_lats)

  meanTempProf=dblarr(150,n_longs,n_lats)
  sdTempProf=dblarr(150,n_longs,n_lats)
  meanPresProf=dblarr(150,n_longs,n_lats)
  sdPresProf=dblarr(150,n_longs,n_lats)
  bpdTemp=dblarr(150,5,n_longs,n_lats)
  bpdPres=dblarr(150,5,n_longs,n_lats)

  for l=0,n_lats-1 do begin
    for m=0,n_longs-1 do begin
      insrat=where(lat[jrat] ge latsran[l] AND lat[jrat] le latsran[l+1] AND lon[jrat] ge longsran[m] AND lon[jrat] le longsran[m+1])
      if insrat[0] ne -1 then begin

        tempT=temp[*,insrat]
        tempP=pres[*,insrat]

        meanTempProf[0,m,l]=mean(tempT,dim=2)
        SDTempProf[0,m,l]=stddev(tempT,dim=2)
        meanPresProf[0,m,l]=mean(tempP,dim=2)
        SDPresProf[0,m,l]=stddev(tempP,dim=2)
        bpdTemp[0,0,m,l]=createboxplotdata(tempT)
        bpdPres[0,0,m,l]=createboxplotdata(tempP)

      endif
      ;stop
      if keyword_set(savPlot) then spagT=boxplot(alt,bpdTemp[*,*,m,l],xtitle='Temperature (C)',ytitle='Altitude (km)',font_size=16,xthick=2,ythick=2,/horizontal,xrange=[150,310],yrange=[0,40],dimensions=[640,720],title=string(mean(lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(lon[jrat[insrat]]),format='(F7.2)')+'                    ',/buffer) else spagT=boxplot(alt,bpdTemp[*,*,m,l],xtitle='Temperature (C)',ytitle='Altitude (km)',font_size=16,xthick=2,ythick=2,/horizontal,xrange=[150,310],yrange=[0,40],dimensions=[640,720],title=string(mean(lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(lon[jrat[insrat]]),format='(F7.2)')+'                    ',/buffer)
      map=map('Geographic',center_longitude=90,limit=[-75,-180,75,180],position=[0.63,0.82,.93,1.06],font_size=18,/current,axis_style=0,linestyle=6,label_show=0)
      mc=mapcontinents(/continents,limit=[-75,-180,75,180])
      coords=[[longsran[m],latsran[l+1]],[longsran[m+1],latsran[l+1]],[longsran[m+1],latsran[l]],[longsran[m],latsran[l]]]
      region=polygon(coords,target=map,/data,fill_color='green')
      spagT.select
      meanT=plot(meanTempProf[*,0,0],alt,'r',/overplot)
      sdT=polygon([meanTempProf[*,0,0]-sdTempProf[*,0,0],reverse(meanTempProf[*,0,0]+sdTempProf[*,0,0])],[alt,reverse(alt)],/data,/fill_background, FILL_COLOR="light steel blue",PATTERN_ORIENTATION=45,pattern_spacing=4,/overplot)
      ;stop


      if keyword_set(savPlot) then spagP=boxplot(alt,bpdPres[*,*,m,l]*1012.,xtitle='Pressure (hPa)',ytitle='Altitude (km)',font_size=16,xthick=2,ythick=2,/horizontal,xrange=[0,300],yrange=[0,30],dimensions=[640,720],title=string(mean(lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(lon[jrat[insrat]]),format='(F7.2)')+'                    ',/buffer) else spagP=boxplot(alt,bpdPres[*,*,m,l]*1012.,xtitle='Pressure (hPa)',ytitle='Altitude (km)',font_size=16,xthick=2,ythick=2,/horizontal,xrange=[0,300],yrange=[0,30],dimensions=[640,720],title=string(mean(lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(lon[jrat[insrat]]),format='(F7.2)')+'                    ')
      map=map('Geographic',center_longitude=90,limit=[-75,-180,75,180],position=[0.63,0.82,.93,1.06],font_size=18,/current,axis_style=0,linestyle=6,label_show=0)
      mc=mapcontinents(/continents,limit=[-75,-180,75,180])
      coords=[[longsran[m],latsran[l+1]],[longsran[m+1],latsran[l+1]],[longsran[m+1],latsran[l]],[longsran[m],latsran[l]]]
      region=polygon(coords,target=map,/data,fill_color='green')
      spagP.select
      meanP=plot(meanPresProf[*,0,0]*1012.,alt,'r',/overplot)
      sdP=polygon([meanPresProf[*,0,0]-sdPresProf[*,0,0],reverse(meanPresProf[*,0,0]+sdPresProf[*,0,0])]*1012.,[alt,reverse(alt)],/data,/fill_background, FILL_COLOR="light steel blue",PATTERN_ORIENTATION=45,pattern_spacing=4,/overplot)

      if keyword_set(savPlot) then begin
        spagT.save,savPlot+mth_str+'_Temp '+string(mean(lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(lon[jrat[insrat]]),format='(F7.2)')+'.png',resolution=300
        spagP.save,savPlot+mth_str+'_Pres '+string(mean(lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(lon[jrat[insrat]]),format='(F7.2)')+'.png',resolution=300
      endif
      spagT.close
      spagP.close

      ;      if keyword_set(spaghetti) then begin
      ;
      ;
      ;        for is=0,nspg-1d0 do begin
      ;          if keyword_set(raw) then gdsp=where(mixrat0[*,jrat[insrat[is]]] ne -999. and mixrat1[*,jrat[insrat[is]]] ne -999. and mixrat0_err[*,jrat[insrat[is]]] ne -999. and mixrat1_err[*,jrat[insrat[is]]] ne -999. and mixrat0_err[*,jrat[insrat[is]]] ne -888. and mixrat1_err[*,jrat[insrat[is]]] ne -888.)
      ;          if keyword_set(acecut) then gdsp=where(mixrat0[*,jrat[insrat[is]]] ne -999. and mixrat1[*,jrat[insrat[is]]] ne -999. and mixrat0_err[*,jrat[insrat[is]]] ne -999. and mixrat1_err[*,jrat[insrat[is]]] ne -999. and mixrat0_err[*,jrat[insrat[is]]] ne -888. and mixrat1_err[*,jrat[insrat[is]]] ne -888. and abs(mixrat0_err[*,jrat[insrat[is]]]) le 0.25 and abs(mixrat1_err[*,jrat[insrat[is]]]) le 0.25)
      ;          if keyword_set(fullcut) then gdsp=where(mixrat0[*,jrat[insrat[is]]] ne -999. and mixrat1[*,jrat[insrat[is]]] ne -999. and mixrat0_err[*,jrat[insrat[is]]] ne -999. and mixrat1_err[*,jrat[insrat[is]]] ne -999. and mixrat0_err[*,jrat[insrat[is]]] ne -888. and mixrat1_err[*,jrat[insrat[is]]] ne -888. and abs(mixrat0_err[*,jrat[insrat[is]]]) le 0.25 and abs(mixrat1_err[*,jrat[insrat[is]]]) le 0.25 and ratio[*,jrat[insrat[is]]] ge -1. and ratio[*,jrat[insrat[is]]] le 2.)
      ;          if gdsp[0] ne -1 then begin
      ;            spag.select
      ;            spag=plot(ratio[gdsp,jrat[insrat[gdsp]]]-1d0,alt[gdsp],'c',/overplot)
      ;            spag0.select
      ;            spag0=plot(mixrat0[gdsp,jrat[insrat[gdsp]]]*1d6,alt[gdsp],'c',/overplot)
      ;            spag1.select
      ;            spag1=plot(mixrat1[gdsp,jrat[insrat[gdsp]]]*1d6,alt[gdsp],'c',/overplot)
      ;
      ;          endif
      ;        endfor
      ;      endif
      ;      nums[m,l]=nspg
    endfor
  endfor



end

PRO ACEData::makeSpotMap,alts,altlims=altlims,outplot=outplot,postfix=postfix,ratIndex=ratIndex

  COMPILE_OPT IDL2

  IF N_ELEMENTS(ratIndex) EQ 0 THEN ratIndex=0

IF NOT KEYWORD_SET(altlims) THEN alims=DINDGEN(28)+12.5 ELSE alims=DINDGEN(altlims[1]-altlims[0]+1)+altlims[0]
IF NOT KEYWORD_SET(postfix) THEN postfix=''


nAlts=N_ELEMENTS(alts)

pnull=plot([0,1],[0,1],layout=[5,nalts/5,1],axis_style=0,dimensions=[1800,640],/nodata)

FOR i=0,nAlts-1 DO BEGIN
  
  map=map('Geographic',center_longitude=90,limit=[-90,-180,90,180],font_size=18,axis_style=0,linestyle=6,label_show=0,layout=[5,nAlts/5,i+1],margin=0.10,/current)
  mc=mapcontinents(/continents,limit=[-90,-180,90,180])
  xAlt=WHERE(alts[i] EQ alims,nxAlt)
  IF nxAlt EQ 0 THEN BEGIN
    PRINT,'Altitude ',alts[i],' not found in altlims. Skipping.'
    CONTINUE
  ENDIF
  xAlt=xAlt[0]

  norbs=(*self.blocks)[0,0,xAlt].npts
  dex=(*(*self.blocks)[0,0,xAlt].orbids)
  IF N_ELEMENTS(dex) EQ 0 THEN BEGIN
    PRINT,'No retrievals in the block at ',alts[i],' km. Skipping.'
    CONTINUE
  ENDIF
  orbs=self.getOcc(dex)

  kAlt=WHERE(alts[i] EQ orbs.alt,nkAlt)
  IF nkAlt EQ 0 THEN BEGIN
    PRINT,'Altitude ',alts[i],' not found on the retrieval grid. Skipping.'
    CONTINUE
  ENDIF
  ;stop
  dd=(orbs.ratio[kAlt,*,ratIndex]-1d0)*1d3
  p=plot(orbs.lon,orbs.lat,'.',sym_filled=1,sym_size=10,vert_colors=BYTSCL(dd,min=-760,max=-300),rgb_table=colortable(72,/reverse),font_size=18,title=STRING(alts[i],FORMAT='(F5.1)')+' km',/overplot)
ENDFOR

cb=colorbar(target=p,title='$\delta D (‰)$',range=[-.760d0,-.300d0]*1d3,position=[0.2,0.005,0.8,0.025],orientation=0,textpos=1,font_size=8)
outplot=map

END


PRO ACEData::stopView

  COMPILE_OPT IDL2

stop
  
END

PRO ACEData__define

  COMPILE_OPT IDL2

  A={paramsA,yearlims:[0.,0.],monthlims:[0.,0.],daylims:[0.,0.],hourlims:[0.,0.],jullims:[0.,0.],latlims:[0.,0.],lonlims:[0.,0.],altlims:[0.,0.],nTimes:0.,nLats:0.,nLons:0.,nAlts:0.,oParamName:'',oParamLims:[0.,0.],oParamNum:0.,mth_str:'',elim:0d0,cuttype:'',tMod:boolean(0),sMod:boolean(0),eMod:boolean(0),bMod:boolean(0),solidBlock:boolean(0),isDMP:boolean(0),dmp_path:''}
  struct={ACEData,raw:PTR_NEW(),hmz:PTR_NEW(),flags:PTR_NEW(),opMask:PTR_NEW(),blocks:PTR_NEW(),dmp:PTR_NEW(),params:A}
  
end