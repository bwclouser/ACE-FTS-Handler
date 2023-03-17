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
;       ACEData::isOpMaskDefined
;       ACEData::setOpMask
;       ACEData::areBlocksDefined
;       ACEData::resetBlocks
;       ACEData::reduceBlock
;       ACEData::calcBlocks
;       ACEData::makeContour
;       ACEData::makeMeridianPlot
;       ACEData::makeTimeSeries
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

  PTR_FREE,self.raw
  PTR_FREE,self.hmz
  IF PTR_VALID(self.flags) THEN PTR_FREE,self.flags
  IF PTR_VALID(self.opMask) THEN PTR_FREE,self.opMask
  IF PTR_VALID(self.blocks) THEN BEGIN    
    IF PTR_VALID((*self.blocks)[0].orbIDs) THEN PTR_FREE,(*self.blocks).orbIDs
    PTR_FREE,self.blocks
  ENDIF
  IF PTR_VALID(self.blocks) THEN PTR_FREE,self.blocks
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

PRO ACEData::loadFiles,fnames,version=version

  IF KEYWORD_SET(version) THEN BEGIN
    CASE version OF
      'p1p2': molnameid=11
      ELSE: molnameid=10
    ENDCASE
  ENDIF ELSE BEGIN
    molnameid=10
  ENDELSE

  n_files=size(fnames,/n_elements)
  mollist=list()
  
  for i=0,n_files-1 do begin

    ;get file id
    id='id'+string(i,format='(I02)')
    cmd='id=ncdf_open(fnames[i])'
    res=execute(cmd)

    ;get molecule name
    cmd0='mol=ncdf_varinq(id,molnameid)'
    res=execute(cmd0)

    ;get altitude
    cmd0='altid=ncdf_varid(id,"altitude")'
    cmd1='ncdf_varget,id,altid,alt'
    res=execute(cmd0+' & '+cmd1)

    ;get sunrise/sunset
    cmd0='srssid=ncdf_varid(id,"sunset_sunrise")'
    cmd1='ncdf_varget,id,srssid,srss'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get year
    cmd0='yearid=ncdf_varid(id,"year")'
    cmd1='ncdf_varget,id,yearid,year'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get month
    cmd0='monid=ncdf_varid(id,"month")'
    cmd1='ncdf_varget,id,monid,month'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get day
    cmd0='dayid=ncdf_varid(id,"day")'
    cmd1='ncdf_varget,id,dayid,day'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get hour
    cmd0='hourid=ncdf_varid(id,"hour")'
    cmd1='ncdf_varget,id,hourid,hour'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get latitude
    cmd0='latid=ncdf_varid(id,"latitude")'
    cmd1='ncdf_varget,id,latid,lat'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get longitude
    cmd0='lonid=ncdf_varid(id,"longitude")'
    cmd1='ncdf_varget,id,lonid,lon'+mol.name
    res=execute(cmd0+' & '+cmd1)

    cmd1='ncdf_varget,id,molnameid,'+mol.name
    res=execute(cmd1)

    ;get orbit id
    cmd0='orbid=ncdf_varid(id,"orbit")'
    cmd1='ncdf_varget,id,orbid,orbit_'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get temperature
    cmd0='tempid=ncdf_varid(id,"temperature")'
    cmd1='ncdf_varget,id,tempid,temp_'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get temperature fits
    cmd0='tempfitid=ncdf_varid(id,"temperature_fit")'
    cmd1='ncdf_varget,id,tempfitid,tempfit_'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get pressure
    cmd0='presid=ncdf_varid(id,"pressure")'
    cmd1='ncdf_varget,id,presid,pres_'+mol.name
    res=execute(cmd0+' & '+cmd1)

    ;get error for molecule
    cmd0='mol_err=ncdf_varinq(id,molnameid+1)'
    res=execute(cmd0)

    cmd1='ncdf_varget,id,molnameid+1,'+mol_err.name
    res=execute(cmd1)

    cmd0='npass_'+mol.name+'=n_elements(lon'+mol.name+')'
    cmd1=mol.name+'_str={name:"'+mol.name+'",npass:npass_'+mol.name+',alt:alt,srss:srss'+mol.name+',orbit:orbit_'+mol.name+',year:year'+mol.name+',month:month'+mol.name+',day:day'+mol.name+',hour:hour'+mol.name+',lat:lat'+mol.name+',lon:lon'+mol.name+',mixrat:'+mol.name+',mixrat_err:'+mol.name+'_error,temp:temp_'+mol.name+',pres:pres_'+mol.name+'}'
    cmd2='mollist.add,'+mol.name+'_str'
    res=execute(cmd0+' & '+cmd1+' & '+cmd2)
    
  endfor
  
  (*self.raw)=mollist

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

;get file id
id=ncdf_open(self.params.dmp_path)

;get sunrise/sunset
cmd0='srssid=ncdf_varid(id,"Sunset Sunrise")'
cmd1='ncdf_varget,id,srssid,srss'
res=execute(cmd0+' & '+cmd1)

;get orbit id
cmd0='orbid=ncdf_varid(id,"Orbit")'
cmd1='ncdf_varget,id,orbid,orbit'
res=execute(cmd0+' & '+cmd1)

;get pressure
cmd0='presid=ncdf_varid(id,"Pressure")'
cmd1='ncdf_varget,id,presid,pres'
res=execute(cmd0+' & '+cmd1)

;get zonal wind
cmd0='uwindid=ncdf_varid(id,"U")'
cmd1='ncdf_varget,id,uwindid,uwind
res=execute(cmd0+' & '+cmd1)

;get meridional wind
cmd0='vwindid=ncdf_varid(id,"V")'
cmd1='ncdf_varget,id,vwindid,vwind
res=execute(cmd0+' & '+cmd1)

;get temperature
cmd0='tempid=ncdf_varid(id,"Temperature")'
cmd1='ncdf_varget,id,tempid,temp
res=execute(cmd0+' & '+cmd1)

;get potential temperature
cmd0='thetaid=ncdf_varid(id,"Theta")'
cmd1='ncdf_varget,id,thetaid,theta
res=execute(cmd0+' & '+cmd1)

;get geopotential height
cmd0='gphid=ncdf_varid(id,"GPH")'
cmd1='ncdf_varget,id,gphid,gph
res=execute(cmd0+' & '+cmd1)

;get altitude
cmd0='altid=ncdf_varid(id,"Altitude")'
cmd1='ncdf_varget,id,altid,alt'
res=execute(cmd0+' & '+cmd1)

;get potential vorticity
cmd0='pvid=ncdf_varid(id,"PV")'
cmd1='ncdf_varget,id,pvid,pv'
res=execute(cmd0+' & '+cmd1)

;get scaled potential vorticity
cmd0='spvid=ncdf_varid(id,"sPV")'
cmd1='ncdf_varget,id,spvid,spv'
res=execute(cmd0+' & '+cmd1)

;get relative vorticity
cmd0='rvid=ncdf_varid(id,"RV")'
cmd1='ncdf_varget,id,rvid,rv'
res=execute(cmd0+' & '+cmd1)

;get lapse rate
cmd0='dtdzid=ncdf_varid(id,"dT_dZ")'
cmd1='ncdf_varget,id,dtdzid,dtdz'
res=execute(cmd0+' & '+cmd1)

;get static stability
cmd0='ssid=ncdf_varid(id,"Static_Stability")'
cmd1='ncdf_varget,id,ssid,ss'
res=execute(cmd0+' & '+cmd1)

;get equivalent latitude
cmd0='eqlid=ncdf_varid(id,"EqL")'
cmd1='ncdf_varget,id,eqlid,eql'
res=execute(cmd0+' & '+cmd1)

;get Normalized Horizontal (Isentropic) sPV Gradient
cmd0='hpvgid=ncdf_varid(id,"HorsPVgrad")'
cmd1='ncdf_varget,id,hpvgid,hpvg'
res=execute(cmd0+' & '+cmd1)

;get Normalizing factor for Horizontal sPV Gradient
cmd0='apvgid=ncdf_varid(id,"avgsPVgrad")'
cmd1='ncdf_varget,id,apvgid,apvg'
res=execute(cmd0+' & '+cmd1)

;get Horizontal (isobaric) Temperature Gradient
cmd0='tgid=ncdf_varid(id,"Tgrad")'
cmd1='ncdf_varget,id,tgid,tg'
res=execute(cmd0+' & '+cmd1)

;get montgomery stream function
cmd0='msfid=ncdf_varid(id,"MSF")'
cmd1='ncdf_varget,id,msfid,msf'
res=execute(cmd0+' & '+cmd1)

;get ozone
cmd0='o3id=ncdf_varid(id,"AssimO3")'
cmd1='ncdf_varget,id,o3id,o3'
res=execute(cmd0+' & '+cmd1)

n_orbs=n_elements(orbit)
n_alts=75.

str={srss:srss,orbit:orbit,pres:reform(pres,n_alts,n_orbs),U:reform(uwind,n_alts,n_orbs),V:reform(vwind,n_alts,n_orbs),temp:reform(temp,n_alts,n_orbs),theta:reform(theta,n_alts,n_orbs),GPH:reform(GPH,n_alts,n_orbs),alt:alt,pv:reform(pv,n_alts,n_orbs),spv:reform(spv,n_alts,n_orbs),rv:reform(rv,n_alts,n_orbs),dtdz:reform(dtdz,n_alts,n_orbs),ss:reform(ss,n_alts,n_orbs),eql:reform(eql,n_alts,n_orbs),hpvg:reform(hpvg,n_alts,n_orbs),apvg:reform(apvg,n_alts,n_orbs),tg:reform(tg,n_alts,n_orbs),msf:reform(msf,n_alts,n_orbs),o3:reform(o3,n_alts,n_orbs)}

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
;       finds retrievals common to ALL molecules and stores them in a new data structure.
;       (Currently, only 2 molecules are supported.) Additionally, the first molecule provided
;       to ACEData::loadFiles is considered to be the 'primary' molecule, and ACEData::harmonize
;       calculates the ratio of each subsequent molecule to the primary molecule.
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
;-

PRO ACEData::harmonize,elements,addDMP=addDMP

mol0=(*self.raw)[elements[0]]
name0=mol0.name

orb0=mol0.orbit
srss0=mol0.srss

n_mols=N_ELEMENTS(elements)

if n_elements(elements) ne 1 then begin

  mol1=(*self.raw)[elements[1]]
  name1=mol1.name

  mndim=min([mol0.npass,mol1.npass],jloc)         ;determines which structure has fewer elements

  orb1=mol1.orbit
  srss1=mol1.srss

  ind0=[]
  ind1=[]
  
  IF KEYWORD_SET(addDMP) THEN BEGIN
    IF PTR_VALID(self.dmp) THEN BEGIN
      IF *self.dmp ne !NULL THEN BEGIN
        isFind=1
        ind2=[]
      ENDIF ELSE BEGIN
        PRINT,'ERROR: DMP NOT LOADED'
      ENDELSE
    ENDIF ELSE BEGIN
      PRINT,'ERROR: POINTER TO DMP NOT VALID'
    ENDELSE
  ENDIF ELSE BEGIN
    isFind=0
  ENDELSE
  
  ;----- This if/else statement and for loop finds all orbits that are present in BOTH structures
  ; and returns their respective indices in the ind0 and ind1 arrays -------;
  if jloc[0] eq 0 then begin
    for i=0,mndim-1 do begin
      x=WHERE(orb0[i] EQ orb1 AND srss0[i] EQ srss1)
      if x[0] NE -1 then begin
        ind0=[ind0,i]
        ind1=[ind1,x]
        IF isFind THEN BEGIN
          y=WHERE(orb0[i] EQ (*self.dmp).orbit AND srss0[i] EQ (*self.dmp).srss)
          ind2=[ind2,y]
        ENDIF
      endif
    endfor
  endif else begin
    for i=0,mndim-1 do begin
      x=where(orb1[i] eq orb0 and srss1[i] eq srss0)
      if x[0] ne -1 then begin
        ind0=[ind0,x]
        ind1=[ind1,i]
        IF isFind THEN BEGIN
          y=WHERE(orb1[i] EQ (*self.dmp).orbit AND srss1[i] EQ (*self.dmp).srss)
          ind2=[ind2,y]
        ENDIF
      endif
    endfor
  endelse

  ;-------- This block defines the structure in the ratio and no ratio cases --------;

;  if KEYWORD_SET(addDMP) THEN BEGIN
;    orbits=mol0.orbit[ind0]
;    srsss=mol0.srss[ind0]
;    npass=n_elements(ind0)
;    if PTR_VALID(self.dmp) THEN BEGIN
;      IF *self.dmp ne !NULL THEN BEGIN
;        n_orbs=n_elements((*self.dmp).orbit)
;        n_alts=75
;        pt=fltarr(150,npass)
;        FOR j=0,npass-1 DO BEGIN
;          k=WHERE((*self.dmp).orbit EQ orbits[j] AND (*self.dmp).srss EQ srsss[j])
;          pt[0:n_alts-1,j]=(*self.dmp).theta[*,k]
;        ENDFOR
;        molout={nmols:n_mols,name0:mol0.name,name1:mol1.name,npass:n_elements(ind0),alt:mol0.alt,srss:mol0.srss[ind0],orbit:mol0.orbit[ind0],year:mol0.year[ind0],month:mol0.month[ind0],day:mol0.day[ind0],hour:mol0.hour[ind0],lat:mol0.lat[ind0],lon:mol0.lon[ind0],mixrat0:mol0.mixrat[*,ind0],mixrat0_err:mol0.mixrat_err[*,ind0],mixrat1:mol1.mixrat[*,ind1],ratio:mol1.mixrat[*,ind1]/mol0.mixrat[*,ind0],mixrat1_err:mol1.mixrat_err[*,ind1],temp:mol0.temp[*,ind0],pres:mol0.pres[*,ind0],pt:pt}
;      ENDIF ELSE BEGIN
;        PRINT,'ERROR: DMP NOT LOADED'
;      ENDELSE
;    ENDIF ELSE BEGIN
;      PRINT,'ERROR: POINTER TO DMP NOT VALID'
;    ENDELSE
;  ENDIF ELSE BEGIN
;    molout={nmols:n_mols,name0:mol0.name,name1:mol1.name,npass:n_elements(ind0),alt:mol0.alt,srss:mol0.srss[ind0],orbit:mol0.orbit[ind0],year:mol0.year[ind0],month:mol0.month[ind0],day:mol0.day[ind0],hour:mol0.hour[ind0],lat:mol0.lat[ind0],lon:mol0.lon[ind0],mixrat0:mol0.mixrat[*,ind0],mixrat0_err:mol0.mixrat_err[*,ind0],mixrat1:mol1.mixrat[*,ind1],ratio:mol1.mixrat[*,ind1]/mol0.mixrat[*,ind0],mixrat1_err:mol1.mixrat_err[*,ind1],temp:mol0.temp[*,ind0],pres:mol0.pres[*,ind0],pt:mol0.temp*(1d0/mol0.pres)^.286d0}
;  ENDELSE  

  juldays=julday(mol0.month[ind0],mol0.day[ind0],mol0.year[ind0],mol0.hour[ind0],0d0,0d0)

  IF isFind THEN BEGIN
    pt=(*self.dmp).theta[*,ind2]
    pt=[pt,FLTARR(75,N_ELEMENTS(ind2))]
    molout={nmols:n_mols,name0:mol0.name,name1:mol1.name,npass:n_elements(ind0),alt:mol0.alt,srss:mol0.srss[ind0],orbit:mol0.orbit[ind0],year:mol0.year[ind0],month:mol0.month[ind0],day:mol0.day[ind0],hour:mol0.hour[ind0],juldays:juldays,lat:mol0.lat[ind0],lon:mol0.lon[ind0],mixrat0:mol0.mixrat[*,ind0],mixrat0_err:mol0.mixrat_err[*,ind0],mixrat1:mol1.mixrat[*,ind1],ratio:mol1.mixrat[*,ind1]/mol0.mixrat[*,ind0],mixrat1_err:mol1.mixrat_err[*,ind1],temp:mol0.temp[*,ind0],pres:mol0.pres[*,ind0],pt:pt}
  ENDIF ELSE BEGIN
    molout={nmols:n_mols,name0:mol0.name,name1:mol1.name,npass:n_elements(ind0),alt:mol0.alt,srss:mol0.srss[ind0],orbit:mol0.orbit[ind0],year:mol0.year[ind0],month:mol0.month[ind0],day:mol0.day[ind0],hour:mol0.hour[ind0],juldays:juldays,lat:mol0.lat[ind0],lon:mol0.lon[ind0],mixrat0:mol0.mixrat[*,ind0],mixrat0_err:mol0.mixrat_err[*,ind0],mixrat1:mol1.mixrat[*,ind1],ratio:mol1.mixrat[*,ind1]/mol0.mixrat[*,ind0],mixrat1_err:mol1.mixrat_err[*,ind1],temp:mol0.temp[*,ind0],pres:mol0.pres[*,ind0],pt:mol0.temp*(1d0/mol0.pres)^.286d0} 
  ENDELSE
  
endif else begin

  juldays=julday(mol0.month,mol0.day,mol0.year,mol0.hour,0d0,0d0)

  molout={nmols:n_mols,name0:mol0.name,npass:mol0.npass,alt:mol0.alt,srss:mol0.srss,orbit:mol0.orbit,year:mol0.year,month:mol0.month,day:mol0.day,hour:mol0.hour,juldays:juldays,lat:mol0.lat,lon:mol0.lon,mixrat0:mol0.mixrat,mixrat0_err:mol0.mixrat_err,temp:mol0.temp,pres:mol0.pres,pt:mol0.temp*(1d0/mol0.pres)^.286d0}

endelse

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

  names=tag_names(self.params)
  isParam=where(names.contains(pName,/FOLD_CASE) eq 1)
  if isParam[0] ne -1 then begin
    cmd='self.params.'+pName+'=pVal'
    res=execute(cmd)
    result=1
    if pName EQ 'latlims' OR pName EQ 'lonlims' OR pName EQ 'altlims' OR pName EQ 'nLats' OR pName EQ 'nLons' OR pName eq 'nAlts' THEN self.params.sMod=1
    if pName EQ 'yearlims' OR pName EQ 'monthlims' OR pName eq 'daylims' OR pName EQ 'hourlims' OR pName eq 'nTimes' THEN self.params.tMod=1
    if pName EQ 'elims' OR pName EQ 'cuttype' THEN self.params.eMod=1
    self.params.bMod=1
  endif else begin
    result=0
    CASE pName OF
      'ctStandard': BEGIN
        self.params.latlims=[-70,70]
        self.params.lonlims=[-180,180]
        self.params.altlims=[16,17]
        self.params.nLats=12
        self.params.nLons=16
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
          self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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
        self.params.dmp_path='/media/benjamin/Data/ACE/ACE-DMPSv4p1p2/ACEFTS_L2_v4p1_GLC_GEOS5MERRA2_DynEqL_jv301.nc4'
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

  names=tag_names(self.params)
  isParam=where(names.contains(pName,/FOLD_CASE) eq 1)
  if isParam[0] ne -1 then begin
    cmd='result=self.params.'+pName
    res=execute(cmd)
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

FUNCTION ACEData::getOcc,orbids

  nOccs=n_elements(orbids)
  inds=orbids

;  inds=ulonarr(nOccs)
;
;  for i=0,nOccs-1 do begin
;    inds[i]=where((*self.hmz).orbit eq orbids[i])
;  endfor
;
;  inds=inds[where(inds ne -1,compliment=k)]
;  nOccs=n_elements(inds)
;
;  if k[0] ne -1 then print,'Orbits ',orbids[k],' not found.'

  CASE (*self.hmz).nmols OF
    1: BEGIN

      out={nmols:(*self.hmz).nmols,name0:(*self.hmz).name0,npass:nOccs,alt:findgen(150)+.5,srss:(*self.hmz).srss[inds],orbit:(*self.hmz).orbit[inds],year:(*self.hmz).year[inds],month:(*self.hmz).month[inds],day:(*self.hmz).day[inds],hour:(*self.hmz).hour[inds],lat:(*self.hmz).lat[inds],lon:(*self.hmz).lon[inds],mixrat0:(*self.hmz).mixrat0[*,inds],mixrat0_err:(*self.hmz).mixrat0_err[*,inds],temp:(*self.hmz).temp[*,inds],pres:(*self.hmz).pres[*,inds],pt:(*self.hmz).pt[*,inds]}

    END
    2: BEGIN

      out={nmols:(*self.hmz).nmols,name0:(*self.hmz).name0,name1:(*self.hmz).name1,npass:nOccs,alt:findgen(150)+.5,srss:(*self.hmz).srss[inds],orbit:(*self.hmz).orbit[inds],year:(*self.hmz).year[inds],month:(*self.hmz).month[inds],day:(*self.hmz).day[inds],hour:(*self.hmz).hour[inds],lat:(*self.hmz).lat[inds],lon:(*self.hmz).lon[inds],mixrat0:(*self.hmz).mixrat0[*,inds],mixrat0_err:(*self.hmz).mixrat0_err[*,inds],mixrat1:(*self.hmz).mixrat1[*,inds],mixrat1_err:(*self.hmz).mixrat1_err[*,inds],ratio:(*self.hmz).ratio[*,inds],temp:(*self.hmz).temp[*,inds],pres:(*self.hmz).pres[*,inds],pt:(*self.hmz).pt[*,inds]}

    END
  ENDCASE

  return,out

END

FUNCTION ACEData::getRawData

  IF PTR_VALID(self.raw) THEN BEGIN
    RETURN,(*self.raw)
  ENDIF ELSE BEGIN
    PRINT,'Raw data not defined.'
    RETURN,!NULL
  ENDELSE

END

FUNCTION ACEData::getHmzData

  IF PTR_VALID(self.hmz) THEN BEGIN
    RETURN,(*self.hmz)
  ENDIF ELSE BEGIN
    PRINT,'Harmonized data not defined.'
    RETURN,!NULL
  ENDELSE


END

FUNCTION ACEData::getBlocks

  IF PTR_VALID(self.blocks) THEN BEGIN
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

  npass=(*self.hmz).npass

  if self.flags eq !null then begin

    self.flags=PTR_NEW(/ALLOCATE)
    B={flags,tflag:boolarr(npass),sflag:boolarr(npass),eflag:boolarr(npass),aflag:boolarr(150,npass)}
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
 
  mth_lab=['J','F','M','A','M','J','J','A','S','O','N','D']
  mth_str=''

  self.areFlagsDefined
  
  IF self.params.solidBlock EQ 0 THEN BEGIN
    
    yearInds=indgen(self.params.yearlims[1]-self.params.yearlims[0]+1)+self.params.yearlims[0]
    yearIncl=intarr((*self.hmz).npass)
    for ii=0,n_elements(yearInds)-1 do begin
      jyr=where((*self.hmz).year eq yearInds[ii])
      yearIncl[jyr]=1
    endfor
    
    if self.params.monthlims[1] ne 0 then begin 
      if self.params.monthlims[0] ge self.params.monthlims[1] then monthInds=[indgen(12.-self.params.monthlims[0]+1)+self.params.monthlims[0],indgen(self.params.monthlims[1])+1] else monthInds=indgen(self.params.monthlims[1]-self.params.monthlims[0]+1)+self.params.monthlims[0]
      monthIncl=intarr((*self.hmz).npass)
      for ii=0,n_elements(monthInds)-1 do begin
        jmth=where((*self.hmz).month eq monthInds[ii])
        monthIncl[jmth]=1
        mth_str+=mth_lab[monthInds[ii]-1]
      endfor
    endif else begin
      monthIncl=intarr((*self.hmz).npass)+1
    endelse
    
    if self.params.daylims[1] ne 0 then begin 
      if self.params.daylims[0] ge self.params.daylims[1] then begin
        ;add code here that can handle rolling over months... maybe this will have to be lumped together with monthlims?
      endif else begin
        dayInds=indgen(self.params.daylims[1]-self.params.daylims[0]+1)+self.params.daylims[0]
      endelse
      dayIncl=intarr((*self.hmz).npass)
      for ii=0,n_elements(dayInds)-1 do begin
        jdy=where((*self.hmz).day eq dayInds[ii])
        dayIncl[jdy]=1
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
    
    ;stop
  
    self.params.mth_str=mth_str
  
  ENDIF ELSE BEGIN
    
    dayLast=julday(self.params.monthlims[1],self.params.daylims[1],self.params.yearlims[1],self.params.hourlims[1],0,0)
    dayFirst=julday(self.params.monthlims[0],self.params.daylims[0],self.params.yearlims[0],self.params.hourlims[0],0,0)

    julIncl=WHERE((*self.hmz).juldays GE self.params.jullims[0] AND (*self.hmz).juldays LE self.params.jullims[1])
    (*self.flags).tflag=julIncl
    
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

  self.areFlagsDefined

  kran=WHERE((*self.hmz).lat GE self.params.latlims[0] AND (*self.hmz).lat LE self.params.latlims[1] AND (*self.hmz).lon ge self.params.lonlims[0] AND (*self.hmz).lon LE self.params.lonlims[1],complement=jran)
  (*self.flags).sflag[kran]=1
  (*self.flags).sflag[jran]=0
  
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
  
  (*self.flags).eflag=1
  IF KEYWORD_SET(cutByTemp) THEN BEGIN
    k=WHERE((*self.hmz).temp LE cutByTemp[0])
    z=k MOD 150
    v=WHERE(z GE cutByTemp[1] AND z LE cutByTemp[2])
    
    (*self.flags).eflag[(floor(k[v]/150))[uniq(floor(k[v]/150))]]=0
    ;stop
  ENDIF
  
  (*self.flags).aflag[*,*]=1
  npass=(*self.hmz).npass
  
  CASE (*self.hmz).nmols OF
    1: BEGIN

      bad1=where((*self.hmz).mixrat0 eq -999. or (*self.hmz).mixrat0_err eq -999. or (*self.hmz).mixrat0_err eq -888.)
      IF self.params.cuttype.contains('pointwise',/FOLD_CASE) THEN BEGIN
        bad2=where((*self.hmz).mixrat0_err ge self.params.eLim)
      ENDIF
      IF self.params.cuttype.contains('scanwise',/FOLD_CASE) THEN BEGIN
        (*self.flags).eflag[FLOOR(WHERE((*self.hmz).mixrat0_err ge self.params.eLim)/150.)]=0
      ENDIF
    END
    2: BEGIN

      bad1=where((*self.hmz).mixrat0 eq -999. or (*self.hmz).mixrat0_err eq -999. or (*self.hmz).mixrat0_err eq -888. or (*self.hmz).mixrat1 eq -999. or (*self.hmz).mixrat1_err eq -999. or (*self.hmz).mixrat1_err eq -888.)
      IF self.params.cuttype.contains('pointwise',/FOLD_CASE) THEN BEGIN
        bad2=where((*self.hmz).mixrat0_err ge self.params.eLim OR (*self.hmz).mixrat1_err ge self.params.eLim)
      ENDIF
      IF self.params.cuttype.contains('scanwise',/FOLD_CASE) THEN BEGIN
        (*self.flags).eflag[FLOOR(WHERE((*self.hmz).mixrat0_err GE self.params.eLim AND (*self.hmz).mixrat1_err GE self.params.elim)/150.)]=0
      ENDIF
    END
  ENDCASE

  (*self.flags).aflag[bad1]=0
  IF self.params.cuttype.contains('pointwise',/FOLD_CASE) THEN BEGIN
    (*self.flags).aflag[bad2]=0
  ENDIF
  ;stop

  
  self.params.eMod=0

END

PRO ACEData::errorByExtrema



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

  if self.opMask eq !null then begin
    self.opMask=PTR_NEW(/ALLOCATE)
    C=boolarr(150,(*self.hmz).npass)
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
    k=WHERE((*self.hmz).juldays GE jdays[0] AND (*self.hmz).juldays LE jdays[1] AND ((*self.hmz).lat GE blats[0]) AND ((*self.hmz).lat LE blats[1]) AND ((*self.hmz).lon GE blons[0]) AND ((*self.hmz).lon LE blons[1]),complement=j)
  ENDIF ELSE BEGIN
    k=WHERE(((*self.hmz).year GE byears[0]) AND ((*self.hmz).year LE byears[1]) AND ((*self.hmz).month GE bmonths[0]) AND ((*self.hmz).month LE bmonths[1]) AND ((*self.hmz).day GE bdays[0]) AND ((*self.hmz).day LE bdays[1]) AND ((*self.hmz).hour GE bhours[0]) AND ((*self.hmz).hour LE bhours[1]) AND ((*self.hmz).lat GE blats[0]) AND ((*self.hmz).lat LE blats[1]) AND ((*self.hmz).lon GE blons[0]) AND ((*self.hmz).lon LE blons[1]),complement=j)
  ENDELSE
  (*self.opMask)[*,k]=1
  (*self.opMask)[*,j]=0
  ;deal with altitudes here

  if balts[0] gt 0 or balts[1] lt 150 then begin
    
    altind0=intarr(150)
    altind0[*]=1
    
    if balts[0] gt 0 then altind0[WHERE((*self.hmz).alt le balts[0])]=0
    if balts[1] lt 150 then altind0[WHERE((*self.hmz).alt ge balts[1])]=0
    altind=rebin(altind0,150,(*self.hmz).npass)
    
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

  if self.blocks eq !null then begin
    self.blocks=PTR_NEW(/ALLOCATE)
    if (*self.hmz).nmols eq 1 then begin
      D={nmols:(*self.hmz).nmols,name0:(*self.hmz).name0,yearlims:[0,0],monthlims:[0,0],daylims:[0,0],hourlims:[0,0],latlims:[0,0],lonlims:[0,0],altlims:[0,0],mean0:0.,sd0:0.,bw0:fltarr(5),meant:0.,sdt:0.,bwt:fltarr(5),meanp:0.,sdp:0.,bwp:fltarr(5),meanLat:0.,sdLat:0.,meanLon:0.,sdLon:0.,nPts:0,orbIDs:PTR_NEW(/ALLOCATE)}
    endif else begin
      D={nmols:(*self.hmz).nmols,name0:(*self.hmz).name0,yearlims:[0,0],monthlims:[0,0],daylims:[0,0],hourlims:[0,0],latlims:[0,0],lonlims:[0,0],altlims:[0,0],mean0:0.,sd0:0.,bw0:fltarr(5),mean1:0.,sd1:0.,bw1:fltarr(5),meanRat:0.,sdRat:0.,bwRat:fltarr(5),meant:0.,sdt:0.,bwt:fltarr(5),meanp:0.,sdp:0.,bwp:fltarr(5),meanLat:0.,sdLat:0.,meanLon:0.,sdLon:0.,nPts:0,orbIDs:PTR_NEW(/ALLOCATE)}
    endelse
      (*self.blocks)=replicate(D,self.params.nLats,self.params.nLons,self.params.nAlts,self.params.nTimes)
      self.params.bMod=0
  endif

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

  self.areBlocksDefined

  if self.blocks ne !null then begin
    PTR_FREE,(*self.blocks).orbIDs
    PTR_FREE,self.blocks
    self.blocks=PTR_NEW(/ALLOCATE)
    if (*self.hmz).nmols eq 1 then begin
      D={nmols:(*self.hmz).nmols,name0:(*self.hmz).name0,yearlims:[0d0,0d0],monthlims:[0d0,0d0],daylims:[0d0,0d0],hourlims:[0d0,0d0],jullims:[0d0,0d0],latlims:[0d0,0d0],lonlims:[0d0,0d0],altlims:[0d0,0d0],mean0:0.,sd0:0.,bw0:fltarr(5),meant:0.,sdt:0.,bwt:fltarr(5),meanp:0.,sdp:0.,bwp:fltarr(5),meanpt:0.,sdpt:0.,bwpt:fltarr(5),meanLat:0.,sdLat:0.,meanLon:0.,sdLon:0.,nPts:0,orbIDs:PTR_NEW()}
    endif else begin
      D={nmols:(*self.hmz).nmols,name0:(*self.hmz).name0,yearlims:[0d0,0d0],monthlims:[0d0,0d0],daylims:[0d0,0d0],hourlims:[0d0,0d0],jullims:[0d0,0d0],latlims:[0d0,0d0],lonlims:[0d0,0d0],altlims:[0d0,0d0],mean0:0.,sd0:0.,bw0:fltarr(5),mean1:0.,sd1:0.,bw1:fltarr(5),meanRat:0.,sdRat:0.,bwRat:fltarr(5),meant:0.,sdt:0.,bwt:fltarr(5),meanp:0.,sdp:0.,bwp:fltarr(5),meanpt:0.,sdpt:0.,bwpt:fltarr(5),meanLat:0.,sdLat:0.,meanLon:0.,sdLon:0.,nPts:0,orbIDs:PTR_NEW()}
    endelse
    (*self.blocks)=replicate(D,self.params.nLats,self.params.nLons,self.params.nAlts,self.params.nTimes)
    nblocks=n_elements(*self.blocks)
    for i=0,nblocks-1 do (*self.blocks)[i].orbIDs=PTR_NEW(/ALLOCATE)
    self.params.bMod=0
  endif

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

PRO ACEData::reduceBlock,i,j,k,l

  npass=(*self.hmz).npass

  goodOrbs=((*self.flags).sflag AND (*self.flags).tflag AND (*self.flags).eflag)
  startOrbs=where(goodOrbs eq 1)
  
  tempAFlags=(*self.flags).aflag[*,startOrbs]
  tempOpMask=(*self.opMask)[*,startOrbs]
  
  fFlag=((*self.flags).aflag[*,startOrbs] AND (*self.opMask)[*,startOrbs])
  
  kGood=where(fFlag eq 1)
  finOrbs=floor(kGood/150)
  ndat=n_elements(kGood)
  
  mr0=(*self.hmz).mixrat0[*,startOrbs]
  mr0g=mr0[kGood]
  temp=(*self.hmz).temp[*,startOrbs]
  tempg=temp[kGood]
  pres=(*self.hmz).pres[*,startOrbs]
  presg=pres[kGood]
  pt=(*self.hmz).pt[*,startOrbs]
  ptg=pt[kGood]
  (*self.blocks)[i,j,k,l].mean0=mean(mr0g)
  (*self.blocks)[i,j,k,l].sd0=stddev(mr0g)
  (*self.blocks)[i,j,k,l].meant=mean(tempg)
  (*self.blocks)[i,j,k,l].sdt=stddev(tempg)
  (*self.blocks)[i,j,k,l].meanp=mean(presg)
  (*self.blocks)[i,j,k,l].sdp=stddev(presg)
  (*self.blocks)[i,j,k,l].meanpt=mean(ptg)
  (*self.blocks)[i,j,k,l].sdpt=stddev(ptg)
  (*self.blocks)[i,j,k,l].meanLat=mean((*self.hmz).lat[startOrbs[finOrbs]])
  (*self.blocks)[i,j,k,l].meanLon=mean((*self.hmz).lon[startOrbs[finOrbs]])
  (*self.blocks)[i,j,k,l].nPts=ndat
  ;stop
  if ndat ge 5 then begin
    (*self.blocks)[i,j,k,l].bw0=createboxplotdata(mr0g)
    (*self.blocks)[i,j,k,l].bwt=createboxplotdata(tempg)
    (*self.blocks)[i,j,k,l].bwp=createboxplotdata(presg)
    (*self.blocks)[i,j,k,l].bwpt=createboxplotdata(ptg)
  endif
  if (*self.hmz).nmols gt 1 then begin
    mr1=(*self.hmz).mixrat0[*,startOrbs]
    mr1g=mr1[kGood]
    rat0=(*self.hmz).ratio[*,startOrbs]
    rat0g=rat0[kGood]
    (*self.blocks)[i,j,k,l].mean1=mean(mr1g)
    (*self.blocks)[i,j,k,l].sd1=stddev(mr1g)
    (*self.blocks)[i,j,k,l].meanRat=mean(rat0g)
    (*self.blocks)[i,j,k,l].sdRat=stddev(rat0g)

    if ndat ge 5 then begin
      (*self.blocks)[i,j,k,l].bw1=createboxplotdata(mr1g)
      (*self.blocks)[i,j,k,l].bwRat=createboxplotdata(rat0g)
    endif
  endif
  ;stop
  ;(*(*self.blocks)[i,j,k,l].orbIDs)=(*self.hmz).orbit[startOrbs[finOrbs]]
  (*(*self.blocks)[i,j,k,l].orbIDs)=startOrbs[finOrbs]
  
END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::calcBlock
;
; PURPOSE:
; 
;       The ACEData::calcBlock calculates the parameters of interest for each lat/lot/alt/time
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

  daysInMonth=[31,28,31,30,31,30,31,31,30,31,30,31]

  IF self.params.sMod EQ 1 OR self.params.tMod EQ 1 OR self.params.eMod EQ 1 OR self.params.bMod EQ 1 THEN BEGIN
    print,'Parameters have Changed, resetting blocks'
    self.setSpaceFlag
    self.setTimeFlag
    self.setErrorFlag
    self.resetBlocks
  ENDIF

  latsran=(self.params.latlims[1]-self.params.latlims[0])/self.params.nlats*dindgen(self.params.nlats+1)+self.params.latlims[0]
  longsran=(self.params.lonlims[1]-self.params.lonlims[0])/self.params.nlons*dindgen(self.params.nlons+1)+self.params.lonlims[0]
  altsran=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts+1)+self.params.altlims[0]
  
  dayLast=julday(self.params.monthlims[1],self.params.daylims[1],self.params.yearlims[1],self.params.hourlims[1],0,0)
  dayFirst=julday(self.params.monthlims[0],self.params.daylims[0],self.params.yearlims[0],self.params.hourlims[0],0,0)
  daysran=(dayLast-dayFirst)/self.params.nTimes*dindgen(self.params.nTimes+1)+dayFirst
  ;stop
  IF self.params.solidBlock EQ 1 THEN BEGIN
    for l=0,self.params.nTimes-1 do begin
      jlims=daysran[l:l+1]
      for k=0,self.params.nAlts-1 do begin
        alims=altsran[k:k+1]
        ;stop
        for j=0,self.params.nLons-1 do begin
          lonlims=longsran[j:j+1]
          for i=0,self.params.nLats-1 do begin
            latlims=latsran[i:i+1]
            ;stop
            self.setOpMask,latlims,lonlims,alims,jdays=jlims
            self.reduceBlock,i,j,k,l
            (*self.blocks)[i,j,k,l].yearlims=[-1,-1]
            (*self.blocks)[i,j,k,l].monthlims=[-1,-1]
            (*self.blocks)[i,j,k,l].daylims=[-1,-1]
            (*self.blocks)[i,j,k,l].hourlims=[-1,-1]
            (*self.blocks)[i,j,k,l].jullims=jlims
            (*self.blocks)[i,j,k,l].latlims=latlims
            (*self.blocks)[i,j,k,l].lonlims=lonlims
            (*self.blocks)[i,j,k,l].altlims=alims
            ;stop

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
    for l=0,self.params.nTimes-1 do begin
      ylims=[timesran[l,0],timesran[l+1,0]]
      mlims=[timesran[l,1],timesran[l+1,1]]
      dlims=[timesran[l,2],timesran[l+1,2]]
      hlims=[timesran[l,3],timesran[l+1,3]]
      for k=0,self.params.nAlts-1 do begin
        alims=altsran[k:k+1]
        ;stop
        for j=0,self.params.nLons-1 do begin
          lonlims=longsran[j:j+1]
          for i=0,self.params.nLats-1 do begin
            latlims=latsran[i:i+1]
            self.setOpMask,latlims,lonlims,alims,byears=ylims,bmonths=mlims,bdays=dlims,bhours=hlims
            self.reduceBlock,i,j,k,l
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
        stop
      ;ENDFOR
      
    ENDFOR
  ENDELSE

END

PRO ACEData::setAxes,xaxis,yaxis



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

PRO ACEData::makeContour,ratPlot,mol0Plot,aLevels,useMedian=useMedian

  alt_map=aLevels             ;currently only one level at a time is supported

  ct=colortable(72,/reverse)
  ;ct=colortable(25)
  
  ;-------- This block defines a few variables for the plots. The case structure
  ;defines the color scale base on altitude, and could probably use further tweaks. --------;
  posit=[0.05,0.1,0.87,0.96]
  dim=[1280,580]
  case (*self.hmz).name0 of
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

  endcase

  ;-------- This section plots the ratio and primary molecule data. If savPlotRat or savPlot0 are set, then the plots are saved
  ;directly to file.
  
  lats=(dindgen(self.params.nLats)+0.5d0)/self.params.nLats*(self.params.latLims[1]-self.params.latLims[0])+self.params.latLims[0]
  lons=(dindgen(self.params.nLons)+0.5d0)/self.params.nLons*(self.params.lonLims[1]-self.params.lonLims[0])+self.params.lonLims[0]
  
  if (*self.hmz).nmols gt 1 then begin
    if keyword_set(useMedian) then rats=(*self.blocks).bwRat[2] else rats=(*self.blocks).meanRat
    if keyword_set(savPlotRat) then mapRat=map('Geographic',center_longitude=180,limit=[-65,-180,65,180],dimensions=dim,position=posit,font_size=18,/buffer) else mapRat=map('Geographic',center_longitude=180,limit=[-65,-180,65,180],dimensions=dim,position=posit,font_size=18)
    ;stop
    ratios=interpolate(transpose((rats-1d0)*1d3),dindgen(self.params.nLons),dindgen(self.params.nLats),/grid)
    ;cc=contour((rats-1d0)*1d3,(*self.blocks).meanLon,(*self.blocks).meanLat,c_value=([dindgen(17)]/16d0*.2d0-.76d0)*1d3,rgb_table=ct,/fill,/overplot)
    cc=contour(ratios,lons,lats,c_value=([dindgen(17)]/16d0*.46d0-.76d0)*1d3,rgb_table=ct,/fill,/overplot)
    mc=mapcontinents(/continents,limit=[-65,-180,65,180])
    grid=mapRat.mapgrid
    grid.linestyle='dotted'
    grid.font_size=16
    grid.label_position=0
    mapRat.title=self.params.mth_str+' $\delta D$, '+string(alt_map,format='(F4.1)')+' km'
      mapRat.font_size=18
    cb=colorbar(target=cc,title='$\delta D (‰)$',range=[-.760d0,-.500d0]*1d3,orientation=1,textpos=1,font_size=14)
    ratPlot=mapRat
    if keyword_set(savPlotRat) then begin
      ratPlot.save,savPlotRat,resolution=300
      ratPlot.close
    endif
    ;stop
  endif
  if keyword_set(savPlot0) then map0=map('Geographic',center_longitude=90,limit=[-65,-180,65,180],dimensions=dim,position=posit,font_size=18,/buffer) else map0=map('Geographic',center_longitude=90,limit=[-65,-180,65,180],dimensions=dim,position=posit,font_size=18)
  if keyword_set(useMedian) then mol0=(*self.blocks).bw0[2] else mol0=(*self.blocks).mean0
  mol0s=interpolate(transpose(mol0*1d6),dindgen(self.params.nLons),dindgen(self.params.nLats),/grid)
  ;cc=contour(mol0*1d6,(*self.blocks).meanLon,(*self.blocks).meanLat,c_value=mol0range,rgb_table=ct,/fill,/overplot)
  cc=contour(mol0s,lons,lats,c_value=mol0range,rgb_table=ct,/fill,/overplot)
  mc=mapcontinents(/continents,limit=[-65,-180,65,180])
  grid=map0.mapgrid
  grid.linestyle='dotted'
  grid.font_size=16
  grid.label_position=0
  map0.title=self.params.mth_str+' $H_2O$, '+string(alt_map,format='(F4.1)')+' km'
    map0.font_size=18
  cb=colorbar(target=cc,title=(*self.hmz).name0+' (ppm)',range=[mol0range[0],mol0range[n_elements(mol0range)-1]],orientation=1,textpos=1,font_size=14)
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

PRO ACEData::makeMeridianPlot,ratPlot,mol0Plot,useMedian=useMedian

  ct=colortable(72,/reverse)

  if keyword_set(useMedian) then rats=(*self.blocks).bwRat[2] else rats=(*self.blocks).meanRat
  n_longs=self.params.nLons
  alts=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts)+self.params.altlims[0]+.5
  alts=transpose(rebin(alts,self.params.nAlts,self.params.nLats))
  ;lats=(self.params.latlims[1]-self.params.latlims[0])/self.params.nLats*dindgen(self.params.nLats)+self.params.latlims[0]
  cbrange=([dindgen(22)]/21d0*0.63d0-.76d0)*1d3
  for i=0,n_longs-1 do begin
    lats=reform((*self.blocks)[*,i,*].meanLat,self.params.nLats,self.params.nAlts)
    ratPlot=contour((reform(rats[*,i,*],self.params.nLats,self.params.nAlts)-1)*1d3,lats,alts,c_value=cbrange,xrange=latlims,yrange=[8,40],font_size=18,axis_style=2,xthick=2,ythick=2,xtitle='Latitude',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=[1024,1024],position=[0.1,0.1,0.85,0.9],/fill)
    cb=colorbar(target=ratPlot,title='$\delta D (‰)$',range=[-.760,-.495]*1d3,orientation=1,textpos=1,font_size=14)

  endfor

  if keyword_set(useMedian) then mol0=(*self.blocks).bw1[2] else mol0=(*self.blocks).mean0
  n_longs=self.params.nLons
  alts=(self.params.altlims[1]-self.params.altlims[0])/self.params.nAlts*dindgen(self.params.nAlts)+self.params.altlims[0]+.5
  alts=transpose(rebin(alts,self.params.nAlts,self.params.nLats))
  ;lats=(self.params.latlims[1]-self.params.latlims[0])/self.params.nLats*dindgen(self.params.nLats)+self.params.latlims[0]
  cbrange=[dindgen(21)]/20d0*12d0+2d0
  for i=0,n_longs-1 do begin
    lats=reform((*self.blocks)[*,i,*].meanLat,self.params.nLats,self.params.nAlts)
    mol0Plot=contour((reform(mol0[*,i,*]*1d6,self.params.nLats,self.params.nAlts)-1),lats,alts,c_value=cbrange,xrange=latlims,yrange=[8,40],font_size=18,axis_style=2,xthick=2,ythick=2,xtitle='Latitude',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=[1024,1024],position=[0.1,0.1,0.85,0.9],/fill)
    cb=colorbar(target=mol0Plot,title='$H_2O (ppm)$',range=[3,15],orientation=1,textpos=1,font_size=14)

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

PRO ACEData::makeTimeSeries,ratPlot,mol0Plot,useMedian=useMedian

  ct=colortable(72,/reverse)

  if keyword_set(useMedian) then rats=(*self.blocks).bwRat[2] else rats=(*self.blocks).meanRat
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

  ct=colortable(72,/reverse)


END

PRO ACEData::makeHeightTimePlot,mol0,rat0,alt_map=alt_map,yearlims=yearlims,monthlims=monthlims,hourlims=hourlims,latlims=latlims,lonlims=lonlims,n_lats=n_lats,n_longs=n_longs,ratio=ratio,ratlims=ratlims,raw=raw,acecut=acecut,fullcut=fullcut,xrange=xrange

  IF NOT KEYWORD_SET(xrange) THEN xrange=[]
  
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
  mol0=contour(transpose(reform((*self.blocks).bw0[2],nAlts,nTimes))*1d6,times,alts,c_value=[dindgen(22)]*.25+3d0,xrange=xrange,yrange=[10,40],font_size=18,axis_style=2,xthick=2,ythick=2,xtitle='Year',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=[1024,1024],position=[0.1,0.1,0.85,0.9],/fill)
  cb=colorbar(target=mol0,title='$H_2O (ppm)$',range=[3,8.25],orientation=1,textpos=1,font_size=14)

  rat0=contour((transpose(reform((*self.blocks).bwrat[2],nAlts,nTimes))-1d0)*1d3,times,alts,c_value=([dindgen(22)]/21d0*0.63d0-.76d0)*1d3,xrange=xrange,yrange=[10,40],font_size=18,axis_style=2,xthick=2,ythick=2,xtitle='Year',ytitle='Altitude (km)',rgb_table=colortable(72,/reverse),dimensions=[1024,1024],position=[0.1,0.1,0.85,0.9],/fill)
  cb=colorbar(target=rat0,title='$\delta D (‰)$',range=[-.760,-.130]*1d3,orientation=1,textpos=1,font_size=14)
  ;stop
END

pro ACEData::makeMeanSD,xname,yname,spag,savPlot=savPlot,spaghetti=spaghetti,colorby=colorby,outcols=outcols

  ct=colortable(72,/reverse)
      nspg=n_elements(insrat)
      
      CASE yname OF
        'alt':BEGIN
          ymean=self.params.altlims[0]+dindgen(self.params.nAlts)+0.5d0
          ytitle='Altitude (km)'
          yrange=self.params.altlims
          cmdy='orbs.alt[jk]'
        END
        'pt':BEGIN
          ymean=reform((*self.blocks).meanpt,self.params.nAlts)
          ytitle='$\Theta (K)$'
          yrange=[300,500]
          cmdy='orbs.pt[jk]*1d6'
        END
        'mol0':BEGIN
          ymean=reform((*self.blocks).mean0*1d6,self.params.nAlts)
          ysd=reform((*self.blocks).sd0*1d6,self.params.nAlts)
          ybw=reform((*self.blocks).bw0,5,self.params.nAlts)
          ytitle='$H_2O (ppm)$'
          yrange=[0,100]
          cmdy='orbs.mixrat0[jk,is]*1d6'
        END
        'temp':BEGIN
          ymean=reform((*self.blocks).meant,self.params.nAlts)
          ysd=reform((*self.blocks).sdt,self.params.nAlts)
          ybw=reform((*self.blocks).bwt,5,self.params.nAlts)
          ytitle='Temp. (K)'
          yrange=[180,270]
          cmdy='orbs.temp[jk,is]'
        END
      ENDCASE

      CASE xname OF
        'delta':BEGIN
          xmean=reform(((*self.blocks).meanrat-1d0)*1d3,self.params.nAlts)
          xsd=reform(((*self.blocks).sdrat)*1d3,self.params.nAlts)
          xbw=reform(((*self.blocks).bwrat-1d0)*1d3,5,self.params.nAlts)
          xtitle='$\delta D (‰)$'
          xrange=[-2000,1000]
          cmdx='(orbs.ratio[jk,is]-1d0)*1d3'
        END
        'mol0':BEGIN
          xmean=reform((*self.blocks).mean0*1d6,self.params.nAlts)
          xsd=reform((*self.blocks).sd0*1d6,self.params.nAlts)
          xbw=reform((*self.blocks).bw0,5,self.params.nAlts)
          xtitle='$H_2O (ppm)$'
          xrange=[0,500]
          cmdx='orbs.mixrat0[jk,is]*1d6'
        END
        'mol1':BEGIN
          xmean=reform((*self.blocks).mean1*1d6,self.params.nAlts)
          xsd=reform((*self.blocks).sd1*1d6,self.params.nAlts)
          xbw=reform((*self.blocks).bw1,5,self.params.nAlts)
          xtitle='$HDO (ppm)$'
          xrange=[0,500]
          cmdx='orbs.mixrat1[jk,is]*1d6'
        END
        'temp':BEGIN
          xmean=reform((*self.blocks).meant,self.params.nAlts)
          xsd=reform((*self.blocks).sdt,self.params.nAlts)
          xbw=reform((*self.blocks).bwt,5,self.params.nAlts)
          xtitle='Temp. (K)'
          xrange=[180,270]
          cmdx='orbs.temp[jk,is]'
        END
        'pres':BEGIN
          xmean=reform((*self.blocks).meanp,self.params.nAlts)
          xsd=reform((*self.blocks).sdp,self.params.nAlts)
          xbw=reform((*self.blocks).bwp,5,self.params.nAlts)
          xtitle='Pres. (hPa)'
          xrange=[0,1]
          cmdx='orbs.pres[jk,is]'
        END
      ENDCASE

      if keyword_set(savPlot) then spag=boxplot(ymean,xbw,xtitle=xtitle,ytitle=ytitle,font_size=16,xthick=2,ythick=2,/horizontal,xrange=xrange,yrange=yrange,dimensions=[640,720],/buffer) else spag=boxplot(ymean,xbw,xtitle=xtitle,ytitle=ytitle,font_size=16,xthick=2,ythick=2,/horizontal,xrange=xrange,yrange=yrange,dimensions=[640,720])
      map=map('Geographic',center_longitude=90,limit=[-75,-180,75,180],position=[0.63,0.82,.93,1.06],font_size=18,/current,axis_style=0,linestyle=6,label_show=0)
      mc=mapcontinents(/continents,limit=[-75,-180,75,180])
      coords=[[(*self.blocks)[0].lonlims[0],(*self.blocks)[0].latlims[1]],[(*self.blocks)[0].lonlims[1],(*self.blocks)[0].latlims[1]],[(*self.blocks)[0].lonlims[1],(*self.blocks)[0].latlims[0]],[(*self.blocks)[0].lonlims[0],(*self.blocks)[0].latlims[0]]]
      region=polygon(coords,target=map,/data,fill_color='green')
      spag.select
      meanRat=plot(xmean,ymean,thick=2,'r',/overplot)
      sdRat=polygon([xmean-xsd,reverse(xmean+xsd)],[ymean,reverse(ymean)],/data,/fill_background, FILL_COLOR="light steel blue",PATTERN_ORIENTATION=45,pattern_spacing=4,/overplot)

      if keyword_set(savPlot) then begin
        spag.save,savPlot+mth_str+'_Rat '+string(mean(mols.lat[jrat[insrat]]),format='(F7.2)')+' X '+string(mean(mols.lon[jrat[insrat]]),format='(F7.2)')+'.png',resolution=300
      endif


  if keyword_set(spaghetti) then begin
    IF KEYWORD_SET(outcols) THEN outcol=LIST()
    maxPts=MAX((*self.blocks).npts,mx)
    norbs=(*self.blocks)[mx].npts
    dex=(*(*self.blocks)[mx].orbids)
    orbs=self.getOcc(dex)
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
      jk=where((*self.flags).aflag[*,dex[is]] eq 1)
      jk1=where(ymean ge orbs.alt[jk[0]])
      
      cmd='p=plot('+cmdx+','+cmdy+',transparency=85,color=color,/overplot)'
      res=execute(cmd)
      
      IF KEYWORD_SET(outcols) THEN BEGIN
        onecall=DBLARR(3,N_ELEMENTS(jk))
        onecall[0,*]=orbs.alt[jk]
        onecall[1,*]=orbs.mixrat0[jk,is]
        onecall[2,*]=(orbs.ratio[jk,is]-1d0)*1d3
        ;stop
        outcol.add,onecall
      ENDIF

    endfor
    
    IF KEYWORD_SET(outcols) THEN outcols=outcol.toarray(dimension=2)
    IF KEYWORD_SET(colorby) THEN cb=colorbar(target=spag,range=[self.params.latlims[0],self.params.latlims[1]],position=[0.03,0.0,0.35,0.03],title=colorby,font_size=10,textpos=1,rgb_table=ct)
  endif

end


pro ACEData::makeMeanSDTP,xname,yname,savPlot=savPlot


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


PRO ACEData::stopView

  stop
  
END

PRO ACEData__define

  A={paramsA,yearlims:[0.,0.],monthlims:[0.,0.],daylims:[0.,0.],hourlims:[0.,0.],jullims:[0.,0.],latlims:[0.,0.],lonlims:[0.,0.],altlims:[0.,0.],nTimes:0.,nLats:0.,nLons:0.,nAlts:0.,oParamName:'',oParamLims:[0.,0.],oParamNum:0.,mth_str:'',elim:0d0,cuttype:'',tMod:boolean(0),sMod:boolean(0),eMod:boolean(0),bMod:boolean(0),solidBlock:boolean(0),dmp_path:''}
  struct={ACEData,raw:PTR_NEW(),hmz:PTR_NEW(),flags:PTR_NEW(),opMask:PTR_NEW(),blocks:PTR_NEW(),dmp:PTR_NEW(),params:A}
  
end