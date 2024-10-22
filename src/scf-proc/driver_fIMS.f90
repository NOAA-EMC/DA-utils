 program driver_fIMS

! calculate fractional snow cover on model grid, from IMS snow cover obs. 
! then calculate SWE from fractional snow cover, assuming the snow 
! depletion curve used by the noah model. 
! 
! Clara Draper, July 2021 (based on code from Tseganeh Gichamo, Youlong Xia)

 use IMSaggregate_mod, only: calculate_scfIMS

 implicit none

 integer             :: idim, jdim
 character(len=11)   :: yyyymmddhh
 character(len=7)    :: jdate
 character(len=20)   :: otype ! orography type, format C$RES (atm) or C${RES}.mx100 (coupled atm/ocean)
 character(len=200)  :: IMS_obs_path, IMS_ind_path, fcst_path
 logical             :: file_exists, skip_SD
 integer             :: io, ierr, lsm, imsformat
 character(len=10)   :: imsversion, imsres

 namelist/fIMS_nml/  idim, jdim, otype, yyyymmddhh, jdate, IMS_obs_path, IMS_ind_path, fcst_path, lsm, imsformat, imsversion, imsres, skip_SD

 ! default to current directory 
 IMS_obs_path="./"
 IMS_ind_path="./"
 fcst_path="./"

 lsm=2 ! 1 - noah, 2 - noah-MP. 
 imsformat=2  ! 1 - ascii format, 2 - netCDF format
 imsversion="1.3" ! 1.2 before Dec 3 2014, 1.3 from Dec 3 2014 onwards.
 imsres="4km"
 skip_SD=.false.

 ! read namelist
 inquire(file='fims.nml', exist=file_exists)

 if (.not. file_exists) then
        print *, 'namelistfile does not exist, exiting' 
        stop 10
 endif

 open (action='read', file='fims.nml', iostat=ierr, newunit=io)
 read (nml=fIMS_nml, iostat=ierr, unit=io) 
 close (io) 
 
 print *, 'lsm', lsm
 call calculate_scfIMS(idim, jdim, otype,yyyymmddhh, jdate,IMS_obs_path, IMS_ind_path, & 
                                fcst_path, lsm, imsformat, imsversion, imsres, skip_SD)

 end program driver_fIMS
