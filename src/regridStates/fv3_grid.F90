 module fv3_grid
 ! ESMF grid-specific routines for GFS regridding program, including IO.
 ! 
 ! Clara Draper, Aug 2024.

 use esmf
 use netcdf
 use ESMF_LogPublicMod
 use utilities, only     : error_handler, netcdf_err

 implicit none

 private

 integer, public, parameter  :: n_tiles=6 ! number tiles in fv3 grid
 integer, public, parameter  :: vtype_water=0, & ! TO DO - which veg classification is this?
                                vtype_landice=15 ! used for soil mask

 public :: setup_grid, &
           read_into_fields, &
           write_from_fields

 contains

!-----------------------------------
! Create ESMF grid objects, with mask if requested
 subroutine setup_grid(localpet, npets, dir_fix, res_atm, & 
                       dir_rst, fname_rst, mask_type, fv3_grid)

 implicit none

 ! INTENT IN
 character(*), intent(in)       :: dir_fix
 integer, intent(in)            :: res_atm
 character(*), intent(in)       :: fname_rst, dir_rst
 character(*), intent(in)       :: mask_type
 integer, intent(in)            :: localpet, npets

 ! INTENT OUT 
 type(esmf_grid)                :: fv3_grid

 ! LOCAL
 type(esmf_field)               :: vtype_field(1)
 real(esmf_kind_r8), pointer    :: ptr_vtype(:,:)
 integer(esmf_kind_i4), pointer :: ptr_mask(:,:)

 character(len=500)             :: fname, dir_fix_res
 integer                        :: extra, ierr, ncid, tile
 integer                        :: decomptile(2,n_tiles)
 character(len=10)              :: variable_list(1)
 character(len=5)               :: rchar


 if (localpet == 0) print*," creating grid for ", res_atm 

! pet distribution
 extra = npets / n_tiles

 do tile = 1, n_tiles
   decomptile(:,tile)=(/1,extra/)
 enddo

 ! mosaic file 
 write(rchar,'(i3.3)') res_atm

 dir_fix_res = dir_fix//"/C"//trim(rchar)//"/"

 fname = trim(dir_fix_res)//"/C"//trim(rchar)// "_mosaic.nc"
 
! create the grid
 fv3_grid = ESMF_GridCreateMosaic(filename=trim(fname), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(dir_fix_res), &
                                  rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateMosaic", ierr)

!--------------------------
! Calcalate and add the mask

 if (mask_type=="soil") then ! mask out ocean and glaciers 
     vtype_field(1) = ESMF_FieldCreate(fv3_grid, &
                                       typekind=ESMF_TYPEKIND_R8, &
                                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                                       name="input veg type for mask", &
                                       rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldCreate, vtype", ierr)

     variable_list(1) = 'vtype     '
     call read_into_fields(localpet, res_atm, res_atm, fname_rst, 1, variable_list(1), &
                           dir_rst, vtype_field(1))

    ! get pointer to vegtype
     call ESMF_FieldGet(vtype_field(1), &
                        farrayPtr=ptr_vtype, &
                        rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGet", ierr)

    ! create and get pointer to the mask
     call ESMF_GridAddItem(fv3_grid, &
                           itemflag=ESMF_GRIDITEM_MASK, &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("in GridAddItem mask", ierr)

     call ESMF_GridGetItem(fv3_grid, & 
                           itemflag=ESMF_GRIDITEM_MASK, &
                           farrayPtr=ptr_mask, &
                           rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("in GridGetItem mask", ierr)

    ! calculate the mask
     ptr_mask = 1 ! initialize land everywhere
     where (nint(ptr_vtype) == vtype_water )   ptr_mask = 0 ! exclude water
     where (nint(ptr_vtype) == vtype_landice ) ptr_mask = 0 ! exclude glaciers

    ! destroy veg type field
     call ESMF_FieldDestroy(vtype_field(1),rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("DESTROYING FIELD", ierr)

  end if ! mask = soil 

 end subroutine setup_grid


 ! read variables from fv3 netcdf restart file into ESMF Fields
 subroutine read_into_fields(localpet, i_dim, j_dim , fname_read, n_vars, variable_list, & 
                             dir_read, fields) 
 
 implicit none 

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim, n_vars
 character(*), intent(in)        :: fname_read
 character(*), intent(in)        :: dir_read

 character(len=10), dimension(n_vars), intent(in)   :: variable_list
 
 ! INTENT OUT 
 type(esmf_field), dimension(n_vars), intent(inout) :: fields

 ! LOCAL
 integer                         :: tt, id_var, ncid, ierr, v
 character(len=1)                :: tchar
 character(len=500)              :: fname
 real(esmf_kind_r8), allocatable :: array2D(:,:)
 real(esmf_kind_r8), allocatable :: array_in(:,:,:)

 if (localpet==0) then
     allocate(array_in(n_vars,i_dim, j_dim))
     allocate(array2D(i_dim, j_dim))
 else 
     allocate(array_in(0,0,0))
     allocate(array2D(0,0))
 end if

 do tt = 1, n_tiles

      ! read from restart
      if (localpet == 0) then

         write(tchar,'(i1)') tt
         fname = dir_read//"/"//fname_read//".tile"//tchar//".nc"

         ierr=nf90_open(trim(fname),NF90_NOWRITE,ncid)
         call netcdf_err(ierr, 'opening: '//trim(fname) )

         do v =1, n_vars
             print *, 'Reading ', trim(variable_list(v)), ' into field, tile', tchar
             ierr=nf90_inq_varid(ncid, trim(variable_list(v)), id_var)
             call netcdf_err(ierr, 'reading variable id' )

             ierr=nf90_get_var(ncid, id_var, array_in(v,:,:))
             call netcdf_err(ierr, 'reading variable' )
         enddo
         ierr = nf90_close(ncid)

      endif

      ! scatter
      do v =1, n_vars
          array2D=array_in(v,:,:) ! scatter misbehaves if given indexed 3D array.
          call ESMF_FieldScatter(fields(v), array2D, rootpet=0, tile=tt, rc=ierr)
          if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldScatter", ierr)
      enddo

 enddo
 
 ! clean up
 deallocate(array_in)

 end subroutine read_into_fields

! write variables from ESMF Fields into netcdf restart-like file 

 subroutine write_from_fields(localpet, i_dim, j_dim , fname_out, n_vars, variable_list, & 
                          dir_out, fields) 

 implicit none 

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim, n_vars
 character(*), intent(in)        :: fname_out
 character(*), intent(in)        :: dir_out
 character(10), dimension(n_vars), intent(in)     :: variable_list
 type(esmf_field), dimension(n_vars), intent(in)  :: fields

 ! LOCAL
 integer                         :: tt, id_var, ncid, ierr, &
                                   id_x, id_y, v
 character(len=1)                :: tchar
 character(len=500)              :: fname
 real(esmf_kind_r8), allocatable :: array2D(:,:)
 real(esmf_kind_r8), allocatable :: array_out(:,:,:)

 do v = 1, n_vars
        if (localpet == 0)  print *, 'Writing ', trim(variable_list(v)), ' into field'
 enddo

 if (localpet==0) then
     allocate(array_out(n_vars, i_dim, j_dim))
     allocate(array2D(i_dim, j_dim))
 else 
     allocate(array_out(0,0,0))
     allocate(array2D(0,0))
 end if

 do tt = 1, n_tiles

      ! fetch data (all PETs)
      do  v=1, n_vars
          call ESMF_FieldGather(fields(v), array2D, rootPet=0, tile=tt, rc=ierr)
          if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldGather", ierr)
          array_out(v,:,:) = array2D
      enddo

      ! write to netcdf 
      if (localpet == 0) then

         ! open file, set dimensions
         write(tchar,'(i1)') tt
         fname = dir_out//"/"//fname_out//".tile"//tchar//".nc"

         ierr = nf90_create(trim(fname), NF90_NETCDF4, ncid)
         call netcdf_err(ierr, 'creating file='//trim(fname) )

         ierr = nf90_def_dim(ncid, 'xaxis_1', i_dim, id_x)
         call netcdf_err(ierr, 'defining xaxis dimension' )

         ierr = nf90_def_dim(ncid, 'yaxis_1', j_dim, id_y)
         call netcdf_err(ierr, 'defining yaxis dimension' )

         do v=1, n_vars

             ierr = nf90_def_var(ncid, trim(variable_list(v)), NF90_DOUBLE, (/id_x, id_y/), id_var)
             call netcdf_err(ierr, 'defining '//variable_list(v) )
              
             ierr = nf90_put_var( ncid, id_var, array_out(v,:,:) ) 
             call netcdf_err(ierr, 'writing '//variable_list(v) ) 

         enddo

         ierr = nf90_close(ncid)

      endif

 enddo
 
 ! clean up
 deallocate(array_out)

 end subroutine write_from_fields

end  module fv3_grid
