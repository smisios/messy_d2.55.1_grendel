MODULE messy_clamsrdfrc
!***********************************************************************!
! CLaMS radiative forcing in EMAC
! Interpolation from CLaMS airparcel postitions to EMAC grid (i3d) 
!***********************************************************************!
  USE messy_clams_global,    ONLY: PREC, nx, ny

  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamsrdfrc'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

  LOGICAL :: use_o3forc, use_h2oforc
  INTEGER :: timestep_rdfrc
  INTEGER :: bound_trop_n, bound_trop_s, lev_bound_trop, lev_bound_extr
  REAL(PREC):: bound_up_zeta

  PUBLIC  :: clamsrdfrc_read_nml
  PUBLIC  :: rdfrc
  PRIVATE :: nc_read_input_data
  PRIVATE :: int3d

!----------- 
CONTAINS
!----------- 

!-----------------------------------------------------------------------
  SUBROUTINE rdfrc(LAT, LON, LEV, CLAMS_DATA, ECHAM_DATA, E5_ZETA, varname)

    USE messy_clamsrdfrc_tools 
    USE messy_clamsrdfrc_global, ONLY: ap, datetype, ctrl_out, max_dist,&
         & ap_g_all, npart_g, ntags, ap_s_all, data_s_all, max_nb, dim_tr,&
         & lev_window, data_g_tot
    USE messy_clams_global, ONLY: init_vertcoorname
    USE messy_clams_tools_utils,  ONLY: lowercase
    USE messy_clamsmix_global,      ONLY: adapt_par
    USE messy_clams_global, ONLY: mdi, nparts_max
    USE messy_clams_tools_triang, ONLY: find_triangle, triang_qhull
    USE messy_clams_tools_dateconv, ONLY: ymds2js
    USE messy_clamsbmix_global, ONLY: nlevs

    IMPLICIT NONE

    ! Arguments
    REAL(PREC), DIMENSION(:),     POINTER :: LAT, LON, LEV, CLAMS_DATA
    REAL(PREC), DIMENSION(:,:,:), INTENT(INOUT) :: ECHAM_DATA
! op_pj_20160830+
! POINTER and INTENT attribute specifier cannot be specified at the same time.
!!$ REAL(PREC), DIMENSION(:,:,:), POINTER,INTENT(IN) :: E5_ZETA
    REAL(PREC), DIMENSION(:,:,:), POINTER :: E5_ZETA
! op_pj_20160830-
    CHARACTER(LEN=*)                      :: varname

    ! Local variables
    TYPE(ap),       DIMENSION(:),   POINTER :: ap_g, ap_g_tot, ap_g_dummy, &
                                               ap_s_down, ap_s_up
    REAL(PREC),     DIMENSION(:,:), POINTER :: vert_g_coor, vert_s_down_coor, vert_s_up_coor
    INTEGER,        DIMENSION(:),   POINTER :: vert_s_down_nb, vert_s_up_nb
    INTEGER,        DIMENSION(:,:), POINTER :: vert_s_down_ind, vert_s_up_ind
    INTEGER,        DIMENSION(:,:), POINTER :: triangle_s_down, triangle_s_up
    REAL(PREC),     DIMENSION(:,:), POINTER :: data_g, data_g_dummy, &
                                               data_s_down, data_s_up
    INTEGER,        DIMENSION(:),   POINTER :: t_tot_down, t_tot_up
    CHARACTER(40)                           :: vertcoorname, vertcoorname_clams
    REAL(PREC),     DIMENSION(:),   POINTER :: lev_grid, lev_grid2
    INTEGER,        DIMENSION(:),   POINTER :: indexarr, indexarr_dummy
    LOGICAL,        DIMENSION(:),   POINTER :: mask
    LOGICAL,        DIMENSION(:),   POINTER :: currentday
    TYPE(datetype)                          :: init_date
    REAL(PREC)                              :: delta_lev_up, delta_lev_down
    INTEGER                                 :: nlevs_loop, ilev, itag, sec, &
                                              index_down, index_up, i, t_act,&
                                              & npartg_act, nparts_act, incday
    LOGICAL                                 :: data_present, firsttime 
    REAL(KIND(0d0))  :: js1, js2
    INTEGER :: nlev, nlat, nlon
    INTEGER :: status

  ctrl_out = .TRUE.
  ntags = 1

  !---------------------------------------------------------
  ! max. distance [km]
  !---------------------------------------------------------
  max_dist = 2000.

  ! vertical coordinate in CLaMS file must be 'zeta'
  IF (TRIM(lowercase(init_vertcoorname))/='zeta') THEN  
     WRITE (*,*) ' Error in clamsrdfrc!!'
     WRITE (*,*) 'Vertical coordinate in CLaMS file not zeta !!! '
     WRITE (*,*) 'Vertical coordinate in CLaMS file: ',TRIM(lowercase(vertcoorname_clams))
     WRITE (*,*)
     STOP
  ENDIF
  vertcoorname_clams = init_vertcoorname

  !---------------------------------------------------------
  ! read lev_grid (CLaMS)
  !---------------------------------------------------------
  lev_grid2 => adapt_par%lev_grid

  IF (nlevs==1) THEN ! 2d interpolation
     ALLOCATE (lev_grid(2))
     lev_grid(1) = 0.
     lev_grid(2) = 10000.
  ELSE
     ALLOCATE (lev_grid(nlevs+2))
     lev_grid(1) = 0.
     lev_grid(2:nlevs+1) = lev_grid2 
     lev_grid(nlevs+2) = 10000.
  ENDIF
  WRITE (*,*) 'in rdfrc: lev_grid=',lev_grid

  !---------------------------------------------------------
  !  read input data (lat, lon, lev)   
  !---------------------------------------------------------
  WRITE (*,*)
  WRITE (*,*) 'read positions (input/output file)'

  CALL nc_read_input_data (nlev, nlat, nlon, npart_g, ap_g_all, E5_ZETA)
     
  ALLOCATE (data_g_tot(ntags, npart_g))
  ALLOCATE (ap_g_tot(npart_g))

  ! initialize  data_g_tot and ap_g_tot  
  data_g_tot = mdi
  ap_g_tot(:)%lon = mdi
  ap_g_tot(:)%lat = mdi
  ap_g_tot(:)%lev = mdi
    
  !---------------------------------------------------------
  ! read positions of APs and their dynamical properties
  !---------------------------------------------------------
  WRITE (*,*)
  WRITE (*,*) 'read positions (CLaMS)'
  CALL nc_read_all_ap (LAT, LON, LEV, ap_s_all)
  WRITE (*,*)
  WRITE (*,*) 'read dynamical properties (CLaMS file)'

! Write O3/H2O to data_s_all
  allocate (data_s_all(ntags,nparts_max))

  data_s_all(1,:) = CLAMS_DATA
   
  firsttime = .TRUE.
  
  NULLIFY (triangle_s_up)
  NULLIFY (ap_s_up)
  NULLIFY (data_s_up)
  NULLIFY (vert_s_up_coor)
  NULLIFY (vert_s_up_nb)
  NULLIFY (vert_s_up_ind)
  
  IF (nlevs==1) THEN ! 2d interpolation
     nlevs_loop = 1
  ELSE
     nlevs_loop = nlevs+1
  ENDIF
  
  !----------------------------------------------------------------
  ! loop over all lev levels
  !----------------------------------------------------------------
  DO ilev = 1, nlevs_loop
        
     !----------------------------------------------------------------
     ! 
     !----------------------------------------------------------------
     IF (ilev == 1) THEN
        index_down = 1 
        index_up = 1
     ELSEIF (ilev == nlevs+1) THEN
        index_down = nlevs 
        index_up = nlevs 
     ELSE
        index_down = ilev-1
        index_up = ilev
     END IF
     
     WRITE (*,*) TRIM(lowercase(vertcoorname_clams)),' layer: ', &
          lev_grid(ilev), lev_grid(ilev+1)
     WRITE (*,*) 'lev window: ',lev_window(index_down,1), lev_window(index_up,2)
     
     delta_lev_down = lev_window(index_down,2) - lev_window(index_down,1)
     delta_lev_up   = lev_window(index_up,2)   - lev_window(index_up,1)
     
     !----------------------------------------------------------------
     ! get points between lev_grid(ilev) and lev_grid(ilev+1)
     !----------------------------------------------------------------
     CALL get_apg (ap_g, indexarr, lev_window(index_down,1), lev_window(index_up,2), &
          data_present)
        
     IF (.NOT. data_present) THEN 
        WRITE (*,*) 'keine Daten !!!'
        firsttime = .TRUE.
        CYCLE
     ENDIF
        
     npartg_act = SIZE(ap_g)
     WRITE (*,*) 'number of points on ', TRIM(lowercase(vertcoorname_clams)),' layer: ', npartg_act
        
     !----------------------------------------------------------------
     ! allocate arrays
     !----------------------------------------------------------------
     ALLOCATE (data_g(ntags, npartg_act))
     ALLOCATE (vert_g_coor(3,npartg_act))
     ALLOCATE (t_tot_down(npartg_act))
     ALLOCATE (t_tot_up  (npartg_act))
        
     !----------------------------------------------------------------
     ! Transf ap-coordinates (lon, lat) on kart. coord. on a unit sph.
     !----------------------------------------------------------------
     CALL set_coor (ap_g, vert_g_coor)
     
     !----------------------------------------------------------------
     ! lower level
     !----------------------------------------------------------------
     IF (firsttime) THEN 
           
        IF (ASSOCIATED(triangle_s_up)) DEALLOCATE (triangle_s_up)
        IF (ASSOCIATED(ap_s_up)) DEALLOCATE (ap_s_up)
        IF (ASSOCIATED(data_s_up)) DEALLOCATE (data_s_up)
        IF (ASSOCIATED(vert_s_up_coor)) DEALLOCATE (vert_s_up_coor)
        IF (ASSOCIATED(vert_s_up_nb)) DEALLOCATE (vert_s_up_nb)
        IF (ASSOCIATED(vert_s_up_ind)) DEALLOCATE (vert_s_up_ind)
        
        CALL get_aps_and_datas (ap_s_down, data_s_down, &
             lev_window(index_down,1), lev_window(index_down,2), &
             data_present)
           
        nparts_act = SIZE(ap_s_down)
           
        write (*,*) 'number of points on lower level: ',nparts_act
           
        IF (nparts_act < 3) THEN 
           WRITE (*,*) 'Not enough CLaMS data in the layer'
           STOP
        ENDIF
        ALLOCATE (vert_s_down_coor(3,nparts_act))
        ALLOCATE (vert_s_down_nb(nparts_act))
        ALLOCATE (vert_s_down_ind(max_nb,nparts_act))
        CALL set_coor (ap_s_down, vert_s_down_coor)
        CALL triang_qhull (status,vert_s_down_coor,vert_s_down_nb,vert_s_down_ind,&
             triangle_s_down,dim_tr)
         
        firsttime = .FALSE.
        
     ELSE
        NULLIFY (ap_s_down)
        NULLIFY (data_s_down) 
        NULLIFY (vert_s_down_coor)
        NULLIFY (vert_s_down_nb)
        NULLIFY (vert_s_down_ind)
        !nullify (triangle_s_down)
        
        ap_s_down => ap_s_up
        vert_s_down_coor => vert_s_up_coor
        vert_s_down_nb   => vert_s_up_nb
        vert_s_down_ind  => vert_s_up_ind
        data_s_down => data_s_up 
        triangle_s_down => triangle_s_up
           
        nparts_act = SIZE(ap_s_down)
        
        !write (*,*) 'number of points on lower level: ',nparts_act
           
     END IF
        
     t_act = INT((dim_tr-1)/2)
     DO i = 1, npartg_act 
        t_tot_down(i) = find_triangle (t_act, vert_s_down_coor, &
                                       vert_s_down_nb, vert_s_down_ind, &
                                       triangle_s_down, dim_tr,  &
                                       vert_g_coor(:,i), max_dist=max_dist, &
                                       write_warnings=ctrl_out)
        IF (t_tot_down(i)<=0) THEN
           IF (ctrl_out) THEN
              WRITE (*,*) 'Wrong triangle found! Try again!'
              WRITE (*,*)
           ENDIF
           t_act = INT((dim_tr-1)/4)
           t_tot_down(i) = find_triangle (t_act, vert_s_down_coor, &
                                          vert_s_down_nb, vert_s_down_ind, &
                                          triangle_s_down, dim_tr, &
                                          vert_g_coor(:,i), max_dist=max_dist, &
                                          write_warnings=ctrl_out)
           IF (t_tot_down(i)<=0 .AND. ctrl_out) THEN
              WRITE (*,*) 'Wrong triangle found! Set missing value! '
              WRITE (*,*)
           ENDIF
        ENDIF
     END DO
        
     !----------------------------------------------------------------
     ! upper level
     !----------------------------------------------------------------
     NULLIFY (ap_s_up)
     NULLIFY (data_s_up) 
     NULLIFY (vert_s_up_coor)
     NULLIFY (vert_s_up_nb)
     NULLIFY (vert_s_up_ind)
     
     CALL get_aps_and_datas (ap_s_up, data_s_up, &
          lev_window(index_up,1), lev_window(index_up,2), data_present)
        
     nparts_act = SIZE(ap_s_up)
     
     !write (*,*) 'number of points on upper level: ',nparts_act
        
     IF (nparts_act < 3) THEN 
        WRITE (*,*) 'Not enough CLaMS data in the layer'
        STOP
     ENDIF
        
     ALLOCATE (vert_s_up_coor(3,nparts_act))
     ALLOCATE (vert_s_up_nb(nparts_act))
     ALLOCATE (vert_s_up_ind(max_nb,nparts_act))
     CALL set_coor (ap_s_up, vert_s_up_coor)
     
     CALL triang_qhull (status,vert_s_up_coor, vert_s_up_nb, vert_s_up_ind, &
             triangle_s_up, dim_tr)
     
     t_act = INT((dim_tr-1)/2)           
     DO i = 1, npartg_act 
        t_tot_up(i) = find_triangle (t_act, vert_s_up_coor, vert_s_up_nb, &
             vert_s_up_ind, triangle_s_up, dim_tr,  &
             vert_g_coor(:,i), max_dist=max_dist, &
             write_warnings=ctrl_out)
        IF (t_tot_up(i)<=0) THEN
           IF (ctrl_out) THEN
              WRITE (*,*) 'Wrong triangle found! Try again!'
              WRITE (*,*)
           ENDIF
           t_act = INT((dim_tr-1)/4)
           t_tot_up(i) = find_triangle (t_act, vert_s_up_coor, vert_s_up_nb, &
                vert_s_up_ind, triangle_s_up, dim_tr, &
                vert_g_coor(:,i), max_dist=max_dist,  &
                write_warnings=ctrl_out)
           IF (t_tot_up(i)<=0 .AND. ctrl_out) THEN
              WRITE (*,*) 'Wrong triangle found! Set missing value! '
              WRITE (*,*)
           ENDIF
        ENDIF
     END DO
        
     ALLOCATE (mask(npartg_act))
     mask = .TRUE.
        
     js1 = ymds2js(init_date%year,init_date%month,init_date%day,init_date%sec)
     js2 = js1
     ! write (*,*) 'js1=',init_date,js1
     ! write (*,*) 'js2=',date,js2
        
     !----------------------------------------------------------------
     ! Main interpolation
     !----------------------------------------------------------------
     CALL int3d (ap_s_down, vert_s_down_coor, triangle_s_down, t_tot_down, data_s_down, &
          ap_s_up,   vert_s_up_coor,   triangle_s_up,   t_tot_up,   data_s_up, &
          ap_g,      vert_g_coor,      data_g,   npartg_act, mask, &
          delta_lev_down, delta_lev_up, ilev, 2, nlevs, js1, js2)
     
     ! pack ap_g
     ALLOCATE (ap_g_dummy(npartg_act))
     ap_g_dummy = ap_g
     DEALLOCATE (ap_g)
     ALLOCATE (ap_g(COUNT(mask)))
     ap_g = PACK(ap_g_dummy, mask)
     DEALLOCATE (ap_g_dummy)
     
     ! pack data_g
     ALLOCATE (data_g_dummy(ntags,npartg_act))
     data_g_dummy = data_g
     DEALLOCATE (data_g)
     ALLOCATE (data_g(ntags,COUNT(mask)))
     DO itag = 1, ntags
        data_g(itag,:) = PACK(data_g_dummy(itag,:),mask)              
     ENDDO
     DEALLOCATE (data_g_dummy)
     
     ! pack indexarr
     ALLOCATE (indexarr_dummy(npartg_act))
     indexarr_dummy = indexarr
     DEALLOCATE (indexarr)
     ALLOCATE (indexarr(COUNT(mask)))
     indexarr = PACK(indexarr_dummy, mask)
     
     npartg_act = COUNT(mask)
     WRITE (*,*) 'number of interpolated points:',npartg_act
     
     DEALLOCATE (mask)
     DEALLOCATE(indexarr_dummy)
     
     !----------------------------------------------------------------
     ! save data from current level
     !----------------------------------------------------------------
     DO i = 1, npartg_act
        ap_g_tot(indexarr(i)) = ap_g(i)
     ENDDO
     DO itag = 1, ntags
        DO i = 1, npartg_act
           data_g_tot(itag,indexarr(i)) = data_g(itag,i)
        ENDDO
     ENDDO
        
     !----------------------------------------------------------------
     ! deallocate arrays
     !----------------------------------------------------------------
     DEALLOCATE (triangle_s_down)
     DEALLOCATE (indexarr)
     DEALLOCATE (ap_g, data_g, vert_g_coor)
     DEALLOCATE (ap_s_down, data_s_down)
     DEALLOCATE (vert_s_down_coor, vert_s_down_nb, vert_s_down_ind)
     !ap_s_up, data_s_up, vert_s_up: werden fuer den naechsten Durchlauf benoetigt!
     DEALLOCATE (t_tot_down, t_tot_up)
     
  ENDDO  ! loop over all lev levels
  
  !----------------------------------------------------------------
  ! write parameters to ECHAM 
  !----------------------------------------------------------------
  write(*,*)'nlev,nlat,nlon',nlev,nlat,nlon

  CALL nc_write_parameters (nlev,nlat,nlon,ECHAM_DATA)
  
  !---------------------------------------------------------
  ! clean up
  !---------------------------------------------------------
  IF (ASSOCIATED(triangle_s_up)) DEALLOCATE (triangle_s_up)
  IF (ASSOCIATED(ap_s_up)) DEALLOCATE (ap_s_up)
  IF (ASSOCIATED(ap_s_up)) DEALLOCATE (data_s_up)
  IF (ASSOCIATED(vert_s_up_coor)) DEALLOCATE (vert_s_up_coor)
  IF (ASSOCIATED(vert_s_up_nb)) DEALLOCATE (vert_s_up_nb)
  IF (ASSOCIATED(vert_s_up_ind)) DEALLOCATE (vert_s_up_ind)
 
  DEALLOCATE (ap_g_all)
  DEALLOCATE (ap_s_all)
  DEALLOCATE (data_s_all)
  IF (ASSOCIATED(currentday)) DEALLOCATE (currentday)
  
  DEALLOCATE (ap_g_tot, data_g_tot)
  DEALLOCATE (lev_grid)

  END SUBROUTINE rdfrc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE nc_read_input_data (nlev, nlat, nlon, nparts, ap_all, E5_ZETA)

    USE messy_clams_global,      ONLY: E5_LAT, E5_LON
    USE messy_clamsrdfrc_global, ONLY: ap

    IMPLICIT NONE

    TYPE(ap), DIMENSION(:), POINTER :: ap_all
    INTEGER                         :: nlev, nlat, nlon, nparts
    REAL(PREC), DIMENSION(:,:,:), POINTER :: E5_ZETA

    REAL(PREC), DIMENSION(:),     POINTER :: lat, lon
    REAL(PREC), DIMENSION(:,:,:), POINTER :: level, level_temp

    INTEGER :: ilat, ilon, ilev, ipart, i,j,k

    !----------------------------------------------------------------------
    ! read dimensions
    !----------------------------------------------------------------------
    nlat = ny
    nlon = nx
    nlev = SIZE(E5_ZETA,2)

    WRITE (*,*) 'nlat, nlon, nlev: ', nlat, nlon, nlev
    !----------------------------------------------------------------------
    ! read variables
    !----------------------------------------------------------------------
    lat => E5_LAT
    lon => E5_LON

    WRITE (*,*) 'lats:',lat
    WRITE (*,*) 'lons:',lon

    ALLOCATE (level(nlev,nlat,nlon))
!!$  (level_temp(nlon,nlat,nlev))

    level_temp=> E5_ZETA
!!$write(*,*) 'SIZE( E5_ZETA,1)',SIZE( E5_ZETA,1)
!!$write(*,*) 'SIZE( E5_ZETA,2)',SIZE( E5_ZETA,2)
!!$write(*,*) 'SIZE( E5_ZETA,3)',SIZE( E5_ZETA,3)

    DO i=1,nlon
       DO j=1,nlat
          DO k=1,nlev
             level(k,j,i)=level_temp(i,k,j)
          END DO
       END DO
    END DO
    !----------------------------------------------------------------------
    ! write to ap_all
    !----------------------------------------------------------------------
    nparts = nlev*nlat*nlon

    ALLOCATE (ap_all(nparts))

    ipart = 0

    DO ilon = 1, nlon
       DO ilat = 1, nlat
          DO ilev = 1, nlev
             ipart = ipart + 1
             ap_all(ipart)%lat   = lat(ilat)
             ap_all(ipart)%lon   = lon(ilon)
             ap_all(ipart)%lev = level(ilev,ilat,ilon) 
!             ap_all(ipart)%time  = time
          ENDDO
       ENDDO
    ENDDO

    ! clean up
    DEALLOCATE (level)
    NULLIFY (level_temp, lat, lon)

  END SUBROUTINE nc_read_input_data
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Subroutine int3d
!******************************************************************************
! PURPOSE: Main procedure for interpolation. For each g-point the triangle
!          from the layer above and the triangle from the layer below are
!          used. Each triangle contains the g-point (on a unit sphere)
!   Input: ap_s, vert_s, triangle_s t_tot (indicies of triangles which
!          contain g-points) and data for up and down layer
!          ap_g and vert_g for g-points
!  Output: data_g - results sored in the matrix (ntags,npartg)
! Remarks: 
!******************************************************************************

SUBROUTINE int3d (ap_s_down,vert_s_down_coor,triangle_s_down,t_tot_down,data_s_down, &
                  ap_s_up,  vert_s_up_coor,  triangle_s_up,  t_tot_up,  data_s_up, &
                  ap_g,     vert_g_coor,     data_g, npartg_act,  mask, &
                  delta_lev_down, delta_lev_up, ilev, wt, nlevs, js1, js2)

  USE messy_clamsrdfrc_global
  USE messy_clams_tools_interpol, ONLY: determ_weight, vert_interpol

  IMPLICIT NONE

  TYPE(ap),       DIMENSION(:),   POINTER :: ap_s_down, ap_s_up, ap_g
  REAL(PREC),           DIMENSION(:,:), POINTER :: vert_s_down_coor, vert_s_up_coor, vert_g_coor
  INTEGER,        DIMENSION(:,:), POINTER :: triangle_s_down, triangle_s_up
  REAL(PREC),           DIMENSION(:,:), POINTER :: data_s_down, data_s_up, data_g
  INTEGER,        DIMENSION(:),   POINTER :: t_tot_down, t_tot_up
  LOGICAL,        DIMENSION(:),   POINTER :: mask
  REAL(PREC)                                    :: delta_lev_down, delta_lev_up
  INTEGER                                 :: npartg_act
  INTEGER                                 :: ilev, wt, nlevs
  REAL(KIND(0d0))                         :: js1, js2
 
  REAL(PREC)                    :: lev_down, lev_up, value_down, value_up
  REAL(PREC)                   :: vert_up, vert_down, ratio
  REAL(PREC),    DIMENSION(3,3):: EckCoor
  REAL(PREC),    DIMENSION(3)  :: weight_up, weight_down
  INTEGER                :: itag,ipart,i
  INTEGER, DIMENSION(3)  :: ind_up, ind_down

  REAL(PREC) , DIMENSION(3) :: lev_arr

  REAL(PREC) , PARAMETER :: eps=1E-6


      
  DO ipart = 1, npartg_act

     !----------------------------------------------------------------
     ! set missing value if triangle not found
     !----------------------------------------------------------------
     IF (t_tot_down(ipart)<0 .OR. t_tot_up(ipart)<0) THEN
        DO itag = 1, ntags
           data_g(itag,ipart)=mdi
        ENDDO
        CYCLE
     ENDIF
        
     !----------------------------------------------------------------
     ! get weights on lower lev level
     !----------------------------------------------------------------
     ind_down = triangle_s_down(:,t_tot_down(ipart))

     DO i=1,3
        EckCoor(:,i) = vert_s_down_coor(:,ind_down(i))
     ENDDO
     
     IF (ap_g(ipart)%lat<lat_down .OR. ap_g(ipart)%lat>lat_up) THEN
        ratio = r_coarse/earth_radius
     ELSE
        ratio = r_high/earth_radius
     END IF

     ! Umspeichern auf lev_arr fuehrt zu drastischer Laufzeitverbesserung!
     lev_arr = ap_s_down(ind_down)%lev
     CALL determ_weight (EckCoor,lev_arr,&
          vert_g_coor(:,ipart),ap_g(ipart)%lev,&
          weight_down,lev_down,delta_lev_down,ratio,wt)
     
     !----------------------------------------------------------------
     ! get weights on upper lev level
     !----------------------------------------------------------------
     ind_up = triangle_s_up(:,t_tot_up(ipart))
     
     DO i=1,3
        EckCoor(:,i) = vert_s_up_coor(:,ind_up(i))
     ENDDO
     
     
     ! Umspeichern auf lev_arr fuehrt zu drastischer Laufzeitverbesserung!
     lev_arr = ap_s_up(ind_up)%lev
     CALL determ_weight (EckCoor,lev_arr,vert_g_coor(:,ipart),&
          ap_g(ipart)%lev, weight_up,lev_up,delta_lev_up,ratio,wt)
     
!!!!!
     IF (nlevs > 1) THEN
        IF ((ap_g(ipart)%lev<lev_down .AND. ilev/=1) .OR. &
             (ap_g(ipart)%lev>lev_up .AND. ilev/=nlevs+1)) THEN
           mask(ipart) = .FALSE.
           CYCLE
        ENDIF
     ENDIF

     DO itag = 1, ntags

!!$        IF (tag_names(itag)=='DELTA_GPH') THEN
!!$           
!!$           value_down = (data_s_down(itag,ind_down(1))+data_s_down(itag,ind_down(2)) &
!!$                        +data_s_down(itag,ind_down(3))) / 3.
!!$
!!$           value_up = (data_s_up(itag,ind_up(1))+data_s_up(itag,ind_up(2)) &
!!$                       +data_s_up(itag,ind_up(3))) / 3.
!!$
!!$           data_g(itag,ipart) = value_up - value_down
!!$
!!$        ELSE
!!$
!!$           !----------------------------------------------------------------
           ! get value on lower lev level
           !----------------------------------------------------------------
           ! wenn eines der Gewichte 1. ist:
           IF (COUNT(ABS(1.-weight_down(:))<eps) > 0) THEN
              IF (ABS(weight_down(1)-1.)<eps) THEN
                 value_down = data_s_down(itag,ind_down(1))
              ELSEIF (ABS(weight_down(2)-1.)<eps) THEN
                 value_down = data_s_down(itag,ind_down(2))
              ELSE
                 value_down = data_s_down(itag,ind_down(3))
              ENDIF
              
           ELSE
              
              ! Ueberpruefe, ob einer der Datenwerte MDI ist
              SELECT CASE (COUNT(ABS(mdi-data_s_down(itag,ind_down(:)))<eps))
              CASE(0) ! kein missing_value vorhanden:
                 value_down = SUM(weight_down(:)*data_s_down(itag,ind_down(:)))
              CASE(1) ! ein missing_value: Gewicht wird auf die anderen beiden verteilt
                 IF (ABS(mdi-data_s_down(itag,ind_down(1)))<eps) THEN
                    value_down = data_s_down(itag,ind_down(2))* &
                              (weight_down(2)+weight_down(1)*weight_down(2)/(weight_down(2)+weight_down(3)))+ &
                           data_s_down(itag,ind_down(3))* &
                              (weight_down(3)+weight_down(1)*weight_down(3)/(weight_down(2)+weight_down(3)))
                 ELSEIF (ABS(mdi-data_s_down(itag,ind_down(2)))<eps) THEN
                    value_down = data_s_down(itag,ind_down(1))* &
                         (weight_down(1)+weight_down(2)*weight_down(1)/(weight_down(1)+weight_down(3)))+ &
                         data_s_down(itag,ind_down(3))* &
                              (weight_down(3)+weight_down(2)*weight_down(3)/(weight_down(1)+weight_down(3)))
                 ELSE
                    value_down = data_s_down(itag,ind_down(1))* &
                         (weight_down(1)+weight_down(3)*weight_down(1)/(weight_down(1)+weight_down(2)))+ &
                           data_s_down(itag,ind_down(2))* &
                              (weight_down(2)+weight_down(3)*weight_down(2)/(weight_down(1)+weight_down(2)))
                 ENDIF

              CASE DEFAULT ! mehr als ein missing_value:
                 value_down = mdi
              END SELECT
              
           ENDIF
           
           !----------------------------------------------------------------
           ! get value on upper lev level
           !----------------------------------------------------------------
           ! wenn eines der Gewichte 1. ist:
           IF (COUNT(ABS(1.-weight_up(:))<eps) > 0) THEN
              IF (ABS(weight_up(1)-1.)<eps) THEN
                 value_up = data_s_up(itag,ind_up(1))
              ELSEIF (ABS(weight_up(2)-1.)<eps) THEN
                 value_up = data_s_up(itag,ind_up(2))
              ELSE
                 value_up = data_s_up(itag,ind_up(3))
              ENDIF
              
           ELSE
              
              ! Ueberpruefe, ob einer der Datenwerte MDI ist
              SELECT CASE (COUNT(ABS(mdi-data_s_up(itag,ind_up(:)))<eps))
              CASE(0) ! kein missing_value vorhanden:
                 value_up = SUM(weight_up(:)*data_s_up(itag,ind_up(:)))
              CASE(1) ! ein missing_value: Gewicht wird auf die anderen beiden verteilt
                 IF (ABS(mdi-data_s_up(itag,ind_up(1)))<eps) THEN
                    value_up = data_s_up(itag,ind_up(2))* &
                         (weight_up(2)+weight_up(1)*weight_up(2)/(weight_up(2)+weight_up(3)))+ &
                         data_s_up(itag,ind_up(3))* &
                         (weight_up(3)+weight_up(1)*weight_up(3)/(weight_up(2)+weight_up(3)))
                 ELSEIF (ABS(mdi-data_s_up(itag,ind_up(2)))<eps) THEN
                    value_up = data_s_up(itag,ind_up(1))* &
                         (weight_up(1)+weight_up(2)*weight_up(1)/(weight_up(1)+weight_up(3)))+ &
                         data_s_up(itag,ind_up(3))* &
                         (weight_up(3)+weight_up(2)*weight_up(3)/(weight_up(1)+weight_up(3)))
                 ELSE
                    value_up = data_s_up(itag,ind_up(1))* &
                         (weight_up(1)+weight_up(3)*weight_up(1)/(weight_up(1)+weight_up(2)))+ &
                         data_s_up(itag,ind_up(2))* &
                         (weight_up(2)+weight_up(3)*weight_up(2)/(weight_up(1)+weight_up(2)))
                 ENDIF
              CASE default ! mehr als ein missing_value:
                 value_up = mdi
              END SELECT
              
!!$           ENDIF
           
           !----------------------------------------------------------------
           ! vertical interpolation
           !----------------------------------------------------------------
           IF (ABS(mdi-value_down)>eps .AND. ABS(mdi-value_up)>eps) THEN
              IF (wt==0) THEN
                 CALL vert_interpol (lev_down, lev_up, ap_g(ipart)%lev, &
                      value_down, value_up, data_g(itag,ipart), vert_up, vert_down, .FALSE.)
              ELSE
                 CALL vert_interpol (lev_down, lev_up, ap_g(ipart)%lev, &
                      value_down, value_up, data_g(itag,ipart), vert_up, vert_down, .TRUE.)
              ENDIF
           ELSEIF (ABS(mdi-value_down)>eps) THEN
              data_g(itag,ipart)=value_down
           ELSEIF (ABS(mdi-value_up)>eps) THEN
              data_g(itag,ipart)=value_up
           ELSE
              data_g(itag,ipart)=mdi
           ENDIF
       
        ENDIF

     END DO

  END DO
  
  
END SUBROUTINE int3d

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamsrdfrc_read_nml(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  
    IMPLICIT NONE

    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamsrdfrc_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    
    NAMELIST /CTRL/ use_o3forc, use_h2oforc, timestep_rdfrc, bound_trop_n,     &
         &    bound_trop_s, lev_bound_trop, lev_bound_extr, bound_up_zeta

    status = 1 !ERROR
    
    ! Read namelist variables:
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 !NO ERROR
 
  END SUBROUTINE clamsrdfrc_read_nml
!-----------------------------------------------------------------------


!***********************************************************************!
END MODULE messy_clamsrdfrc
