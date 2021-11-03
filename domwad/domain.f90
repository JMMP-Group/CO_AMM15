MODULE domain
   !!==============================================================================
   !!                       ***  MODULE domain   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================
   !! History :  OPA  !  1990-10  (C. Levy - G. Madec)  Original code
   !!                 !  1992-01  (M. Imbard) insert time step initialization
   !!                 !  1996-06  (G. Madec) generalized vertical coordinate 
   !!                 !  1997-02  (G. Madec) creation of domwri.F
   !!                 !  2001-05  (E.Durand - G. Madec) insert closed sea
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.3  !  2010-11  (G. Madec)  initialisation in C1D configuration
   !!            3.6  !  2013     ( J. Simeon, C. Calone, G. Madec, C. Ethe ) Online coarsening of outputs
   !!            3.7  !  2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dom_init       : initialize the space and time domain
   !!   dom_nam        : read and contral domain namelists
   !!   dom_ctl        : control print for the ocean domain
   !!----------------------------------------------------------------------
   USE oce             ! ocean variables
   USE dom_oce         ! domain: ocean
   USE phycst          ! physical constants
   USE closea          ! closed seas
   USE domhgr          ! domain: set the horizontal mesh
   USE domzgr          ! domain: set the vertical mesh
   USE domstp          ! domain: set the time-step
   USE dommsk          ! domain: set the mask system
   USE domwri          ! domain: write the meshmask file
   USE domvvl          ! variable volume
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! 
   USE wrk_nemo        ! Memory Allocation
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary condition (or mpp link)
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_init   ! called by opa.F90

   !!-------------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: domain.F90 6140 2015-12-21 11:35:23Z timgraham $
   !! Software governed by the CeCILL licence        (./LICENSE)
   !!-------------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_init  ***
      !!                    
      !! ** Purpose :   Domain initialization. Call the routines that are 
      !!              required to create the arrays which define the space 
      !!              and time domain of the ocean model.
      !!
      !! ** Method  : - dom_msk: compute the masks from the bathymetry file
      !!              - dom_hgr: compute or read the horizontal grid-point position
      !!                         and scale factors, and the coriolis factor
      !!              - dom_zgr: define the vertical coordinate and the bathymetry
      !!              - dom_stp: defined the model time step
      !!              - dom_wri: create the meshmask file if nmsh=1
      !!              - 1D configuration, move Coriolis, u and v at T-point
      !!----------------------------------------------------------------------
      INTEGER ::   jk          ! dummy loop indices
      INTEGER ::   iconf = 0   ! local integers
      REAL(wp), POINTER, DIMENSION(:,:) ::   z1_hu_0, z1_hv_0
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dom_init')
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_init : domain initialization'
         WRITE(numout,*) '~~~~~~~~'
      ENDIF
      !
      !                       !==  Reference coordinate system  ==!
      !
                     CALL dom_nam               ! read namelist ( namrun, namdom )
                     CALL dom_clo               ! Closed seas and lake
                     CALL dom_hgr               ! Horizontal mesh
                     CALL dom_zgr               ! Vertical mesh and bathymetry
                     CALL dom_msk               ! Masks
      !
      ht_0(:,:) = e3t_0(:,:,1) * tmask(:,:,1)   ! Reference ocean thickness
      hu_0(:,:) = e3u_0(:,:,1) * umask(:,:,1)
      hv_0(:,:) = e3v_0(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpk
         ht_0(:,:) = ht_0(:,:) + e3t_0(:,:,jk) * tmask(:,:,jk)
         hu_0(:,:) = hu_0(:,:) + e3u_0(:,:,jk) * umask(:,:,jk)
         hv_0(:,:) = hv_0(:,:) + e3v_0(:,:,jk) * vmask(:,:,jk)
      END DO
      !
      !              !==  time varying part of coordinate system  ==!
      !
      IF( ln_linssh ) THEN          ! Fix in time : set to the reference one for all
         !       before        !          now          !       after         !
         ;  gdept_b = gdept_0  ;   gdept_n = gdept_0   !        ---          ! depth of grid-points
         ;  gdepw_b = gdepw_0  ;   gdepw_n = gdepw_0   !        ---          !
         ;                     ;   gde3w_n = gde3w_0   !        ---          !
         !                                                                  
         ;    e3t_b =   e3t_0  ;     e3t_n =   e3t_0   ;   e3t_a =  e3t_0    ! scale factors
         ;    e3u_b =   e3u_0  ;     e3u_n =   e3u_0   ;   e3u_a =  e3u_0    !
         ;    e3v_b =   e3v_0  ;     e3v_n =   e3v_0   ;   e3v_a =  e3v_0    !
         ;                     ;     e3f_n =   e3f_0   !        ---          !
         ;    e3w_b =   e3w_0  ;     e3w_n =   e3w_0   !        ---          !
         ;   e3uw_b =  e3uw_0  ;    e3uw_n =  e3uw_0   !        ---          !
         ;   e3vw_b =  e3vw_0  ;    e3vw_n =  e3vw_0   !        ---          !
         !
         CALL wrk_alloc( jpi,jpj,   z1_hu_0, z1_hv_0 )
         !
         z1_hu_0(:,:) = ssumask(:,:) / ( hu_0(:,:) + 1._wp - ssumask(:,:) )     ! _i mask due to ISF
         z1_hv_0(:,:) = ssvmask(:,:) / ( hv_0(:,:) + 1._wp - ssvmask(:,:) )
         !
         !        before       !          now          !       after         !
         ;                     ;      ht_n =    ht_0   !                     ! water column thickness
         ;     hu_b =    hu_0  ;      hu_n =    hu_0   ;    hu_a =    hu_0   ! 
         ;     hv_b =    hv_0  ;      hv_n =    hv_0   ;    hv_a =    hv_0   !
         ;  r1_hu_b = z1_hu_0  ;   r1_hu_n = z1_hu_0   ; r1_hu_a = z1_hu_0   ! inverse of water column thickness
         ;  r1_hv_b = z1_hv_0  ;   r1_hv_n = z1_hv_0   ; r1_hv_a = z1_hv_0   !
         !
         CALL wrk_dealloc( jpi,jpj,   z1_hu_0, z1_hv_0 )
         !
      ELSE                         ! time varying : initialize before/now/after variables
         !
         CALL dom_vvl_init 
         !
      ENDIF
      !
      IF( nmsh /= 0 ) CALL dom_wri ! DB create mesh_mask
      CALL cfg_write               ! create the configuration file
      !
      IF( nn_timing == 1 )   CALL timing_stop('dom_init')
      !
   END SUBROUTINE dom_init


   SUBROUTINE dom_nam
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read domaine namelists and print the variables.
      !!
      !! ** input   : - namrun namelist
      !!              - namdom namelist
      !!              - namnc4 namelist   ! "key_netcdf4" only
      !!----------------------------------------------------------------------
      USE ioipsl
      NAMELIST/namrun/ cn_ocerst_indir, cn_ocerst_outdir, nn_stocklist, ln_rst_list,                 &
                       nn_no   , cn_exp   , cn_ocerst_in, cn_ocerst_out, ln_rstart , nn_rstctl ,     &
         &             nn_it000, nn_itend , nn_date0    , nn_time0     , nn_leapy  , nn_istate ,     &
         &             nn_stock, nn_write , ln_mskland  , ln_clobber   , nn_chunksz, nn_euler  ,     &
         &             ln_cfmeta, ln_iscpl
      NAMELIST/namdom/ nn_bathy, rn_bathy , rn_e3zps_min, rn_e3zps_rat, nn_msh, rn_hmin, rn_isfhmin, &
         &             rn_atfp , rn_rdt   , nn_closea   , ln_crs      , jphgr_msh ,                  &
         &             ppglam0, ppgphi0, ppe1_deg, ppe2_deg, ppe1_m, ppe2_m,                         &
         &             ppsur, ppa0, ppa1, ppkth, ppacr, ppdzmin, pphmax, ldbletanh,                  &
         &             ppa2, ppkth2, ppacr2, rn_wd_ref_depth



      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namrun in reference namelist : Parameters of the run
      READ  ( numnam_ref, namrun, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namrun in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namrun in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namrun, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namrun in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namrun )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_nam  : domain initialization through namelist read'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namrun'
         WRITE(numout,*) '      job number                      nn_no      = ', nn_no
         WRITE(numout,*) '      experiment name for output      cn_exp     = ', cn_exp
         WRITE(numout,*) '      file prefix restart input       cn_ocerst_in= ', cn_ocerst_in
         WRITE(numout,*) '      restart input directory         cn_ocerst_indir= ', cn_ocerst_indir
         WRITE(numout,*) '      file prefix restart output      cn_ocerst_out= ', cn_ocerst_out
         WRITE(numout,*) '      restart output directory        cn_ocerst_outdir= ', cn_ocerst_outdir
         WRITE(numout,*) '      restart logical                 ln_rstart  = ', ln_rstart
         WRITE(numout,*) '      start with forward time step    nn_euler   = ', nn_euler
         WRITE(numout,*) '      control of time step            nn_rstctl  = ', nn_rstctl
         WRITE(numout,*) '      number of the first time step   nn_it000   = ', nn_it000
         WRITE(numout,*) '      number of the last time step    nn_itend   = ', nn_itend
         WRITE(numout,*) '      initial calendar date aammjj    nn_date0   = ', nn_date0
         WRITE(numout,*) '      initial time of day in hhmm     nn_time0   = ', nn_time0
         WRITE(numout,*) '      leap year calendar (0/1)        nn_leapy   = ', nn_leapy
         WRITE(numout,*) '      initial state output            nn_istate  = ', nn_istate
         IF( ln_rst_list ) THEN
            WRITE(numout,*) '      list of restart dump times      nn_stocklist   =', nn_stocklist
         ELSE
            WRITE(numout,*) '      frequency of restart file       nn_stock   = ', nn_stock
         ENDIF
         WRITE(numout,*) '      frequency of output file        nn_write   = ', nn_write
         WRITE(numout,*) '      mask land points                ln_mskland = ', ln_mskland
         WRITE(numout,*) '      additional CF standard metadata ln_cfmeta  = ', ln_cfmeta
         WRITE(numout,*) '      overwrite an existing file      ln_clobber = ', ln_clobber
         WRITE(numout,*) '      NetCDF chunksize (bytes)        nn_chunksz = ', nn_chunksz
         WRITE(numout,*) '      IS coupling at the restart step ln_iscpl   = ', ln_iscpl
         WRITE(numout,*) '      rn_wd_ref_depth =', rn_wd_ref_depth
      ENDIF

      no = nn_no                    ! conversion DOCTOR names into model names (this should disappear soon)
      cexper = cn_exp
      nrstdt = nn_rstctl
      nit000 = nn_it000
      nitend = nn_itend
      ndate0 = nn_date0
      nleapy = nn_leapy
      ninist = nn_istate
      nstock = nn_stock
      nstocklist = nn_stocklist
      nwrite = nn_write
      neuler = nn_euler
      IF ( neuler == 1 .AND. .NOT. ln_rstart ) THEN
         WRITE(ctmp1,*) 'ln_rstart =.FALSE., nn_euler is forced to 0 '
         CALL ctl_warn( ctmp1 )
         neuler = 0
      ENDIF

      !                             ! control of output frequency
      IF ( nstock == 0 .OR. nstock > nitend ) THEN
         WRITE(ctmp1,*) 'nstock = ', nstock, ' it is forced to ', nitend
         CALL ctl_warn( ctmp1 )
         nstock = nitend
      ENDIF
      IF ( nwrite == 0 ) THEN
         WRITE(ctmp1,*) 'nwrite = ', nwrite, ' it is forced to ', nitend
         CALL ctl_warn( ctmp1 )
         nwrite = nitend
      ENDIF




      SELECT CASE ( nleapy )        ! Choose calendar for IOIPSL
      CASE (  1 ) 
         CALL ioconf_calendar('gregorian')
         IF(lwp) WRITE(numout,*) '   The IOIPSL calendar is "gregorian", i.e. leap year'
      CASE (  0 )
         CALL ioconf_calendar('noleap')
         IF(lwp) WRITE(numout,*) '   The IOIPSL calendar is "noleap", i.e. no leap year'
      CASE ( 30 )
         CALL ioconf_calendar('360d')
         IF(lwp) WRITE(numout,*) '   The IOIPSL calendar is "360d", i.e. 360 days in a year'
      END SELECT




      REWIND( numnam_ref )              ! Namelist namdom in reference namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_ref, namdom, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdom in reference namelist', lwp )
  
      !
      REWIND( numnam_cfg )              ! Namelist namdom in configuration namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_cfg, namdom, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdom in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namdom )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namdom : space & time domain'
         WRITE(numout,*) '      flag read/compute bathymetry      nn_bathy     = ', nn_bathy
         WRITE(numout,*) '      Depth (if =0 bathy=jpkm1)         rn_bathy     = ', rn_bathy
         WRITE(numout,*) '      min depth of the ocean    (>0) or    rn_hmin   = ', rn_hmin
         WRITE(numout,*) '      min number of ocean level (<0)       '
         WRITE(numout,*) '      treshold to open the isf cavity   rn_isfhmin   = ', rn_isfhmin, ' (m)'
         WRITE(numout,*) '      minimum thickness of partial      rn_e3zps_min = ', rn_e3zps_min, ' (m)'
         WRITE(numout,*) '         step level                     rn_e3zps_rat = ', rn_e3zps_rat
         WRITE(numout,*) '      create mesh/mask file(s)          nn_msh       = ', nn_msh
         WRITE(numout,*) '           = 0   no file created           '
         WRITE(numout,*) '           = 1   mesh_mask                 '
         WRITE(numout,*) '           = 2   mesh and mask             '
         WRITE(numout,*) '           = 3   mesh_hgr, msh_zgr and mask'
         WRITE(numout,*) '      ocean time step                       rn_rdt    = ', rn_rdt
         WRITE(numout,*) '      asselin time filter parameter         rn_atfp   = ', rn_atfp
         WRITE(numout,*) '      suppression of closed seas (=0)       nn_closea = ', nn_closea
         WRITE(numout,*) '      online coarsening of dynamical fields ln_crs    = ', ln_crs
         WRITE(numout,*) '      type of horizontal mesh jphgr_msh           = ', jphgr_msh
         WRITE(numout,*) '      longitude of first raw and column T-point ppglam0 = ', ppglam0
         WRITE(numout,*) '      latitude  of first raw and column T-point ppgphi0 = ', ppgphi0
         WRITE(numout,*) '      zonal      grid-spacing (degrees) ppe1_deg        = ', ppe1_deg
         WRITE(numout,*) '      meridional grid-spacing (degrees) ppe2_deg        = ', ppe2_deg
         WRITE(numout,*) '      zonal      grid-spacing (degrees) ppe1_m          = ', ppe1_m
         WRITE(numout,*) '      meridional grid-spacing (degrees) ppe2_m          = ', ppe2_m
         WRITE(numout,*) '      ORCA r4, r2 and r05 coefficients  ppsur           = ', ppsur
         WRITE(numout,*) '                                        ppa0            = ', ppa0
         WRITE(numout,*) '                                        ppa1            = ', ppa1
         WRITE(numout,*) '                                        ppkth           = ', ppkth
         WRITE(numout,*) '                                        ppacr           = ', ppacr
         WRITE(numout,*) '      Minimum vertical spacing ppdzmin                  = ', ppdzmin
         WRITE(numout,*) '      Maximum depth pphmax                              = ', pphmax
         WRITE(numout,*) '      Use double tanf function for vertical coordinates ldbletanh = ', ldbletanh
         WRITE(numout,*) '      Double tanh function parameters ppa2              = ', ppa2
         WRITE(numout,*) '                                      ppkth2            = ', ppkth2
         WRITE(numout,*) '                                      ppacr2            = ', ppacr2
      ENDIF
      !
      ntopo     = nn_bathy          ! conversion DOCTOR names into model names (this should disappear soon)
      e3zps_min = rn_e3zps_min
      e3zps_rat = rn_e3zps_rat
      nmsh      = nn_msh
      atfp      = rn_atfp
      rdt       = rn_rdt

      snc4set%luse = .FALSE.        ! No NetCDF 4 case
      !
   END SUBROUTINE dom_nam


   SUBROUTINE dom_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_ctl  ***
      !!
      !! ** Purpose :   Domain control.
      !!
      !! ** Method  :   compute and print extrema of masked scale factors
      !!----------------------------------------------------------------------
      INTEGER ::   iimi1, ijmi1, iimi2, ijmi2, iima1, ijma1, iima2, ijma2
      INTEGER, DIMENSION(2) ::   iloc   ! 
      REAL(wp) ::   ze1min, ze1max, ze2min, ze2max
      !!----------------------------------------------------------------------
      !
      IF(lk_mpp) THEN
         CALL mpp_minloc( e1t(:,:), tmask_i(:,:), ze1min, iimi1,ijmi1 )
         CALL mpp_minloc( e2t(:,:), tmask_i(:,:), ze2min, iimi2,ijmi2 )
         CALL mpp_maxloc( e1t(:,:), tmask_i(:,:), ze1max, iima1,ijma1 )
         CALL mpp_maxloc( e2t(:,:), tmask_i(:,:), ze2max, iima2,ijma2 )
      ELSE
         ze1min = MINVAL( e1t(:,:), mask = tmask_i(:,:) == 1._wp )    
         ze2min = MINVAL( e2t(:,:), mask = tmask_i(:,:) == 1._wp )    
         ze1max = MAXVAL( e1t(:,:), mask = tmask_i(:,:) == 1._wp )    
         ze2max = MAXVAL( e2t(:,:), mask = tmask_i(:,:) == 1._wp )    

         iloc  = MINLOC( e1t(:,:), mask = tmask_i(:,:) == 1._wp )
         iimi1 = iloc(1) + nimpp - 1
         ijmi1 = iloc(2) + njmpp - 1
         iloc  = MINLOC( e2t(:,:), mask = tmask_i(:,:) == 1._wp )
         iimi2 = iloc(1) + nimpp - 1
         ijmi2 = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e1t(:,:), mask = tmask_i(:,:) == 1._wp )
         iima1 = iloc(1) + nimpp - 1
         ijma1 = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e2t(:,:), mask = tmask_i(:,:) == 1._wp )
         iima2 = iloc(1) + nimpp - 1
         ijma2 = iloc(2) + njmpp - 1
      ENDIF
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_ctl : extrema of the masked scale factors'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,"(14x,'e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, iima1, ijma1
         WRITE(numout,"(14x,'e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, iimi1, ijmi1
         WRITE(numout,"(14x,'e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, iima2, ijma2
         WRITE(numout,"(14x,'e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, iimi2, ijmi2
      ENDIF
      !
   END SUBROUTINE dom_ctl


   SUBROUTINE cfg_write
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cfg_write  ***
      !!                   
      !! ** Purpose :   Create the "domain_cfg" file, a NetCDF file which 
      !!              contains all the ocean domain informations required to 
      !!              define an ocean configuration.
      !!
      !! ** Method  :   Write in a file all the arrays required to set up an
      !!              ocean configuration.
      !!
      !! ** output file :   domain_cfg.nc : domain size, characteristics,horizontal mesh,
      !!                              Coriolis parameter, and vertical scale factors
      !!                              NB: also contains ORCA family information (if cp_cfg = "ORCA")
      !!                              and depths (ln_e3_dep=F) 
      !!----------------------------------------------------------------------
      INTEGER                          ::   ji, jj, jk             ! dummy loop indices
      INTEGER                          ::   izco, izps, isco, icav
      INTEGER                          ::   inum                   ! temporary units for 
                                                                   ! 'domain_cfg.nc' file
      CHARACTER(len=21)                ::   clnam                  ! filename (mesh and 
                                                                   ! mask informations)
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d                    ! workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z3d                    ! workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cfg_write : create the "domain_cfg.nc" file containing all required configuration information'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'
      !
      !                       ! ============================= !
      !                       !  create 'domain_cfg.nc' file  !
      !                       ! ============================= !
      !         
      clnam = 'domain_cfg'  ! filename (configuration information)
      CALL iom_open( TRIM(clnam), inum, ldwrt = .TRUE., kiolib = jprstlib )
      
      !
      !                             !==  ORCA family specificities  ==!
      IF( cp_cfg == "ORCA" ) THEN
         CALL iom_rstput( 0, 0, inum, 'ORCA'      , 1._wp            , ktype = jp_i4 )
         CALL iom_rstput( 0, 0, inum, 'ORCA_index', REAL( jp_cfg, wp), ktype = jp_i4 )         
      ENDIF
      !                             !==  global domain size  ==!
      !
      CALL iom_rstput( 0, 0, inum, 'jpiglo', REAL( jpiglo, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'jpjglo', REAL( jpjglo, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'jpkglo', REAL( jpk   , wp), ktype = jp_i4 )
      !
      !                             !==  domain characteristics  ==!
      !
      !                                   ! lateral boundary of the global
      !                                   domain
      CALL iom_rstput( 0, 0, inum, 'jperio', REAL( jperio, wp), ktype = jp_i4 )
      !
      !                                   ! type of vertical coordinate
      IF( ln_zco             ) THEN   ;   izco = 1   ;   ELSE   ;   izco = 0   ;   ENDIF
      IF( ln_zps             ) THEN   ;   izps = 1   ;   ELSE   ;   izps = 0   ;   ENDIF
      IF( ln_sco .OR. ln_mes ) THEN   ;   isco = 1   ;   ELSE   ;   isco = 0   ;   ENDIF
      CALL iom_rstput( 0, 0, inum, 'ln_zco'   , REAL( izco, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'ln_zps'   , REAL( izps, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'ln_sco'   , REAL( isco, wp), ktype = jp_i4 )
      !
      !                                   ! ocean cavities under iceshelves
      IF( ln_isfcav ) THEN   ;   icav = 1   ;   ELSE   ;   icav = 0   ;   ENDIF
      CALL iom_rstput( 0, 0, inum, 'ln_isfcav', REAL( icav, wp), ktype = jp_i4 )
      !
      !                             !==  horizontal mesh  !
      !

      !CEOD Just force it to be this for now
      z2d(:,:) = hbatt(:,:) ! add back on reference height to get appox dep
                                 !this is later corrected for with specified min depth bg user for above greoid
                                 ! WAD points
      !where (z2d   (:,:).lte.1e-5)  z2d(:,:) = -10.0
      where (tmask  (:,:,1).eq.0)  z2d(:,:) = 0.0
      !CEOD CALL iom_rstput( 0, 0, inum, 'rn_wd_ref_depth'   , rn_wd_ref_depth, ktype = jp_i4  ) ! replace this later with variable
      CALL iom_rstput( 0, 0, inum, 'rn_wd_ref_depth'   , rn_wd_ref_depth  ) ! replace this later with variable
      CALL iom_rstput( 0, 0, inum, 'ht_wd', z2d )        !    ht_wd
      !CEOD

      CALL iom_rstput( 0, 0, inum, 'glamt', glamt, ktype = jp_r8 )   ! latitude
      CALL iom_rstput( 0, 0, inum, 'glamu', glamu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamv', glamv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamf', glamf, ktype = jp_r8 )
      !                                
      CALL iom_rstput( 0, 0, inum, 'gphit', gphit, ktype = jp_r8 )   ! longitude
      CALL iom_rstput( 0, 0, inum, 'gphiu', gphiu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphiv', gphiv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphif', gphif, ktype = jp_r8 )
      !                                
      CALL iom_rstput( 0, 0, inum, 'e1t'  , e1t  , ktype = jp_r8 )   ! i-scale factors (e1.)
      CALL iom_rstput( 0, 0, inum, 'e1u'  , e1u  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1v'  , e1v  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1f'  , e1f  , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'e2t'  , e2t  , ktype = jp_r8 )   ! j-scale factors (e2.)
      CALL iom_rstput( 0, 0, inum, 'e2u'  , e2u  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2v'  , e2v  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2f'  , e2f  , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'ff_f' , ff_f , ktype = jp_r8 )   ! coriolis factor
      CALL iom_rstput( 0, 0, inum, 'ff_t' , ff_t , ktype = jp_r8 )
      !
      !                             !==  vertical mesh  ==!
      !                                                     
      CALL iom_rstput( 0, 0, inum, 'e3t_1d'  , e3t_1d  , ktype = jp_r8 )   !  reference 1D-coordinate
      CALL iom_rstput( 0, 0, inum, 'e3w_1d'  , e3w_1d  , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'e3t_0'   , e3t_0   , ktype = jp_r8 )   !  vertical scale factors (e
      CALL iom_rstput( 0, 0, inum, 'e3u_0'   , e3u_0   , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3v_0'   , e3v_0   , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3f_0'   , e3f_0   , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3w_0'   , e3w_0   , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3uw_0'  , e3uw_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3vw_0'  , e3vw_0  , ktype = jp_r8 )
      !
      IF(.NOT.ln_e3_dep ) THEN                                             !  depth (t- & w-points)
         CALL iom_rstput( 0, 0, inum, 'gdept_1d', gdept_1d, ktype = jp_r8 )   ! required only with  
         CALL iom_rstput( 0, 0, inum, 'gdepw_1d', gdepw_1d, ktype = jp_r8 )   ! the old e3. definition
         CALL iom_rstput( 0, 0, inum, 'gdept_0' , gdept_0 , ktype = jp_r8 )
         CALL iom_rstput( 0, 0, inum, 'gdepw_0' , gdepw_0 , ktype = jp_r8 )
      ENDIF
      !                                         
      !                             !==  ocean top and bottom level  ==!
      !
      CALL iom_rstput( 0, 0, inum, 'bottom_level' , REAL( mbkt, wp )*ssmask , ktype = jp_i4 )   ! nb of ocean T-points
      CALL iom_rstput( 0, 0, inum, 'top_level'    , REAL( mikt, wp )*ssmask , ktype = jp_i4 )   ! nb of ocean T-points (ISF)
      DO jj = 1,jpj
         DO ji = 1,jpi
            z2d (ji,jj) = SUM ( e3t_0(ji,jj, 1:mbkt(ji,jj) ) ) * ssmask(ji,jj) 
         END DO
      END DO
      CALL iom_rstput( 0, 0, inum, 'bathy_meter'   , z2d , ktype = jp_r4 )

      !
      IF( ln_sco .OR. ln_mes ) THEN             ! s-coordinate: store grid stiffness ratio  (Not required anyway)
         CALL dom_stiff( z2d )
         CALL iom_rstput( 0, 0, inum, 'stiffness', z2d )        ! Max. grid stiffness ratio
         CALL dom_stiff_3D( z3d )
         CALL iom_rstput( 0, 0, inum, 'stiff3D', z3d )          ! 3D stiffness ratio
         CALL saw_tooth( z2d )
         CALL iom_rstput( 0, 0, inum, 'saw_tooth', z2d )        ! Saw_tooth diagnostic
      ENDIF
      !
      !                                ! ============================
      !                                !        close the files 
      !                                ! ============================
      CALL iom_close( inum )
      !
   END SUBROUTINE cfg_write



   !!======================================================================
END MODULE domain
