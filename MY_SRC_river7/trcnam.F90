MODULE trcnam
   !!======================================================================
   !!                       ***  MODULE trcnam  ***
   !! TOP :   Read and print options for the passive tracer run (namelist)
   !!======================================================================
   !! History :    -   !  1996-11  (M.A. Foujols, M. Levy)  original code
   !!              -   !  1998-04  (M.A Foujols, L. Bopp) ahtrb0 for isopycnal mixing
   !!              -   !  1999-10  (M.A. Foujols, M. Levy) separation of sms
   !!              -   !  2000-07  (A. Estublier) add TVD and MUSCL : Tests on ndttrc
   !!              -   !  2000-11  (M.A Foujols, E Kestenare) trcrat, ahtrc0 and aeivtr0
   !!              -   !  2001-01 (E Kestenare) suppress ndttrc=1 for CEN2 and TVD schemes
   !!             1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_nam    :  Read and print options for the passive tracer run (namelist)
   !!----------------------------------------------------------------------
   USE oce_trc     ! shared variables between ocean and passive tracers
   USE trc         ! passive tracers common variables
   USE trd_oce     !       
   USE trdtrc_oce  !
   USE iom         ! I/O manager
#if defined key_mpp_mpi
   USE lib_mpp, ONLY: ncom_dttrc
#endif

   IMPLICIT NONE
   PRIVATE 

   PUBLIC   trc_nam_run  ! called in trcini
   PUBLIC   trc_nam      ! called in trcini

   TYPE(PTRACER), DIMENSION(jpmaxtrc), PUBLIC  :: sn_tracer  !: type of tracer for saving if not key_iomput

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnam.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_nam
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam  ***
      !!
      !! ** Purpose :   READ and PRINT options for the passive tracer run (namelist) 
      !!
      !! ** Method  : - read passive tracer namelist 
      !!              - read namelist of each defined SMS model
      !!                ( (PISCES, CFC, MY_TRC )
      !!---------------------------------------------------------------------
      INTEGER  ::   jn   ! dummy loop indice
      !!---------------------------------------------------------------------
      !
      IF( .NOT.l_offline )   CALL trc_nam_run     ! Parameters of the run                                  
      !               
      CALL trc_nam_trc                            ! passive tracer informations
      !                                        
      ! NB : we actually want to keep realising tracer even after a restart 
      !IF( ln_rsttr                     )   ln_trcdta = .FALSE.   ! restart : no need of clim data
      ! END NB
      !
      IF( ln_trcdmp .OR. ln_trcdmp_clo )   ln_trcdta = .TRUE.    ! damping : need to have clim data
      !
      !
      IF(lwp) THEN                   ! control print
         IF( ln_rsttr ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   ==>>>   Read a restart file for passive tracer : ', TRIM( cn_trcrst_in )
         ENDIF
         IF( ln_trcdta .AND. .NOT.ln_rsttr ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   ==>>>   Some of the passive tracers are initialised from climatologies '
         ENDIF
         IF( .NOT.ln_trcdta ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   ==>>>   All the passive tracers are initialised with constant values '
         ENDIF
      ENDIF
      !
      rdttrc = rdt * FLOAT( nn_dttrc )          ! passive tracer time-step      
      ! 
      IF(lwp) THEN                              ! control print
        WRITE(numout,*) 
        WRITE(numout,*) '   ==>>>   Passive Tracer  time step    rdttrc = nn_dttrc*rdt = ', rdttrc
      ENDIF
      !
      IF( l_trdtrc )        CALL trc_nam_trd    ! Passive tracer trends
      !
   END SUBROUTINE trc_nam


   SUBROUTINE trc_nam_run
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam  ***
      !!
      !! ** Purpose :   read options for the passive tracer run (namelist) 
      !!
      !!---------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer
      !!
      NAMELIST/namtrc_run/ nn_dttrc, ln_rsttr, nn_rsttr, ln_top_euler, &
        &                  cn_trcrst_indir, cn_trcrst_outdir, cn_trcrst_in, cn_trcrst_out
      !!---------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_nam_run : read the passive tracer namelists'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      CALL ctl_opn( numnat_ref, 'namelist_top_ref'   , 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnat_cfg, 'namelist_top_cfg'   , 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numont, 'output.namelist.top', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., 1 )
      !
      REWIND( numnat_ref )              ! Namelist namtrc in reference namelist : Passive tracer variables
      READ  ( numnat_ref, namtrc_run, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc in reference namelist' )
      REWIND( numnat_cfg )              ! Namelist namtrc in configuration namelist : Passive tracer variables
      READ  ( numnat_cfg, namtrc_run, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc in configuration namelist' )
      IF(lwm) WRITE( numont, namtrc_run )

      nittrc000 = nit000 + nn_dttrc - 1      ! first time step of tracer model

      IF(lwp) THEN                   ! control print
         WRITE(numout,*) '   Namelist : namtrc_run'
         WRITE(numout,*) '      time step freq. for passive tracer           nn_dttrc      = ', nn_dttrc
         WRITE(numout,*) '      restart  for passive tracer                  ln_rsttr      = ', ln_rsttr
         WRITE(numout,*) '      control of time step for passive tracer      nn_rsttr      = ', nn_rsttr
         WRITE(numout,*) '      first time step for pass. trac.              nittrc000     = ', nittrc000
         WRITE(numout,*) '      Use euler integration for TRC (y/n)          ln_top_euler  = ', ln_top_euler
      ENDIF
      !
#if defined key_mpp_mpi
      ncom_dttrc = nn_dttrc    ! make nn_fsbc available for lib_mpp
#endif
      !
   END SUBROUTINE trc_nam_run


   SUBROUTINE trc_nam_trc
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam  ***
      !!
      !! ** Purpose :   read options for the passive tracer run (namelist) 
      !!
      !!---------------------------------------------------------------------
      INTEGER ::   ios, ierr, icfc       ! Local integer
      !!
      NAMELIST/namtrc/jp_bgc, ln_pisces, ln_my_trc, ln_age, ln_cfc11, ln_cfc12, ln_sf6, ln_c14, &
         &            sn_tracer, ln_trcdta, ln_trcdmp, ln_trcdmp_clo, jp_dia3d, jp_dia2d
      !!---------------------------------------------------------------------
      ! Dummy settings to fill tracers data structure
      !                  !   name   !   title   !   unit   !   init  !   sbc   !   cbc   !   obc  !
      sn_tracer = PTRACER( 'NONAME' , 'NOTITLE' , 'NOUNIT' , .false. , .false. , .false. , .false.)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_nam_trc : read the passive tracer namelists'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      REWIND( numnat_ref )              ! Namelist namtrc in reference namelist : Passive tracer variables
      READ  ( numnat_ref, namtrc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc in reference namelist' )
      REWIND( numnat_cfg )              ! Namelist namtrc in configuration namelist : Passive tracer variables
      READ  ( numnat_cfg, namtrc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc in configuration namelist' )
      IF(lwm) WRITE( numont, namtrc )

      ! Control settings
      IF( ln_pisces .AND. ln_my_trc )   CALL ctl_stop( 'Choose only ONE BGC model - PISCES or MY_TRC' )
      IF( .NOT. ln_pisces .AND. .NOT. ln_my_trc )   jp_bgc = 0
      ll_cfc = ln_cfc11 .OR. ln_cfc12 .OR. ln_sf6
      !
      jptra       =  0
      jp_pisces   =  0    ;   jp_pcs0  =  0    ;   jp_pcs1  = 0
      jp_my_trc   =  0    ;   jp_myt0  =  0    ;   jp_myt1  = 0
      jp_cfc      =  0    ;   jp_cfc0  =  0    ;   jp_cfc1  = 0
      jp_age      =  0    ;   jp_c14   =  0
      !
      IF( ln_pisces )  THEN
         jp_pisces = jp_bgc
         jp_pcs0   = 1
         jp_pcs1   = jp_pisces
      ENDIF
      IF( ln_my_trc )  THEN
          jp_my_trc = jp_bgc
          jp_myt0   = 1
          jp_myt1   = jp_my_trc
      ENDIF
      !
      jptra  = jp_bgc
      !
      IF( ln_age )    THEN
      ! NB: we want an age per tracer
      !   jptra     = jptra + 1
         jp_age    = jptra
         !jp_age    =2 
         jptra     = jptra + jp_age
         
      ! END NB
      ENDIF
      IF( ln_cfc11 )  jp_cfc = jp_cfc + 1
      IF( ln_cfc12 )  jp_cfc = jp_cfc + 1
      IF( ln_sf6   )  jp_cfc = jp_cfc + 1
      IF( ll_cfc )    THEN
          jptra     = jptra + jp_cfc
          jp_cfc0   = jptra - jp_cfc + 1
          jp_cfc1   = jptra
      ENDIF
      IF( ln_c14 )    THEN
           jptra     = jptra + 1
           jp_c14    = jptra
      ENDIF
      !
      IF( jptra == 0 )   CALL ctl_stop( 'All TOP tracers disabled: change namtrc setting or check if key_top is active' )
      !
      IF(lwp) THEN                   ! control print
         WRITE(numout,*) '   Namelist : namtrc'
         WRITE(numout,*) '      Total number of passive tracers              jptra         = ', jptra
         WRITE(numout,*) '      Total number of BGC tracers                  jp_bgc        = ', jp_bgc
         WRITE(numout,*) '      Simulating PISCES model                      ln_pisces     = ', ln_pisces
         WRITE(numout,*) '      Simulating MY_TRC  model                     ln_my_trc     = ', ln_my_trc
         WRITE(numout,*) '      Simulating water mass age                    ln_age        = ', ln_age
         WRITE(numout,*) '      Simulating CFC11 passive tracer              ln_cfc11      = ', ln_cfc11
         WRITE(numout,*) '      Simulating CFC12 passive tracer              ln_cfc12      = ', ln_cfc12
         WRITE(numout,*) '      Simulating SF6 passive tracer                ln_sf6        = ', ln_sf6
         WRITE(numout,*) '      Total number of CFCs tracers                 jp_cfc        = ', jp_cfc
         WRITE(numout,*) '      Simulating C14   passive tracer              ln_c14        = ', ln_c14
         WRITE(numout,*) '      Read inputs data from file (y/n)             ln_trcdta     = ', ln_trcdta
         WRITE(numout,*) '      Damping of passive tracer (y/n)              ln_trcdmp     = ', ln_trcdmp
         WRITE(numout,*) '      Restoring of tracer on closed seas           ln_trcdmp_clo = ', ln_trcdmp_clo
      ENDIF
      !
      IF( ll_cfc .OR. ln_c14 ) THEN
        !                             ! Open namelist files
        CALL ctl_opn( numtrc_ref, 'namelist_trc_ref'   ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
        CALL ctl_opn( numtrc_cfg, 'namelist_trc_cfg'   ,     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
        IF(lwm) CALL ctl_opn( numonr, 'output.namelist.trc', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
        !
      ENDIF
      !
   END SUBROUTINE trc_nam_trc


   SUBROUTINE trc_nam_trd
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam_dia  ***
      !!
      !! ** Purpose :   read options for the passive tracer diagnostics
      !!
      !! ** Method  : - read passive tracer namelist 
      !!              - read namelist of each defined SMS model
      !!                ( (PISCES, CFC, MY_TRC )
      !!---------------------------------------------------------------------
#if defined key_trdmxl_trc  || defined key_trdtrc
      INTEGER  ::   ios, ierr                 ! Local integer
      !!
      NAMELIST/namtrc_trd/ nn_trd_trc, nn_ctls_trc, rn_ucf_trc, &
         &                ln_trdmxl_trc_restart, ln_trdmxl_trc_instant, &
         &                cn_trdrst_trc_in, cn_trdrst_trc_out, ln_trdtrc
      !!---------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_nam_trd : read the passive tracer diagnostics options'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      ALLOCATE( ln_trdtrc(jptra) ) 
      !
      REWIND( numnat_ref )              ! Namelist namtrc_trd in reference namelist : Passive tracer trends
      READ  ( numnat_ref, namtrc_trd, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc_trd in reference namelist' )
      REWIND( numnat_cfg )              ! Namelist namtrc_trd in configuration namelist : Passive tracer trends
      READ  ( numnat_cfg, namtrc_trd, IOSTAT = ios, ERR = 906 )
906   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_trd in configuration namelist' )
      IF(lwm) WRITE( numont, namtrc_trd )

      IF(lwp) THEN
         WRITE(numout,*) '   Namelist : namtrc_trd                    '
         WRITE(numout,*) '      frequency of trends diagnostics   nn_trd_trc             = ', nn_trd_trc
         WRITE(numout,*) '      control surface type              nn_ctls_trc            = ', nn_ctls_trc
         WRITE(numout,*) '      restart for ML diagnostics        ln_trdmxl_trc_restart  = ', ln_trdmxl_trc_restart
         WRITE(numout,*) '      instantantaneous or mean trends   ln_trdmxl_trc_instant  = ', ln_trdmxl_trc_instant
         WRITE(numout,*) '      unit conversion factor            rn_ucf_trc             = ', rn_ucf_trc
         DO jn = 1, jptra
            IF( ln_trdtrc(jn) ) WRITE(numout,*) '      compute ML trends for tracer number :', jn
         END DO
      ENDIF
#endif
      !
   END SUBROUTINE trc_nam_trd

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam                      ! Empty routine   
   END SUBROUTINE trc_nam
   SUBROUTINE trc_nam_run                      ! Empty routine   
   END SUBROUTINE trc_nam_run
#endif

   !!======================================================================
END MODULE trcnam
