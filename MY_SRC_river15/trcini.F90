MODULE trcini
   !!======================================================================
   !!                         ***  MODULE trcini  ***
   !! TOP :   Manage the passive tracer initialization
   !!======================================================================
   !! History :   -   ! 1991-03 (O. Marti)  original code
   !!            1.0  ! 2005-03 (O. Aumont, A. El Moussaoui) F90
   !!            2.0  ! 2005-10 (C. Ethe, G. Madec) revised architecture
   !!            4.0  ! 2011-01 (A. R. Porter, STFC Daresbury) dynamical allocation
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_init  :   Initialization for passive tracer
   !!   top_alloc :   allocate the TOP arrays
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables
   USE trcnam          ! Namelist read
   USE daymod          ! calendar manager
   USE prtctl_trc      ! Print control passive tracers (prt_ctl_trc_init routine)
   USE trcsub          ! variables to substep passive tracers
   USE trcrst
   USE lib_mpp         ! distribued memory computing library
   USE trcice          ! tracers in sea ice
   USE trcbc,   only : trc_bc_ini ! generalized Boundary Conditions
 
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_init   ! called by opa

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini.F90 10817 2019-03-29 17:23:45Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE trc_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_init  ***
      !!
      !! ** Purpose :   Initialization of the passive tracer fields 
      !!
      !! ** Method  : - read namelist
      !!              - control the consistancy 
      !!              - compute specific initialisations
      !!              - set initial tracer fields (either read restart 
      !!                or read data or analytical formulation
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_init')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_init : initial set up of the passive tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      !
      CALL trc_ini_ctl   ! control 
      CALL trc_nam       ! read passive tracers namelists
      CALL top_alloc()   ! allocate TOP arrays
      !
      IF(.NOT.ln_trcdta )   ln_trc_ini(:) = .FALSE.
      !
      IF(lwp) WRITE(numout,*)
      IF( ln_rsttr .AND. .NOT. l_offline ) CALL trc_rst_cal( nit000, 'READ' )   ! calendar
      IF(lwp) WRITE(numout,*)
      !
      CALL trc_ini_sms   ! SMS
      CALL trc_ini_trp   ! passive tracers transport
      CALL trc_ice_ini   ! Tracers in sea ice
      !
      IF( lwm .AND. sn_cfctl%l_trcstat ) THEN
         CALL ctl_opn( numstr, 'tracer.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp , narea )
      ENDIF
      !
      CALL trc_ini_state  !  passive tracers initialisation : from a restart or from clim
      IF( nn_dttrc /= 1 ) &
      CALL trc_sub_ini    ! Initialize variables for substepping passive tracers
      !
      CALL trc_ini_inv   ! Inventories
      !


      IF( ln_timing )   CALL timing_stop('trc_init')
      !
   END SUBROUTINE trc_init


   SUBROUTINE trc_ini_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ini_ctl  ***
      !! ** Purpose :        Control  + ocean volume
      !!----------------------------------------------------------------------
      INTEGER ::   jk    ! dummy loop indices
      !
      ! Define logical parameter ton control dirunal cycle in TOP
      l_trcdm2dc = ln_dm2dc .OR. ( ln_cpl .AND. ncpl_qsr_freq /= 1 )
      l_trcdm2dc = l_trcdm2dc  .AND. .NOT. l_offline
      IF( l_trcdm2dc .AND. lwp )   CALL ctl_warn( 'Coupling with passive tracers and used of diurnal cycle.',   &
         &                           'Computation of a daily mean shortwave for some biogeochemical models ' )
      !
   END SUBROUTINE trc_ini_ctl


   SUBROUTINE trc_ini_inv
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ini_stat  ***
      !! ** Purpose :      passive tracers inventories at initialsation phase
      !!----------------------------------------------------------------------
      INTEGER ::  jk, jn    ! dummy loop indices
      CHARACTER (len=25) :: charout
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_ini_inv : initial passive tracers inventories'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !                          ! masked grid volume
      DO jk = 1, jpk
         cvol(:,:,jk) = e1e2t(:,:) * e3t_n(:,:,jk) * tmask(:,:,jk)
      END DO
      !                          ! total volume of the ocean 
      areatot = glob_sum( 'trcini', cvol(:,:,:) )
      !
      trai(:) = 0._wp            ! initial content of all tracers
      DO jn = 1, jptra
         trai(jn) = trai(jn) + glob_sum( 'trcini', trn(:,:,:,jn) * cvol(:,:,:)   )
      END DO

      IF(lwp) THEN               ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   ==>>>   Total number of passive tracer jptra = ', jptra
         WRITE(numout,*) '           Total volume of ocean                = ', areatot
         WRITE(numout,*) '           Total inital content of all tracers '
         WRITE(numout,*)
         DO jn = 1, jptra
            WRITE(numout,9000) jn, TRIM( ctrcnm(jn) ), trai(jn)
         ENDDO
         WRITE(numout,*)
      ENDIF
      IF(lwp) WRITE(numout,*)
      IF(ln_ctl) THEN            ! print mean trends (used for debugging)
         CALL prt_ctl_trc_init
         WRITE(charout, FMT="('ini ')")
         CALL prt_ctl_trc_info( charout )
         CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF
9000  FORMAT('      tracer nb : ',i2,'      name :',a10,'      initial content :',e18.10)
      !
   END SUBROUTINE trc_ini_inv


   SUBROUTINE trc_ini_sms
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ini_sms  ***
      !! ** Purpose :   SMS initialisation
      !!----------------------------------------------------------------------
      USE trcini_pisces  ! PISCES   initialisation
      USE trcini_cfc     ! CFC      initialisation
      USE trcini_c14     ! C14  initialisation
      USE trcini_age     ! age initialisation
      USE trcini_my_trc  ! MY_TRC   initialisation
      !
      INTEGER :: jn
      !!----------------------------------------------------------------------
      !
      ! Pass sn_tracer fields to specialized arrays 
      DO jn = 1, jp_bgc
         ctrcnm    (jn) = TRIM( sn_tracer(jn)%clsname )
         ctrcln    (jn) = TRIM( sn_tracer(jn)%cllname )
         ctrcun    (jn) = TRIM( sn_tracer(jn)%clunit  )
         ln_trc_ini(jn) =       sn_tracer(jn)%llinit
         ln_trc_sbc(jn) =       sn_tracer(jn)%llsbc
         ln_trc_cbc(jn) =       sn_tracer(jn)%llcbc
         ln_trc_obc(jn) =       sn_tracer(jn)%llobc
      END DO
      !    
      IF( ln_pisces      )   CALL trc_ini_pisces     !  PISCES model
      IF( ln_my_trc      )   CALL trc_ini_my_trc     !  MY_TRC model
      IF( ll_cfc         )   CALL trc_ini_cfc        !  CFC's
      IF( ln_c14         )   CALL trc_ini_c14        !  C14 model
      IF( ln_age         )   CALL trc_ini_age        !  AGE
      !
      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'trc_init_sms : Summary for selected passive tracers'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '    ID     NAME     INI  SBC  CBC  OBC'
         DO jn = 1, jptra
            WRITE(numout,9001) jn, TRIM(ctrcnm(jn)), ln_trc_ini(jn), ln_trc_sbc(jn),ln_trc_cbc(jn),ln_trc_obc(jn)
         END DO
      ENDIF
9001  FORMAT(3x,i3,1x,a10,3x,l2,3x,l2,3x,l2,3x,l2)
      !
   END SUBROUTINE trc_ini_sms


   SUBROUTINE trc_ini_trp
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ini_trp  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!----------------------------------------------------------------------
      USE trcdmp , ONLY:  trc_dmp_ini
      USE trcadv , ONLY:  trc_adv_ini
      USE trcldf , ONLY:  trc_ldf_ini
      USE trcrad , ONLY:  trc_rad_ini
      USE trcsink, ONLY:  trc_sink_ini
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      IF( ln_trcdmp )  CALL  trc_dmp_ini          ! damping
                       CALL  trc_adv_ini          ! advection
                       CALL  trc_ldf_ini          ! lateral diffusion
                       !                          ! vertical diffusion: always implicit time stepping scheme
                       CALL  trc_rad_ini          ! positivity of passive tracers 
                       CALL  trc_sink_ini         ! Vertical sedimentation of particles
      !
   END SUBROUTINE trc_ini_trp


   SUBROUTINE trc_ini_state
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ini_state ***
      !! ** Purpose :          Initialisation of passive tracer concentration 
      !!----------------------------------------------------------------------
      USE zpshde          ! partial step: hor. derivative   (zps_hde routine)
      USE trcrst          ! passive tracers restart
      USE trcdta          ! initialisation from files
      !
      INTEGER :: jn, jl   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_trcdta )   CALL trc_dta_ini( jptra )      ! set initial tracers values
      !
      IF( ln_my_trc )   CALL trc_bc_ini ( jptra )      ! set tracers Boundary Conditions
      !
      !
      !-- NB comment restart read here as we want to 
      !      load flags at every restart to know where 
      !      to apply mass flux (river mouth or random places)
!      IF( ln_rsttr ) THEN              ! restart from a file
!        !
!        CALL trc_rst_read
!        !
!      ELSE                             ! Initialisation of tracer from a file that may also be used for damping
!!gm BUG ?   if damping and restart, what's happening ?
        IF( ln_trcdta .AND. nb_trcdta > 0 ) THEN
            ! update passive tracers arrays with input data read from file
            DO jn = 1, jptra
               IF( ln_trc_ini(jn) ) THEN
                  jl = n_trc_index(jn) 
                  CALL trc_dta( nit000, sf_trcdta(jl), rf_trfac(jl), trn(:,:,:,jn) )
                  !
                  ! deallocate data structure if data are not used for damping
                  IF( .NOT.ln_trcdmp .AND. .NOT.ln_trcdmp_clo ) THEN
                     IF(lwp) WRITE(numout,*) 'trc_ini_state: deallocate data arrays as they are only used to initialize the run'
                                                  DEALLOCATE( sf_trcdta(jl)%fnow )
                     IF( sf_trcdta(jl)%ln_tint )  DEALLOCATE( sf_trcdta(jl)%fdta )
                     !
                  ENDIF
               ENDIF
            END DO
            !
        ENDIF
        !
        !--- NB -- STORE LOCATION WHERE TO APPLY THE MASS FLUX
        flagRIV(:,:,:,:) = 0.              ! init flag to 0
        flagOTH(:,:,:,:) = 0.              ! init flag to 0
        WHERE( trn .gt. 0. ) flagRIV = 1.  ! 1 to state where there are rivers (runoff used for mass flux)
        WHERE( trn .lt. 0. ) flagOTH = 1.  ! 1 to state where there are other sources (velocity used for mass flux)
        trn(:,:,:,:) = 0.                  ! turn off input - now mass flux
        !--- END NB 
        !
        trb(:,:,:,:) = trn(:,:,:,:)
        !
      !-- NB - READ RESTART HERE AFTER INITIATING THE CONSTANT TRACER
      IF( ln_rsttr ) THEN
          CALL trc_rst_read
      ENDIF
      !-- END NB -----------------------------------------------------
! END NB      ENDIF
      !
      tra(:,:,:,:) = 0._wp
      !                                                         ! Partial top/bottom cell: GRADh(trn)
   END SUBROUTINE trc_ini_state


   SUBROUTINE top_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE top_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!----------------------------------------------------------------------
      USE trc           , ONLY:   trc_alloc
      USE trdtrc_oce    , ONLY:   trd_trc_oce_alloc
#if defined key_trdmxl_trc 
      USE trdmxl_trc    , ONLY:   trd_mxl_trc_alloc
#endif
      !
      INTEGER ::   ierr   ! local integer
      !!----------------------------------------------------------------------
      !
      ierr =        trc_alloc()
      ierr = ierr + trd_trc_oce_alloc()
#if defined key_trdmxl_trc 
      ierr = ierr + trd_mxl_trc_alloc()
#endif
      !
      CALL mpp_sum( 'trcini', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'top_alloc : unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE top_alloc

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_init                      ! Dummy routine   
   END SUBROUTINE trc_init
#endif

   !!======================================================================
END MODULE trcini
