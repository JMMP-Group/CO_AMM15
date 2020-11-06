MODULE trcbc
   !!======================================================================
   !!                     ***  MODULE  trcbc  ***
   !! TOP :  module for passive tracer boundary conditions
   !!=====================================================================
   !! History :  3.5 !  2014 (M. Vichi, T. Lovato)  Original
   !!            3.6 !  2015 (T . Lovato) Revision and BDY support
   !!            4.0 !  2016 (T . Lovato) Include application of sbc and cbc
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP model 
   !!----------------------------------------------------------------------
   !!   trc_bc       :  Apply tracer Boundary Conditions
   !!----------------------------------------------------------------------
   USE par_trc       !  passive tracers parameters
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables
   USE iom           !  I/O manager
   USE lib_mpp       !  MPP library
   USE fldread       !  read input fields
   USE bdy_oce,  ONLY: ln_bdy, nb_bdy , idx_bdy, ln_coords_file, rn_time_dmp, rn_time_dmp_out

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_bc         ! called in trcstp.F90 or within TOP modules
   PUBLIC   trc_bc_ini     ! called in trcini.F90 

   INTEGER  , SAVE, PUBLIC                             :: nb_trcobc    ! number of tracers with open BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trcsbc    ! number of tracers with surface BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trccbc    ! number of tracers with coastal BC
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indobc ! index of tracer with OBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indsbc ! index of tracer with SBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indcbc ! index of tracer with CBC data
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trsfac    ! multiplicative factor for SBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcsbc    ! structure of data input SBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trcfac    ! multiplicative factor for CBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trccbc    ! structure of data input CBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trofac    ! multiplicative factor for OBCtracer values
#if defined key_agrif
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcobc    ! structure of data input OBC (file informations, fields read)
#else
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:), TARGET  :: sf_trcobc
#endif

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcbc.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_bc_ini( ntrc )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_bc_ini  ***
      !!                    
      !! ** Purpose :   initialisation of passive tracer BC data 
      !! 
      !! ** Method  : - Read namtsd namelist
      !!              - allocates passive tracer BC data structure 
      !!----------------------------------------------------------------------
      INTEGER,INTENT(in) :: ntrc                           ! number of tracers
      !
      INTEGER            :: jl, jn , ib, ibd, ii, ij, ik   ! dummy loop indices
      INTEGER            :: ierr0, ierr1, ierr2, ierr3     ! temporary integers
      INTEGER            :: ios                            ! Local integer output status for namelist read
      INTEGER            :: nblen, igrd                    ! support arrays for BDY
      CHARACTER(len=100) :: clndta, clntrc
      !
      CHARACTER(len=100) :: cn_dir_sbc, cn_dir_cbc, cn_dir_obc
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) :: slf_i  ! local array of namelist informations on the fields to read
      TYPE(FLD_N), DIMENSION(jpmaxtrc) :: sn_trcobc    ! open
      TYPE(FLD_N), DIMENSION(jpmaxtrc) :: sn_trcsbc    ! surface
      TYPE(FLD_N), DIMENSION(jpmaxtrc) :: sn_trccbc    ! coastal
      REAL(wp)   , DIMENSION(jpmaxtrc) :: rn_trofac    ! multiplicative factor for tracer values
      REAL(wp)   , DIMENSION(jpmaxtrc) :: rn_trsfac    ! multiplicative factor for tracer values
      REAL(wp)   , DIMENSION(jpmaxtrc) :: rn_trcfac    ! multiplicative factor for tracer values
      !!
      NAMELIST/namtrc_bc/ cn_dir_obc, sn_trcobc, rn_trofac, cn_dir_sbc, sn_trcsbc, rn_trsfac, & 
                        & cn_dir_cbc, sn_trccbc, rn_trcfac, ln_rnf_ctl, rn_bc_time
      NAMELIST/namtrc_bdy/ cn_trc_dflt, cn_trc, nn_trcdmp_bdy
      !!----------------------------------------------------------------------
      !
      IF( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_bc_ini : Tracers Boundary Conditions (BC)'
         WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !  Initialisation and local array allocation
      ierr0 = 0   ;   ierr1 = 0   ;   ierr2 = 0   ;   ierr3 = 0  
      ALLOCATE( slf_i(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_ini: unable to allocate local slf_i' )   ;   RETURN
      ENDIF

      ! Compute the number of tracers to be initialised with open, surface and boundary data
      ALLOCATE( n_trc_indobc(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_ini: unable to allocate n_trc_indobc' )   ;   RETURN
      ENDIF
      nb_trcobc       = 0
      n_trc_indobc(:) = 0
      !
      ALLOCATE( n_trc_indsbc(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_ini: unable to allocate n_trc_indsbc' )   ;   RETURN
      ENDIF
      nb_trcsbc       = 0
      n_trc_indsbc(:) = 0
      !
      ALLOCATE( n_trc_indcbc(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_ini: unable to allocate n_trc_indcbc' )   ;   RETURN
      ENDIF
      nb_trccbc       = 0
      n_trc_indcbc(:) = 0
      !
      ! Read Boundary Conditions Namelists
      REWIND( numnat_ref )              ! Namelist namtrc_bc in reference namelist : Passive tracer data structure
      READ  ( numnat_ref, namtrc_bc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc_bc in reference namelist' )
      REWIND( numnat_cfg )              ! Namelist namtrc_bc in configuration namelist : Passive tracer data structure
      READ  ( numnat_cfg, namtrc_bc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_bc in configuration namelist' )
      IF(lwm) WRITE ( numont, namtrc_bc )

      IF ( ln_bdy ) THEN
         REWIND( numnat_ref )              ! Namelist namtrc_bdy in reference namelist : Passive tracer data structure
         READ  ( numnat_ref, namtrc_bdy, IOSTAT = ios, ERR = 903)
903      IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc_bdy in reference namelist' )
         ! make sur that all elements of the namelist variables have a default definition from namelist_ref
         cn_trc     (2:jp_bdy) = cn_trc     (1)
         cn_trc_dflt(2:jp_bdy) = cn_trc_dflt(1)
         REWIND( numnat_cfg )              ! Namelist namtrc_bdy in configuration namelist : Passive tracer data structure
         READ  ( numnat_cfg, namtrc_bdy, IOSTAT = ios, ERR = 904 )
904      IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_bdy in configuration namelist' )
         IF(lwm) WRITE ( numont, namtrc_bdy )
      
         ! setup up preliminary informations for BDY structure
         DO jn = 1, ntrc
            DO ib = 1, nb_bdy
               ! Set type of obc in BDY data structure (around here we may plug user override of obc type from nml)
               IF ( ln_trc_obc(jn) ) THEN   ;   trcdta_bdy(jn,ib)%cn_obc = TRIM( cn_trc     (ib) )
               ELSE                         ;   trcdta_bdy(jn,ib)%cn_obc = TRIM( cn_trc_dflt(ib) )
               ENDIF
               ! set damping use in BDY data structure
               trcdta_bdy(jn,ib)%dmp = .false.
               IF(nn_trcdmp_bdy(ib) == 1 .AND. ln_trc_obc(jn) )   trcdta_bdy(jn,ib)%dmp = .true.
               IF(nn_trcdmp_bdy(ib) == 2                      )   trcdta_bdy(jn,ib)%dmp = .true.
               IF(trcdta_bdy(jn,ib)%cn_obc == 'frs' .AND. nn_trcdmp_bdy(ib) /= 0 )  &
                   & CALL ctl_stop( 'trc_bc_ini: Use FRS OR relaxation' )
               !! NB change the < for <=
               IF(  .NOT.( 0 <= nn_trcdmp_bdy(ib)  .AND.  nn_trcdmp_bdy(ib) <= 2 )  )  THEN
                   CALL ctl_stop( 'trc_bc_ini: Not a valid option for nn_trcdmp_bdy. Allowed: 0,1,2.' )
               ENDIF
            END DO
         END DO
      ELSE
         ! Force all tracers OBC to false if bdy not used
         ln_trc_obc = .false.
      ENDIF

      ! compose BC data indexes
      DO jn = 1, ntrc
         IF( ln_trc_obc(jn) ) THEN
             nb_trcobc       = nb_trcobc + 1   ;   n_trc_indobc(jn) = nb_trcobc
         ENDIF
         IF( ln_trc_sbc(jn) ) THEN
             nb_trcsbc       = nb_trcsbc + 1   ;   n_trc_indsbc(jn) = nb_trcsbc
         ENDIF
         IF( ln_trc_cbc(jn) ) THEN
             nb_trccbc       = nb_trccbc + 1   ;   n_trc_indcbc(jn) = nb_trccbc
         ENDIF
      END DO

      ! Print summmary of Boundary Conditions
      IF( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,'(a,i3)') '   Total tracers to be initialized with SURFACE BCs data:', nb_trcsbc
         IF ( nb_trcsbc > 0 ) THEN
            WRITE(numout,*) '   #trc        NAME        Boundary     Mult.Fact. '
            DO jn = 1, ntrc
               IF ( ln_trc_sbc(jn) ) WRITE(numout,9001) jn, TRIM( sn_trcsbc(jn)%clvar ), 'SBC', rn_trsfac(jn)
            END DO
         ENDIF
         WRITE(numout,'(2a)') '   SURFACE BC data repository : ', TRIM(cn_dir_sbc)
         !
         WRITE(numout,*)
         WRITE(numout,'(a,i3)') '   Total tracers to be initialized with COASTAL BCs data:', nb_trccbc
         IF( nb_trccbc > 0 ) THEN
            WRITE(numout,*) '   #trc        NAME        Boundary     Mult.Fact. '
            DO jn = 1, ntrc
               IF ( ln_trc_cbc(jn) ) WRITE(numout, 9001) jn, TRIM( sn_trccbc(jn)%clvar ), 'CBC', rn_trcfac(jn)
            END DO
         ENDIF
         WRITE(numout,'(2a)') '   COASTAL BC data repository : ', TRIM(cn_dir_cbc)
         IF( .NOT.ln_rnf .OR. .NOT.ln_linssh )   ln_rnf_ctl = .FALSE.
         IF( ln_rnf_ctl )  WRITE(numout,'(a)') &
              &            ' -> Remove runoff dilution effect on tracers with absent river load (ln_rnf_ctl = .TRUE.)'
         WRITE(numout,*)
         WRITE(numout,'(a,i3)') '   Total tracers to be initialized with OPEN BCs data:', nb_trcobc

         IF( ln_bdy .AND. nb_trcobc > 0 ) THEN
            WRITE(numout,*) '   #trc        NAME        Boundary     Mult.Fact.   OBC Settings'
            DO jn = 1, ntrc
               IF (       ln_trc_obc(jn) )  WRITE(numout, 9001) jn, TRIM( sn_trcobc(jn)%clvar ), 'OBC', rn_trofac(jn), &
                    &                                           (trcdta_bdy(jn,ib)%cn_obc,ib=1,nb_bdy)
               IF ( .NOT. ln_trc_obc(jn) )  WRITE(numout, 9002) jn, 'Set data to IC and use default condition'       , &
                    &                                           (trcdta_bdy(jn,ib)%cn_obc,ib=1,nb_bdy)
            END DO
            WRITE(numout,*) ' '
            DO ib = 1, nb_bdy
               IF(nn_trcdmp_bdy(ib) == 0) WRITE(numout,9003) '   Boundary ', ib, &
                  &                                          ' -> NO damping of tracers'
               IF(nn_trcdmp_bdy(ib) == 1) WRITE(numout,9003) '   Boundary ', ib, &
                  &                                          ' -> damping ONLY for tracers with external data provided'
               IF(nn_trcdmp_bdy(ib) == 2) WRITE(numout,9003) '   Boundary ', ib, &
                  &                                          ' -> damping of ALL tracers'
               IF(nn_trcdmp_bdy(ib) >  0) THEN
                   WRITE(numout,9003) '     USE damping parameters from nambdy for boundary ', ib,' : '
                   WRITE(numout,'(a,f10.2,a)') '     - Inflow damping time scale  : ',rn_time_dmp    (ib),' days'
                   WRITE(numout,'(a,f10.2,a)') '     - Outflow damping time scale : ',rn_time_dmp_out(ib),' days'
               ENDIF
            END DO
         ENDIF
         !
         WRITE(numout,'(2a)') '   OPEN BC data repository : ', TRIM(cn_dir_obc)
      ENDIF
9001  FORMAT(2x,i5, 3x, a15, 3x, a5, 6x, e11.3, 4x, 10a13)
9002  FORMAT(2x,i5, 3x, a41, 3x, 10a13)
9003  FORMAT(a, i5, a)
      !
      !
      ! OPEN Lateral boundary conditions
      IF( ln_bdy .AND. nb_trcobc > 0 ) THEN 
         ALLOCATE ( sf_trcobc(nb_trcobc), rf_trofac(nb_trcobc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_ini: unable to allocate sf_trcobc structure' )   ;   RETURN
         ENDIF
         !
         igrd = 1                       ! Everything is at T-points here
         !
         DO jn = 1, ntrc
            DO ib = 1, nb_bdy
               !
               nblen = idx_bdy(ib)%nblen(igrd)
               !
               IF( ln_trc_obc(jn) ) THEN     !* Initialise from external data *!
                  jl = n_trc_indobc(jn)
                  slf_i(jl)    = sn_trcobc(jn)
                  rf_trofac(jl) = rn_trofac(jn)
                                                ALLOCATE( sf_trcobc(jl)%fnow(nblen,1,jpk)   , STAT=ierr2 )
                  IF( sn_trcobc(jn)%ln_tint )   ALLOCATE( sf_trcobc(jl)%fdta(nblen,1,jpk,2) , STAT=ierr3 )
                  IF( ierr2 + ierr3 > 0 ) THEN
                    CALL ctl_stop( 'trc_bc_ini : unable to allocate passive tracer OBC data arrays' )   ;   RETURN
                  ENDIF
                  trcdta_bdy(jn,ib)%trc => sf_trcobc(jl)%fnow(:,1,:)
                  trcdta_bdy(jn,ib)%rn_fac = rf_trofac(jl)
               ELSE                          !* Initialise obc arrays from initial conditions *!
                  ALLOCATE ( trcdta_bdy(jn,ib)%trc(nblen,jpk) )
                  DO ibd = 1, nblen
                     DO ik = 1, jpkm1
                        ii = idx_bdy(ib)%nbi(ibd,igrd)
                        ij = idx_bdy(ib)%nbj(ibd,igrd)
                        trcdta_bdy(jn,ib)%trc(ibd,ik) = trn(ii,ij,ik,jn) * tmask(ii,ij,ik)
                     END DO
                  END DO
                  trcdta_bdy(jn,ib)%rn_fac = 1._wp
               ENDIF
            END DO
         END DO
         !
         CALL fld_fill( sf_trcobc, slf_i, cn_dir_obc, 'trc_bc_ini', 'Passive tracer OBC data', 'namtrc_bc' )
         DO jn = 1, ntrc   ! define imap pointer, must be done after the call to fld_fill
            DO ib = 1, nb_bdy
               IF( ln_trc_obc(jn) ) THEN     !* Initialise from external data *!
                  jl = n_trc_indobc(jn)
                  sf_trcobc(jl)%imap => idx_bdy(ib)%nbmap(1:idx_bdy(ib)%nblen(igrd),igrd)
               ENDIF
            END DO
         END DO
         !
      ENDIF

      ! SURFACE Boundary conditions
      IF( nb_trcsbc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trcsbc(nb_trcsbc), rf_trsfac(nb_trcsbc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_ini: unable to allocate  sf_trcsbc structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, ntrc
            IF( ln_trc_sbc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indsbc(jn)
               slf_i(jl)    = sn_trcsbc(jn)
               rf_trsfac(jl) = rn_trsfac(jn)
                                            ALLOCATE( sf_trcsbc(jl)%fnow(jpi,jpj,1)   , STAT=ierr2 )
               IF( sn_trcsbc(jn)%ln_tint )  ALLOCATE( sf_trcsbc(jl)%fdta(jpi,jpj,1,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc_ini : unable to allocate passive tracer SBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         END DO
         !                         ! fill sf_trcsbc with slf_i and control print
         CALL fld_fill( sf_trcsbc, slf_i, cn_dir_sbc, 'trc_bc_ini', 'Passive tracer SBC data', 'namtrc_bc' )
         !
      ENDIF
      !
      ! COSTAL Boundary conditions
      IF( nb_trccbc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trccbc(nb_trccbc), rf_trcfac(nb_trccbc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_ini: unable to allocate  sf_trccbc structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, ntrc
            IF( ln_trc_cbc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indcbc(jn)
               slf_i(jl)    = sn_trccbc(jn)
               rf_trcfac(jl) = rn_trcfac(jn)
                                            ALLOCATE( sf_trccbc(jl)%fnow(jpi,jpj,1)   , STAT=ierr2 )
               IF( sn_trccbc(jn)%ln_tint )  ALLOCATE( sf_trccbc(jl)%fdta(jpi,jpj,1,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc_ini : unable to allocate passive tracer CBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         END DO
         !                         ! fill sf_trccbc with slf_i and control print
         CALL fld_fill( sf_trccbc, slf_i, cn_dir_cbc, 'trc_bc_ini', 'Passive tracer CBC data', 'namtrc_bc' )
         !
      ENDIF
      !
      DEALLOCATE( slf_i )          ! deallocate local field structure
      !
   END SUBROUTINE trc_bc_ini


   SUBROUTINE trc_bc(kt, jit)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_bc  ***
      !!
      !! ** Purpose :  Apply Boundary Conditions data to tracers
      !!
      !! ** Method  :  1) Read BC inputs and update data structures using fldread
      !!               2) Apply Boundary Conditions to tracers
      !!----------------------------------------------------------------------
      USE fldread
      !!      
      INTEGER, INTENT(in)           ::   kt    ! ocean time-step index
      INTEGER, INTENT(in), OPTIONAL ::   jit   ! subcycle time-step index (for timesplitting option)
      !!
      INTEGER  :: ji, jj, jk, jn, jl             ! Loop index
      REAL(wp) :: zfact, zrnf
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_bc')

      IF( kt == nit000 .AND. lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_bc : Surface boundary conditions for passive tracers.'
         WRITE(numout,*) '~~~~~~~ '
      ENDIF

      ! 1. Update Boundary conditions data
      IF( PRESENT(jit) ) THEN 
         !
         ! OPEN boundary conditions (use kt_offset=+1 as they are applied at the end of the step)
         IF( nb_trcobc > 0 ) THEN
           if (lwp) write(numout,'(a,i5,a,i10)') '   reading OBC data for ', nb_trcobc ,' variable(s) at step ', kt
           CALL fld_read( kt=kt, kn_fsbc=1, sd=sf_trcobc, kit=jit, kt_offset=+1)
         ENDIF
         !
         ! SURFACE boundary conditions
         IF( nb_trcsbc > 0 ) THEN
           if (lwp) write(numout,'(a,i5,a,i10)') '   reading SBC data for ', nb_trcsbc ,' variable(s) at step ', kt
           CALL fld_read( kt=kt, kn_fsbc=1, sd=sf_trcsbc, kit=jit)
         ENDIF
         !
         ! COASTAL boundary conditions
         IF( nb_trccbc > 0 ) THEN
           if (lwp) write(numout,'(a,i5,a,i10)') '   reading CBC data for ', nb_trccbc ,' variable(s) at step ', kt
           CALL fld_read( kt=kt, kn_fsbc=1, sd=sf_trccbc, kit=jit)
         ENDIF
         !
      ELSE
         !
         ! OPEN boundary conditions (use kt_offset=+1 as they are applied at the end of the step)
         IF( nb_trcobc > 0 ) THEN
           if (lwp) write(numout,'(a,i5,a,i10)') '   reading OBC data for ', nb_trcobc ,' variable(s) at step ', kt
           CALL fld_read( kt=kt, kn_fsbc=1, sd=sf_trcobc, kt_offset=+1)
         ENDIF
         !
         ! SURFACE boundary conditions
         IF( nb_trcsbc > 0 ) THEN
           if (lwp) write(numout,'(a,i5,a,i10)') '   reading SBC data for ', nb_trcsbc ,' variable(s) at step ', kt
           CALL fld_read( kt=kt, kn_fsbc=1, sd=sf_trcsbc )
         ENDIF
         !
         ! COASTAL boundary conditions
         IF( nb_trccbc > 0 ) THEN
           if (lwp) write(numout,'(a,i5,a,i10)') '   reading CBC data for ', nb_trccbc ,' variable(s) at step ', kt
           CALL fld_read( kt=kt, kn_fsbc=1, sd=sf_trccbc )
         ENDIF
         !
      ENDIF

      ! 2. Apply Boundary conditions data
      ! 
      DO jn = 1 , jptra
         !
         ! Remove river dilution for tracers with absent river load
         IF( ln_rnf_ctl .AND. .NOT.ln_trc_cbc(jn) ) THEN
            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1
                  DO jk = 1, nk_rnf(ji,jj)
                     zrnf = (rnf(ji,jj) + rnf_b(ji,jj)) * 0.5_wp * r1_rau0 / h_rnf(ji,jj)
                     tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn)  + (trn(ji,jj,jk,jn) * zrnf)
                  END DO
               END DO
            END DO
         ENDIF
         !
         ! OPEN boundary conditions: trcbdy is called in trcnxt !
         !
         ! SURFACE boundary conditions
         IF( ln_trc_sbc(jn) ) THEN
            jl = n_trc_indsbc(jn)
            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zfact = 1. / ( e3t_n(ji,jj,1) * rn_bc_time )
                  tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + rf_trsfac(jl) * sf_trcsbc(jl)%fnow(ji,jj,1) * zfact
               END DO
            END DO
         ENDIF
         !
         ! COASTAL boundary conditions
         IF( ln_rnf .AND. ln_trc_cbc(jn) ) THEN
            jl = n_trc_indcbc(jn)
            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  DO jk = 1, nk_rnf(ji,jj)
                     zfact = rn_rfact / ( e1e2t(ji,jj) * h_rnf(ji,jj) * rn_bc_time ) 
                     tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + rf_trcfac(jl) * sf_trccbc(jl)%fnow(ji,jj,1) * zfact
                  END DO
               END DO
            END DO
         ENDIF
         !                                                       ! ===========
      END DO                                                     ! tracer loop
      !                                                          ! ===========
      IF( ln_timing )   CALL timing_stop('trc_bc')
      !
   END SUBROUTINE trc_bc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_bc_ini( ntrc )        ! Empty routine
      INTEGER,INTENT(IN) :: ntrc                           ! number of tracers
      WRITE(*,*) 'trc_bc_ini: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bc_ini
   SUBROUTINE trc_bc( kt )        ! Empty routine
      WRITE(*,*) 'trc_bc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bc
#endif

   !!======================================================================
END MODULE trcbc
