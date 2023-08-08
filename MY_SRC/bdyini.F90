MODULE bdyini
   !!======================================================================
   !!                       ***  MODULE  bdyini  ***
   !! Unstructured open boundaries : initialisation
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-01  (D. Storkey) Tidal forcing
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) updates for Shelf configurations
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.4  !  2012     (J. Chanut) straight open boundary case update
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) optimization of BDY communications
   !!            3.7  !  2016     (T. Lovato) Remove bdy macro, call here init for dta and tides
   !!----------------------------------------------------------------------
   !!   bdy_init      : Initialization of unstructured open boundaries
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce , ONLY: nn_ice
   USE bdy_oce        ! unstructured open boundary conditions
   USE bdydta         ! open boundary cond. setting   (bdy_dta_init routine)
   USE bdytides       ! open boundary cond. setting   (bdytide_init routine)
   USE tide_mod, ONLY: ln_tide ! tidal forcing
   USE phycst  , ONLY: rday
   !
   USE in_out_manager ! I/O units
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! for mpp_sum  
   USE iom            ! I/O

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_init    ! routine called in nemo_init
   PUBLIC   find_neib   ! routine called in bdy_nmn

   INTEGER, PARAMETER ::   jp_nseg = 100   ! 
   ! Straight open boundary segment parameters:
   INTEGER  ::   nbdysege, nbdysegw, nbdysegn, nbdysegs 
   INTEGER, DIMENSION(jp_nseg) ::   jpieob, jpjedt, jpjeft, npckge   !
   INTEGER, DIMENSION(jp_nseg) ::   jpiwob, jpjwdt, jpjwft, npckgw   !
   INTEGER, DIMENSION(jp_nseg) ::   jpjnob, jpindt, jpinft, npckgn   !
   INTEGER, DIMENSION(jp_nseg) ::   jpjsob, jpisdt, jpisft, npckgs   !
   
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdyini.F90 15368 2021-10-14 08:25:34Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_init  ***
      !!
      !! ** Purpose :   Initialization of the dynamics and tracer fields with
      !!              unstructured open boundaries.
      !!
      !! ** Method  :   Read initialization arrays (mask, indices) to identify
      !!              an unstructured open boundary
      !!
      !! ** Input   :  bdy_init.nc, input file for unstructured open boundaries
      !!----------------------------------------------------------------------
      NAMELIST/nambdy/ ln_bdy, nb_bdy, ln_coords_file, cn_coords_file,         &
         &             ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta,     &
         &             cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta,             &
         &             ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, &
         &             cn_ice, nn_ice_dta,                                     &
         &             ln_vol, nn_volctl, nn_rimwidth
         !
      ! RDP boundary shift of ssh
      NAMELIST/nambdy_ssh/ ln_ssh_bdy, rn_ssh_shift
      INTEGER  ::   ib_bdy              ! dummy loop indices
      ! END RDP
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      ! ------------------------
      ! Read namelist parameters
      ! ------------------------
      READ  ( numnam_ref, nambdy, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy in reference namelist' )
      ! make sur that all elements of the namelist variables have a default definition from namelist_ref
      ln_coords_file (2:jp_bdy) = ln_coords_file (1)
      cn_coords_file (2:jp_bdy) = cn_coords_file (1)
      cn_dyn2d       (2:jp_bdy) = cn_dyn2d       (1)
      nn_dyn2d_dta   (2:jp_bdy) = nn_dyn2d_dta   (1)
      cn_dyn3d       (2:jp_bdy) = cn_dyn3d       (1)
      nn_dyn3d_dta   (2:jp_bdy) = nn_dyn3d_dta   (1)
      cn_tra         (2:jp_bdy) = cn_tra         (1)
      nn_tra_dta     (2:jp_bdy) = nn_tra_dta     (1)    
      ln_tra_dmp     (2:jp_bdy) = ln_tra_dmp     (1)
      ln_dyn3d_dmp   (2:jp_bdy) = ln_dyn3d_dmp   (1)
      rn_time_dmp    (2:jp_bdy) = rn_time_dmp    (1)
      rn_time_dmp_out(2:jp_bdy) = rn_time_dmp_out(1)
      cn_ice         (2:jp_bdy) = cn_ice         (1)
      nn_ice_dta     (2:jp_bdy) = nn_ice_dta     (1)
      nn_rimwidth    (2:jp_bdy) = nn_rimwidth    (1)
      READ  ( numnam_cfg, nambdy, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy in configuration namelist' )
      IF(lwm) WRITE ( numond, nambdy )

      ! RDP boundary shift of ssh
      READ  ( numnam_ref, nambdy_ssh, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambdy_ssh in reference namelist' )

      READ  ( numnam_cfg, nambdy_ssh, IOSTAT = ios, ERR = 906)
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambdy_ssh in configuration namelist' )
      IF(lwm) WRITE ( numond, nambdy_ssh )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'nambdy_ssh : use of ssh boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      IF(lwp) WRITE(numout,*) '      ln_ssh_bdy: '
      DO ib_bdy = 1,nb_bdy
        IF(lwp) WRITE(numout,*) '      ln_ssh_bdy  (',ib_bdy,'): ',ln_ssh_bdy(ib_bdy)
      IF(lwp) WRITE(numout,*) '      rn_ssh_shift: '
      ENDDO
      DO ib_bdy = 1,nb_bdy
        IF(lwp) WRITE(numout,*) '      rn_ssh_shift(',ib_bdy,'): ',rn_ssh_shift(ib_bdy)
      ENDDO
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      ! END RDP


      IF( .NOT. Agrif_Root() ) ln_bdy = .FALSE.   ! forced for Agrif children

      IF( nb_bdy == 0 ) ln_bdy = .FALSE.
      
      ! -----------------------------------------
      ! unstructured open boundaries use control
      ! -----------------------------------------
      IF ( ln_bdy ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'bdy_init : initialization of open boundaries'
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         !
         ! Open boundaries definition (arrays and masks)
         CALL bdy_def
         IF( ln_meshmask )   CALL bdy_meshwri()
         !
         ! Open boundaries initialisation of external data arrays
         CALL bdy_dta_init
         !
         ! Open boundaries initialisation of tidal harmonic forcing
         IF( ln_tide ) CALL bdytide_init
         !
      ELSE
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'bdy_init : open boundaries not used (ln_bdy = F)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         !
      ENDIF
      !
   END SUBROUTINE bdy_init


   SUBROUTINE bdy_def
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_init  ***
      !!         
      !! ** Purpose :   Definition of unstructured open boundaries.
      !!
      !! ** Method  :   Read initialization arrays (mask, indices) to identify 
      !!              an unstructured open boundary
      !!
      !! ** Input   :  bdy_init.nc, input file for unstructured open boundaries
      !!----------------------------------------------------------------------      
      INTEGER  ::   ji, jj                                 ! dummy loop indices
      INTEGER  ::   ib_bdy, ii, ij, igrd, ib, ir, iseg     ! dummy loop indices
      INTEGER  ::   icount, icountr, icountr0, ibr_max     ! local integers
      INTEGER  ::   ilen1                           !   -       -
      INTEGER  ::   iiRst, iiRnd, iiSst, iiSnd, iiSstdiag, iiSnddiag, iiSstsono, iiSndsono
      INTEGER  ::   ijRst, ijRnd, ijSst, ijSnd, ijSstdiag, ijSnddiag, ijSstsono, ijSndsono
      INTEGER  ::   iiout, ijout, iioutdir, ijoutdir, icnt
      INTEGER  ::   iRnei, iRdiag, iRsono
      INTEGER  ::   iSnei, iSdiag, iSsono                  !   -       -
      INTEGER  ::   iwe, ies, iso, ino, inum, id_dummy     !   -       -
      INTEGER  ::   jpbdta                                 !   -       -
      INTEGER  ::   ib_bdy1, ib_bdy2, ib1, ib2             !   -       -
      INTEGER  ::   ii1, ii2, ii3, ij1, ij2, ij3           !   -       -
      INTEGER  ::   iibe, ijbe, iibi, ijbi                 !   -       -
      INTEGER  ::   flagu, flagv                           ! short cuts
      INTEGER  ::   nbdyind, nbdybeg, nbdyend
      INTEGER  ::   itanh                                  ! reference scaling for FRS tanh profile (RDP)
      INTEGER              , DIMENSION(4)             ::   kdimsz
      INTEGER              , DIMENSION(jpbgrd,jp_bdy) ::   nblendta          ! Length of index arrays 
      INTEGER,  ALLOCATABLE, DIMENSION(:,:,:)         ::   nbidta, nbjdta    ! Index arrays: i and j indices of bdy dta
      INTEGER,  ALLOCATABLE, DIMENSION(:,:,:)         ::   nbrdta            ! Discrete distance from rim points
      CHARACTER(LEN=1)     , DIMENSION(jpbgrd)        ::   cgrid
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     ::   zz_read                 ! work space for 2D global boundary data
      REAL(wp), POINTER    , DIMENSION(:,:)     ::   zmask                   ! pointer to 2D mask fields
      REAL(wp)             , DIMENSION(jpi,jpj) ::   zfmask   ! temporary fmask array excluding coastal boundary condition (shlat)
      REAL(wp)             , DIMENSION(jpi,jpj) ::   ztmask, zumask, zvmask  ! temporary u/v mask array
      REAL(wp)             , DIMENSION(jpi,jpj) ::   zzbdy
      !!----------------------------------------------------------------------
      !
      cgrid = (/'t','u','v'/)

      ! -----------------------------------------
      ! Check and write out namelist parameters
      ! -----------------------------------------
      
      IF(lwp) WRITE(numout,*) 'Number of open boundary sets : ', nb_bdy

      DO ib_bdy = 1,nb_bdy

         IF(lwp) THEN
            WRITE(numout,*) ' ' 
            WRITE(numout,*) '------ Open boundary data set ',ib_bdy,' ------' 
            IF( ln_coords_file(ib_bdy) ) THEN
               WRITE(numout,*) 'Boundary definition read from file '//TRIM(cn_coords_file(ib_bdy))
            ELSE
               WRITE(numout,*) 'Boundary defined in namelist.'
            ENDIF
            WRITE(numout,*)
         ENDIF

         ! barotropic bdy
         !----------------
         IF(lwp) THEN
            WRITE(numout,*) 'Boundary conditions for barotropic solution:  '
            SELECT CASE( cn_dyn2d(ib_bdy) )                  
            CASE( 'none' )           ;   WRITE(numout,*) '      no open boundary condition'        
            CASE( 'frs' )            ;   WRITE(numout,*) '      Flow Relaxation Scheme'
            CASE( 'flather' )        ;   WRITE(numout,*) '      Flather radiation condition'
            CASE( 'orlanski' )       ;   WRITE(numout,*) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
            CASE( 'orlanski_npo' )   ;   WRITE(numout,*) '      Orlanski (NPO) radiation condition with adaptive nudging'
            CASE DEFAULT             ;   CALL ctl_stop( 'unrecognised value for cn_dyn2d' )
            END SELECT
         ENDIF

         dta_bdy(ib_bdy)%lneed_ssh   = cn_dyn2d(ib_bdy) == 'flather'
         dta_bdy(ib_bdy)%lneed_dyn2d = cn_dyn2d(ib_bdy) /= 'none'

         ! RDP override dta_bdy(ib_bdy)%ll_ssh with namelist value (ln_ssh_bdy)
         dta_bdy(ib_bdy)%lforced_ssh = ln_ssh_bdy(ib_bdy)
         IF(lwp) WRITE(numout,*) 'nambdy_ssh : use of ssh boundaries'
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         IF(lwp) WRITE(numout,*) '      ib_bdy: ',ib_bdy
         IF(lwp) WRITE(numout,*) '      dta_bdy(ib_bdy)%lneed_ssh  : ',dta_bdy(ib_bdy)%lneed_ssh
         IF(lwp) WRITE(numout,*) '      dta_bdy(ib_bdy)%lforced_ssh: ',dta_bdy(ib_bdy)%lforced_ssh
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         ! END RDP

         IF( lwp .AND. dta_bdy(ib_bdy)%lneed_dyn2d ) THEN
            SELECT CASE( nn_dyn2d_dta(ib_bdy) )                   ! 
            CASE( 0 )      ;   WRITE(numout,*) '      initial state used for bdy data'        
            CASE( 1 )      ;   WRITE(numout,*) '      boundary data taken from file'
            CASE( 2 )      ;   WRITE(numout,*) '      tidal harmonic forcing taken from file'
            CASE( 3 )      ;   WRITE(numout,*) '      boundary data AND tidal harmonic forcing taken from files'
            CASE DEFAULT   ;   CALL ctl_stop( 'nn_dyn2d_dta must be between 0 and 3' )
            END SELECT
         ENDIF
         IF ( dta_bdy(ib_bdy)%lneed_dyn2d .AND. nn_dyn2d_dta(ib_bdy) .GE. 2  .AND. .NOT.ln_tide ) THEN
            CALL ctl_stop( 'You must activate with ln_tide to add tidal forcing at open boundaries' )
         ENDIF
         IF(lwp) WRITE(numout,*)

         ! baroclinic bdy
         !----------------
         IF(lwp) THEN
            WRITE(numout,*) 'Boundary conditions for baroclinic velocities:  '
            SELECT CASE( cn_dyn3d(ib_bdy) )                  
            CASE('none')           ;   WRITE(numout,*) '      no open boundary condition'        
            CASE('frs')            ;   WRITE(numout,*) '      Flow Relaxation Scheme'
            CASE('specified')      ;   WRITE(numout,*) '      Specified value'
            CASE('neumann')        ;   WRITE(numout,*) '      Neumann conditions'
            CASE('zerograd')       ;   WRITE(numout,*) '      Zero gradient for baroclinic velocities'
            CASE('zero')           ;   WRITE(numout,*) '      Zero baroclinic velocities (runoff case)'
            CASE('orlanski')       ;   WRITE(numout,*) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
            CASE('orlanski_npo')   ;   WRITE(numout,*) '      Orlanski (NPO) radiation condition with adaptive nudging'
            CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for cn_dyn3d' )
            END SELECT
         ENDIF

         dta_bdy(ib_bdy)%lneed_dyn3d = cn_dyn3d(ib_bdy) == 'frs'      .OR. cn_dyn3d(ib_bdy) == 'specified'   &
            &                     .OR. cn_dyn3d(ib_bdy) == 'orlanski' .OR. cn_dyn3d(ib_bdy) == 'orlanski_npo'

         IF( lwp .AND. dta_bdy(ib_bdy)%lneed_dyn3d ) THEN
            SELECT CASE( nn_dyn3d_dta(ib_bdy) )                   ! 
            CASE( 0 )      ;   WRITE(numout,*) '      initial state used for bdy data'        
            CASE( 1 )      ;   WRITE(numout,*) '      boundary data taken from file'
            CASE DEFAULT   ;   CALL ctl_stop( 'nn_dyn3d_dta must be 0 or 1' )
            END SELECT
         END IF

         IF ( ln_dyn3d_dmp(ib_bdy) ) THEN
            IF ( cn_dyn3d(ib_bdy) == 'none' ) THEN
               IF(lwp) WRITE(numout,*) 'No open boundary condition for baroclinic velocities: ln_dyn3d_dmp is set to .false.'
               ln_dyn3d_dmp(ib_bdy) = .false.
            ELSEIF ( cn_dyn3d(ib_bdy) == 'frs' ) THEN
               CALL ctl_stop( 'Use FRS OR relaxation' )
            ELSE
               IF(lwp) WRITE(numout,*) '      + baroclinic velocities relaxation zone'
               IF(lwp) WRITE(numout,*) '      Damping time scale: ',rn_time_dmp(ib_bdy),' days'
               IF(rn_time_dmp(ib_bdy)<0) CALL ctl_stop( 'Time scale must be positive' )
               dta_bdy(ib_bdy)%lneed_dyn3d = .TRUE.
            ENDIF
         ELSE
            IF(lwp) WRITE(numout,*) '      NO relaxation on baroclinic velocities'
         ENDIF
         IF(lwp) WRITE(numout,*)

         !    tra bdy
         !----------------
         IF(lwp) THEN
            WRITE(numout,*) 'Boundary conditions for temperature and salinity:  '
            SELECT CASE( cn_tra(ib_bdy) )                  
            CASE('none')           ;   WRITE(numout,*) '      no open boundary condition'        
            CASE('frs')            ;   WRITE(numout,*) '      Flow Relaxation Scheme'
            CASE('specified')      ;   WRITE(numout,*) '      Specified value'
            CASE('neumann')        ;   WRITE(numout,*) '      Neumann conditions'
            CASE('runoff')         ;   WRITE(numout,*) '      Runoff conditions : Neumann for T and specified to 0.1 for salinity'
            CASE('orlanski')       ;   WRITE(numout,*) '      Orlanski (fully oblique) radiation condition with adaptive nudging'
            CASE('orlanski_npo')   ;   WRITE(numout,*) '      Orlanski (NPO) radiation condition with adaptive nudging'
            CASE DEFAULT           ;   CALL ctl_stop( 'unrecognised value for cn_tra' )
            END SELECT
         ENDIF

         dta_bdy(ib_bdy)%lneed_tra = cn_tra(ib_bdy) == 'frs'       .OR. cn_tra(ib_bdy) == 'specified'   &
            &                   .OR. cn_tra(ib_bdy) == 'orlanski'  .OR. cn_tra(ib_bdy) == 'orlanski_npo' 

         IF( lwp .AND. dta_bdy(ib_bdy)%lneed_tra ) THEN
            SELECT CASE( nn_tra_dta(ib_bdy) )                   ! 
            CASE( 0 )      ;   WRITE(numout,*) '      initial state used for bdy data'        
            CASE( 1 )      ;   WRITE(numout,*) '      boundary data taken from file'
            CASE DEFAULT   ;   CALL ctl_stop( 'nn_tra_dta must be 0 or 1' )
            END SELECT
         ENDIF

         IF ( ln_tra_dmp(ib_bdy) ) THEN
            IF ( cn_tra(ib_bdy) == 'none' ) THEN
               IF(lwp) WRITE(numout,*) 'No open boundary condition for tracers: ln_tra_dmp is set to .false.'
               ln_tra_dmp(ib_bdy) = .false.
            ELSEIF ( cn_tra(ib_bdy) == 'frs' ) THEN
               CALL ctl_stop( 'Use FRS OR relaxation' )
            ELSE
               IF(lwp) WRITE(numout,*) '      + T/S relaxation zone'
               IF(lwp) WRITE(numout,*) '      Damping time scale: ',rn_time_dmp(ib_bdy),' days'
               IF(lwp) WRITE(numout,*) '      Outflow damping time scale: ',rn_time_dmp_out(ib_bdy),' days'
               IF(lwp.AND.rn_time_dmp(ib_bdy)<0) CALL ctl_stop( 'Time scale must be positive' )
               dta_bdy(ib_bdy)%lneed_tra = .TRUE.
            ENDIF
         ELSE
            IF(lwp) WRITE(numout,*) '      NO T/S relaxation'
         ENDIF
         IF(lwp) WRITE(numout,*)

#if defined key_si3
         IF(lwp) THEN
            WRITE(numout,*) 'Boundary conditions for sea ice:  '
            SELECT CASE( cn_ice(ib_bdy) )                  
            CASE('none')   ;   WRITE(numout,*) '      no open boundary condition'        
            CASE('frs')    ;   WRITE(numout,*) '      Flow Relaxation Scheme'
            CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for cn_ice' )
            END SELECT
         ENDIF

         dta_bdy(ib_bdy)%lneed_ice = cn_ice(ib_bdy) /= 'none'

         IF( dta_bdy(ib_bdy)%lneed_ice .AND. nn_ice /= 2 ) THEN
            WRITE(ctmp1,*) 'bdy number ', ib_bdy,', needs ice model but nn_ice = ', nn_ice
            CALL ctl_stop( ctmp1 )
         ENDIF

         IF( lwp .AND. dta_bdy(ib_bdy)%lneed_ice ) THEN 
            SELECT CASE( nn_ice_dta(ib_bdy) )                   ! 
            CASE( 0 )      ;   WRITE(numout,*) '      initial state used for bdy data'        
            CASE( 1 )      ;   WRITE(numout,*) '      boundary data taken from file'
            CASE DEFAULT   ;   CALL ctl_stop( 'nn_ice_dta must be 0 or 1' )
            END SELECT
         ENDIF
#else
         dta_bdy(ib_bdy)%lneed_ice = .FALSE.
#endif
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      Width of relaxation zone = ', nn_rimwidth(ib_bdy)
         IF(lwp) WRITE(numout,*)
         !
      END DO   ! nb_bdy

      IF( lwp ) THEN
         IF( ln_vol ) THEN                     ! check volume conservation (nn_volctl value)
            WRITE(numout,*) 'Volume correction applied at open boundaries'
            WRITE(numout,*)
            SELECT CASE ( nn_volctl )
            CASE( 1 )      ;   WRITE(numout,*) '      The total volume will be constant'
            CASE( 0 )      ;   WRITE(numout,*) '      The total volume will vary according to the surface E-P flux'
            CASE DEFAULT   ;   CALL ctl_stop( 'nn_volctl must be 0 or 1' )
            END SELECT
            WRITE(numout,*)
            !
            ! sanity check if used with tides        
            IF( ln_tide ) THEN 
               WRITE(numout,*) ' The total volume correction is not working with tides. '
               WRITE(numout,*) ' Set ln_vol to .FALSE. '
               WRITE(numout,*) ' or '
               WRITE(numout,*) ' equilibriate your bdy input files '
               CALL ctl_stop( 'The total volume correction is not working with tides.' )
            END IF
         ELSE
            WRITE(numout,*) 'No volume correction applied at open boundaries'
            WRITE(numout,*)
         ENDIF
      ENDIF

      ! -------------------------------------------------
      ! Initialise indices arrays for open boundaries
      ! -------------------------------------------------

      nblendta(:,:) = 0
      nbdysege = 0
      nbdysegw = 0
      nbdysegn = 0
      nbdysegs = 0

      ! Define all boundaries 
      ! ---------------------
      DO ib_bdy = 1, nb_bdy
         !
         IF( .NOT. ln_coords_file(ib_bdy) ) THEN     ! build bdy coordinates with segments defined in namelist

            CALL bdy_read_seg( ib_bdy, nblendta(:,ib_bdy) )

         ELSE                                        ! Read size of arrays in boundary coordinates file.
            
            CALL iom_open( cn_coords_file(ib_bdy), inum )
            DO igrd = 1, jpbgrd
               id_dummy = iom_varid( inum, 'nbi'//cgrid(igrd), kdimsz=kdimsz )  
               nblendta(igrd,ib_bdy) = MAXVAL(kdimsz)
            END DO
            CALL iom_close( inum )
         ENDIF
         !
      END DO ! ib_bdy

      ! Now look for crossings in user (namelist) defined open boundary segments:
      IF( nbdysege > 0 .OR. nbdysegw > 0 .OR. nbdysegn > 0 .OR. nbdysegs > 0)   CALL bdy_ctl_seg
      
      ! Allocate arrays
      !---------------
      jpbdta = MAXVAL(nblendta(1:jpbgrd,1:nb_bdy))
      ALLOCATE( nbidta(jpbdta, jpbgrd, nb_bdy), nbjdta(jpbdta, jpbgrd, nb_bdy), nbrdta(jpbdta, jpbgrd, nb_bdy) )
      nbrdta(:,:,:) = 0   ! initialize nbrdta as it may not be completely defined for each bdy
      
      ! Calculate global boundary index arrays or read in from file
      !------------------------------------------------------------               
      ! 1. Read global index arrays from boundary coordinates file.
      DO ib_bdy = 1, nb_bdy
         !
         IF( ln_coords_file(ib_bdy) ) THEN
            !
            ALLOCATE( zz_read( MAXVAL(nblendta), 1 ) )          
            CALL iom_open( cn_coords_file(ib_bdy), inum )
            !
            DO igrd = 1, jpbgrd
               CALL iom_get( inum, jpdom_unknown, 'nbi'//cgrid(igrd), zz_read(1:nblendta(igrd,ib_bdy),:) )
               DO ii = 1,nblendta(igrd,ib_bdy)
                  nbidta(ii,igrd,ib_bdy) = NINT( zz_read(ii,1) ) + nn_hls
               END DO
               CALL iom_get( inum, jpdom_unknown, 'nbj'//cgrid(igrd), zz_read(1:nblendta(igrd,ib_bdy),:) )
               DO ii = 1,nblendta(igrd,ib_bdy)
                  nbjdta(ii,igrd,ib_bdy) = NINT( zz_read(ii,1) ) + nn_hls
               END DO
               CALL iom_get( inum, jpdom_unknown, 'nbr'//cgrid(igrd), zz_read(1:nblendta(igrd,ib_bdy),:) )
               DO ii = 1,nblendta(igrd,ib_bdy)
                  nbrdta(ii,igrd,ib_bdy) = NINT( zz_read(ii,1) )
               END DO
               !
               ibr_max = MAXVAL( nbrdta(:,igrd,ib_bdy) )
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' Maximum rimwidth in file is ', ibr_max
               IF(lwp) WRITE(numout,*) ' nn_rimwidth from namelist is ', nn_rimwidth(ib_bdy)
               IF (ibr_max < nn_rimwidth(ib_bdy))   &
                  CALL ctl_stop( 'nn_rimwidth is larger than maximum rimwidth in file',cn_coords_file(ib_bdy) )
            END DO
            !
            CALL iom_close( inum )
            DEALLOCATE( zz_read )
            !
         ENDIF
         !
      END DO

      ! 2. Now fill indices corresponding to straight open boundary arrays:
      CALL bdy_coords_seg( nbidta, nbjdta, nbrdta )

      !  Deal with duplicated points
      !-----------------------------
      ! We assign negative indices to duplicated points (to remove them from bdy points to be updated)
      ! if their distance to the bdy is greater than the other
      ! If their distance are the same, just keep only one to avoid updating a point twice
      DO igrd = 1, jpbgrd
         DO ib_bdy1 = 1, nb_bdy
            DO ib_bdy2 = 1, nb_bdy
               IF (ib_bdy1/=ib_bdy2) THEN
                  DO ib1 = 1, nblendta(igrd,ib_bdy1)
                     DO ib2 = 1, nblendta(igrd,ib_bdy2)
                        IF ((nbidta(ib1, igrd, ib_bdy1)==nbidta(ib2, igrd, ib_bdy2)).AND. &
                           &   (nbjdta(ib1, igrd, ib_bdy1)==nbjdta(ib2, igrd, ib_bdy2))) THEN
                           !                           IF ((lwp).AND.(igrd==1)) WRITE(numout,*) ' found coincident point ji, jj:', & 
                           !                                                       &              nbidta(ib1, igrd, ib_bdy1),      & 
                           !                                                       &              nbjdta(ib2, igrd, ib_bdy2)
                           ! keep only points with the lowest distance to boundary:
                           IF (nbrdta(ib1, igrd, ib_bdy1)<nbrdta(ib2, igrd, ib_bdy2)) THEN
                              nbidta(ib2, igrd, ib_bdy2) =-ib_bdy2
                              nbjdta(ib2, igrd, ib_bdy2) =-ib_bdy2
                           ELSEIF (nbrdta(ib1, igrd, ib_bdy1)>nbrdta(ib2, igrd, ib_bdy2)) THEN
                              nbidta(ib1, igrd, ib_bdy1) =-ib_bdy1
                              nbjdta(ib1, igrd, ib_bdy1) =-ib_bdy1
                              ! Arbitrary choice if distances are the same:
                           ELSE
                              nbidta(ib1, igrd, ib_bdy1) =-ib_bdy1
                              nbjdta(ib1, igrd, ib_bdy1) =-ib_bdy1
                           ENDIF
                        END IF
                     END DO
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      !
      ! Find lenght of boundaries and rim on local mpi domain
      !------------------------------------------------------
      !
      iwe = mig(1)
      ies = mig(jpi)
      iso = mjg(1) 
      ino = mjg(jpj) 
      !
      DO ib_bdy = 1, nb_bdy
         DO igrd = 1, jpbgrd
            icount   = 0   ! initialization of local bdy length
            icountr  = 0   ! initialization of local rim 0 and rim 1 bdy length
            icountr0 = 0   ! initialization of local rim 0 bdy length
            idx_bdy(ib_bdy)%nblen(igrd)     = 0
            idx_bdy(ib_bdy)%nblenrim(igrd)  = 0
            idx_bdy(ib_bdy)%nblenrim0(igrd) = 0
            DO ib = 1, nblendta(igrd,ib_bdy)
               ! check that data is in correct order in file
               IF( ib > 1 ) THEN
                  IF( nbrdta(ib,igrd,ib_bdy) < nbrdta(ib-1,igrd,ib_bdy) ) THEN
                     CALL ctl_stop('bdy_segs : ERROR : boundary data in file must be defined ', &
                        &        ' in order of distance from edge nbr A utility for re-ordering ', &
                        &        ' boundary coordinates and data files exists in the TOOLS/OBC directory')
                  ENDIF
               ENDIF
               ! check if point is in local domain
               IF(  nbidta(ib,igrd,ib_bdy) >= iwe .AND. nbidta(ib,igrd,ib_bdy) <= ies .AND.   &
                  & nbjdta(ib,igrd,ib_bdy) >= iso .AND. nbjdta(ib,igrd,ib_bdy) <= ino      ) THEN
                  !
                  icount = icount + 1
                  IF( nbrdta(ib,igrd,ib_bdy) == 1 .OR. nbrdta(ib,igrd,ib_bdy) == 0 )   icountr = icountr + 1
                  IF( nbrdta(ib,igrd,ib_bdy) == 0 )   icountr0 = icountr0 + 1
               ENDIF
            END DO
            idx_bdy(ib_bdy)%nblen    (igrd) = icount   !: length of boundary data on each proc
            idx_bdy(ib_bdy)%nblenrim (igrd) = icountr  !: length of rim 0 and rim 1 boundary data on each proc   
            idx_bdy(ib_bdy)%nblenrim0(igrd) = icountr0 !: length of rim 0 boundary data on each proc     
         END DO   ! igrd

         ! Allocate index arrays for this boundary set
         !--------------------------------------------
         ilen1 = MAXVAL( idx_bdy(ib_bdy)%nblen(:) )
         ALLOCATE( idx_bdy(ib_bdy)%nbi   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbj   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbr   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbd   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbdout(ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%ntreat(ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbmap (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%nbw   (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%flagu (ilen1,jpbgrd) ,   &
            &      idx_bdy(ib_bdy)%flagv (ilen1,jpbgrd) )

         ! Dispatch mapping indices and discrete distances on each processor
         ! -----------------------------------------------------------------
         DO igrd = 1, jpbgrd
            icount  = 0
            ! Outer loop on rimwidth to ensure outermost points come first in the local arrays.
            DO ir = 0, nn_rimwidth(ib_bdy)
               DO ib = 1, nblendta(igrd,ib_bdy)
                  ! check if point is in local domain and equals ir
                  IF(  nbidta(ib,igrd,ib_bdy) >= iwe .AND. nbidta(ib,igrd,ib_bdy) <= ies .AND.   &
                     & nbjdta(ib,igrd,ib_bdy) >= iso .AND. nbjdta(ib,igrd,ib_bdy) <= ino .AND.   &
                     & nbrdta(ib,igrd,ib_bdy) == ir  ) THEN
                     !
                     icount = icount  + 1
                     idx_bdy(ib_bdy)%nbi(icount,igrd)   = nbidta(ib,igrd,ib_bdy) - mig(1) + 1   ! global to local indexes
                     idx_bdy(ib_bdy)%nbj(icount,igrd)   = nbjdta(ib,igrd,ib_bdy) - mjg(1) + 1   ! global to local indexes
                     idx_bdy(ib_bdy)%nbr(icount,igrd)   = nbrdta(ib,igrd,ib_bdy)
                     idx_bdy(ib_bdy)%nbmap(icount,igrd) = ib
                  ENDIF
               END DO
            END DO
         END DO   ! igrd

      END DO   ! ib_bdy

      ! Initialize array indicating communications in bdy
      ! -------------------------------------------------
      ALLOCATE( lsend_bdyolr(nb_bdy,jpbgrd,8,0:1), lrecv_bdyolr(nb_bdy,jpbgrd,8,0:1) )
      lsend_bdyolr(:,:,:,:) = .false.
      lrecv_bdyolr(:,:,:,:) = .false. 

      DO ib_bdy = 1, nb_bdy
         DO igrd = 1, jpbgrd
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)   ! only the rim triggers communications, see bdy routines
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               IF( ib .LE. idx_bdy(ib_bdy)%nblenrim0(igrd) ) THEN   ;   ir = 0
               ELSE                                                 ;   ir = 1
               END IF
               !
               ! check if point has to be sent     to   a neighbour
               IF( ii >= Nis0 .AND. ii < Nis0 + nn_hls .AND. ij >= Njs0 .AND. ij <= Nje0         ) THEN   ! we inner side
                  IF( mpiSnei(nn_hls,jpwe) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpwe,ir) = .TRUE.
               ENDIF
               IF( ii <= Nie0 .AND. ii > Nie0 - nn_hls .AND. ij >= Njs0 .AND. ij <= Nje0         ) THEN   ! ea inner side
                  IF( mpiSnei(nn_hls,jpea) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpea,ir) = .TRUE.
               ENDIF
               IF( ii >= Nis0 .AND. ii <= Nie0         .AND. ij >= Njs0 .AND. ij < Njs0 + nn_hls ) THEN   ! so inner side
                  IF( mpiSnei(nn_hls,jpso) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpso,ir) = .TRUE.
               ENDIF
               IF( ii  < Nis0                          .AND. ij >= Njs0 .AND. ij < Njs0 + nn_hls ) THEN   ! so side we-halo
                  IF( mpiSnei(nn_hls,jpso) > -1 .AND. nn_comm == 1 )   lsend_bdyolr(ib_bdy,igrd,jpso,ir) = .TRUE.
               ENDIF
               IF( ii  > Nie0                          .AND. ij >= Njs0 .AND. ij < Njs0 + nn_hls ) THEN   ! so side ea-halo 
                  IF( mpiSnei(nn_hls,jpso) > -1 .AND. nn_comm == 1 )   lsend_bdyolr(ib_bdy,igrd,jpso,ir) = .TRUE.
               ENDIF
               IF( ii >= Nis0 .AND. ii <= Nie0         .AND. ij <= Nje0 .AND. ij > Nje0 - nn_hls ) THEN   ! no inner side
                  IF( mpiSnei(nn_hls,jpno) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpno,ir) = .TRUE.
               ENDIF
               IF( ii  < Nis0                          .AND. ij <= Nje0 .AND. ij > Nje0 - nn_hls ) THEN   ! no side we-halo
                  IF( mpiSnei(nn_hls,jpno) > -1 .AND. nn_comm == 1 )   lsend_bdyolr(ib_bdy,igrd,jpno,ir) = .TRUE.
               ENDIF
               IF( ii  > Nie0                          .AND. ij <= Nje0 .AND. ij > Nje0 - nn_hls ) THEN   ! no side ea-halo
                  IF( mpiSnei(nn_hls,jpno) > -1 .AND. nn_comm == 1 )   lsend_bdyolr(ib_bdy,igrd,jpno,ir) = .TRUE.
               ENDIF
               IF( ii >= Nis0 .AND. ii < Nis0 + nn_hls .AND. ij >= Njs0 .AND. ij < Njs0 + nn_hls ) THEN   ! sw inner corner
                  IF( mpiSnei(nn_hls,jpsw) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpsw,ir) = .TRUE.
               ENDIF
               IF( ii <= Nie0 .AND. ii > Nie0 - nn_hls .AND. ij >= Njs0 .AND. ij < Njs0 + nn_hls ) THEN   ! se inner corner
                  IF( mpiSnei(nn_hls,jpse) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpse,ir) = .TRUE.
               ENDIF
               IF( ii >= Nis0 .AND. ii < Nis0 + nn_hls .AND. ij <= Nje0 .AND. ij > Nje0 - nn_hls ) THEN   ! nw inner corner
                  IF( mpiSnei(nn_hls,jpnw) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpnw,ir) = .TRUE.
               ENDIF
               IF( ii <= Nie0 .AND. ii > Nie0 - nn_hls .AND. ij <= Nje0 .AND. ij > Nje0 - nn_hls ) THEN   ! ne inner corner
                  IF( mpiSnei(nn_hls,jpne) > -1                    )   lsend_bdyolr(ib_bdy,igrd,jpne,ir) = .TRUE.
               ENDIF
               !
               ! check if point has to be received from a neighbour
               IF( ii  < Nis0                  .AND. ij >= Njs0 .AND. ij <= Nje0 ) THEN   ! we side
                  IF( mpiRnei(nn_hls,jpwe) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpwe,ir) = .TRUE.
               ENDIF
               IF( ii  > Nie0                  .AND. ij >= Njs0 .AND. ij <= Nje0 ) THEN   ! ea side
                  IF( mpiRnei(nn_hls,jpea) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpea,ir) = .TRUE.
               ENDIF
               IF( ii >= Nis0 .AND. ii <= Nie0 .AND. ij  < Njs0                  ) THEN   ! so side
                  IF( mpiRnei(nn_hls,jpso) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpso,ir) = .TRUE.
               ENDIF
               IF( ii >= Nis0 .AND. ii <= Nie0 .AND. ij  > Nje0                  ) THEN   ! no side
                  IF( mpiRnei(nn_hls,jpno) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpno,ir) = .TRUE.
               ENDIF
               IF( ii  < Nis0                  .AND. ij  < Njs0                  ) THEN   ! sw corner
                  IF( mpiRnei(nn_hls,jpsw) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpsw,ir) = .TRUE.
                  IF( mpiRnei(nn_hls,jpso) > -1 .AND. nn_comm == 1 )   lrecv_bdyolr(ib_bdy,igrd,jpso,ir) = .TRUE.
               ENDIF
               IF( ii  > Nie0                  .AND. ij  < Njs0                  ) THEN   ! se corner
                  IF( mpiRnei(nn_hls,jpse) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpse,ir) = .TRUE.
                  IF( mpiRnei(nn_hls,jpso) > -1 .AND. nn_comm == 1 )   lrecv_bdyolr(ib_bdy,igrd,jpso,ir) = .TRUE.
               ENDIF
               IF( ii  < Nis0                  .AND. ij  > Nje0                  ) THEN   ! nw corner
                  IF( mpiRnei(nn_hls,jpnw) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpnw,ir) = .TRUE.
                  IF( mpiRnei(nn_hls,jpno) > -1 .AND. nn_comm == 1 )   lrecv_bdyolr(ib_bdy,igrd,jpno,ir) = .TRUE.
               ENDIF
               IF( ii  > Nie0                  .AND. ij  > Nje0                  ) THEN   ! ne corner
                  IF( mpiRnei(nn_hls,jpne) > -1                    )   lrecv_bdyolr(ib_bdy,igrd,jpne,ir) = .TRUE.
                  IF( mpiRnei(nn_hls,jpno) > -1 .AND. nn_comm == 1 )   lrecv_bdyolr(ib_bdy,igrd,jpno,ir) = .TRUE.
               ENDIF
               !
            END DO
         END DO   !   igrd
         
         ! Comment out for debug
!!$         DO ir = 0,1
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'T', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyolr(ib_bdy,1,:,ir), lrecv = lrecv_bdyolr(ib_bdy,1,:,ir) )
!!$            IF(lwp) WRITE(numout,*) ' seb bdy debug olr T', ir ; CALL FLUSH(numout)
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'U', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyolr(ib_bdy,2,:,ir), lrecv = lrecv_bdyolr(ib_bdy,2,:,ir) )
!!$            IF(lwp) WRITE(numout,*) ' seb bdy debug olr U', ir ; CALL FLUSH(numout)
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'V', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyolr(ib_bdy,3,:,ir), lrecv = lrecv_bdyolr(ib_bdy,3,:,ir) )    
!!$            IF(lwp) WRITE(numout,*) ' seb bdy debug olr V', ir ; CALL FLUSH(numout)
!!$         END DO
         
         ! Compute rim weights for FRS scheme
         ! ----------------------------------
         itanh = 10 ! RDP reference length scale for TanH profile (from JG)
         DO igrd = 1, jpbgrd
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ir = MAX( 1, idx_bdy(ib_bdy)%nbr(ib,igrd) )   ! both rim 0 and rim 1 have the same weights
               ! RDP Set TanH profile equivalent regardless of rimwidth, according to reference length scale (itanh)
               idx_bdy(ib_bdy)%nbw(ib,igrd) = 1.- TANH( REAL( ir - 1 ) *0.5 &
               &          *(FLOAT(itanh)/FLOAT(nn_rimwidth(ib_bdy))) )      ! tanh formulation
               ! END RDP
               !               idx_bdy(ib_bdy)%nbw(ib,igrd) = (REAL(nn_rimwidth(ib_bdy)+1-ir)/REAL(nn_rimwidth(ib_bdy)))**2.  ! quadratic
               !               idx_bdy(ib_bdy)%nbw(ib,igrd) =  REAL(nn_rimwidth(ib_bdy)+1-ir)/REAL(nn_rimwidth(ib_bdy))       ! linear
            END DO
         END DO

         ! Compute damping coefficients
         ! ----------------------------
         DO igrd = 1, jpbgrd
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ir = MAX( 1, idx_bdy(ib_bdy)%nbr(ib,igrd) )   ! both rim 0 and rim 1 have the same damping coefficients
               idx_bdy(ib_bdy)%nbd(ib,igrd) = 1. / ( rn_time_dmp(ib_bdy) * rday ) & 
                  & *(REAL(nn_rimwidth(ib_bdy)+1-ir)/REAL(nn_rimwidth(ib_bdy)))**2.   ! quadratic
               idx_bdy(ib_bdy)%nbdout(ib,igrd) = 1. / ( rn_time_dmp_out(ib_bdy) * rday ) & 
                  & *(REAL(nn_rimwidth(ib_bdy)+1-ir)/REAL(nn_rimwidth(ib_bdy)))**2.   ! quadratic
            END DO
         END DO

      END DO ! ib_bdy

      ! ------------------------------------------------------
      ! Initialise masks and find normal/tangential directions
      ! ------------------------------------------------------

      ! ------------------------------------------
      ! handle rim0, do as if rim 1 was free ocean
      ! ------------------------------------------

      ztmask(:,:) = tmask(:,:,1)   ;   zumask(:,:) = umask(:,:,1)   ;   zvmask(:,:) = vmask(:,:,1)
      ! For the flagu/flagv calculation below we require a version of fmask without
      ! the land boundary condition (shlat) included:
      DO_2D( 0, 0, 0, 0 )
         zfmask(ji,jj) =  ztmask(ji,jj  ) * ztmask(ji+1,jj  )   &
            &           * ztmask(ji,jj+1) * ztmask(ji+1,jj+1)
      END_2D
      CALL lbc_lnk( 'bdyini', zfmask, 'F', 1.0_wp )

      ! Read global 2D mask at T-points: bdytmask
      ! -----------------------------------------
      ! bdytmask = 1  on the computational domain but not on open boundaries
      !          = 0  elsewhere   

      bdytmask(:,:) = ssmask(:,:)

      ! Derive mask on U and V grid from mask on T grid
      DO_2D( 0, 0, 0, 0 )
            bdyumask(ji,jj) = bdytmask(ji,jj) * bdytmask(ji+1,jj  )
            bdyvmask(ji,jj) = bdytmask(ji,jj) * bdytmask(ji  ,jj+1)  
      END_2D
      CALL lbc_lnk( 'bdyini', bdyumask, 'U', 1.0_wp , bdyvmask, 'V', 1.0_wp )   ! Lateral boundary cond. 

      ! bdy masks are now set to zero on rim 0 points:
      DO ib_bdy = 1, nb_bdy
         DO ib = 1, idx_bdy(ib_bdy)%nblenrim0(1)   ! extent of rim 0
            bdytmask(idx_bdy(ib_bdy)%nbi(ib,1), idx_bdy(ib_bdy)%nbj(ib,1)) = 0._wp
         END DO
         DO ib = 1, idx_bdy(ib_bdy)%nblenrim0(2)   ! extent of rim 0
            bdyumask(idx_bdy(ib_bdy)%nbi(ib,2), idx_bdy(ib_bdy)%nbj(ib,2)) = 0._wp
         END DO
         DO ib = 1, idx_bdy(ib_bdy)%nblenrim0(3)   ! extent of rim 0
            bdyvmask(idx_bdy(ib_bdy)%nbi(ib,3), idx_bdy(ib_bdy)%nbj(ib,3)) = 0._wp
         END DO
      END DO

      CALL bdy_rim_treat( zumask, zvmask, zfmask, .true. )   ! compute flagu, flagv, ntreat on rim 0

      ! ------------------------------------
      ! handle rim1, do as if rim 0 was land
      ! ------------------------------------
      
      ! z[tuv]mask are now set to zero on rim 0 points:
      DO ib_bdy = 1, nb_bdy
         DO ib = 1, idx_bdy(ib_bdy)%nblenrim0(1)   ! extent of rim 0
            ztmask(idx_bdy(ib_bdy)%nbi(ib,1), idx_bdy(ib_bdy)%nbj(ib,1)) = 0._wp
         END DO
         DO ib = 1, idx_bdy(ib_bdy)%nblenrim0(2)   ! extent of rim 0
            zumask(idx_bdy(ib_bdy)%nbi(ib,2), idx_bdy(ib_bdy)%nbj(ib,2)) = 0._wp
         END DO
         DO ib = 1, idx_bdy(ib_bdy)%nblenrim0(3)   ! extent of rim 0
            zvmask(idx_bdy(ib_bdy)%nbi(ib,3), idx_bdy(ib_bdy)%nbj(ib,3)) = 0._wp
         END DO
      END DO

      ! Recompute zfmask
      DO_2D( 0, 0, 0, 0 )
         zfmask(ji,jj) =  ztmask(ji,jj  ) * ztmask(ji+1,jj  )   &
            &           * ztmask(ji,jj+1) * ztmask(ji+1,jj+1)
      END_2D
      CALL lbc_lnk( 'bdyini', zfmask, 'F', 1.0_wp )

      ! bdy masks are now set to zero on rim1 points:
      DO ib_bdy = 1, nb_bdy
         DO ib = idx_bdy(ib_bdy)%nblenrim0(1) + 1,  idx_bdy(ib_bdy)%nblenrim(1)   ! extent of rim 1
            bdytmask(idx_bdy(ib_bdy)%nbi(ib,1), idx_bdy(ib_bdy)%nbj(ib,1)) = 0._wp
         END DO
         DO ib = idx_bdy(ib_bdy)%nblenrim0(2) + 1,  idx_bdy(ib_bdy)%nblenrim(2)   ! extent of rim 1
            bdyumask(idx_bdy(ib_bdy)%nbi(ib,2), idx_bdy(ib_bdy)%nbj(ib,2)) = 0._wp
         END DO
         DO ib = idx_bdy(ib_bdy)%nblenrim0(3) + 1,  idx_bdy(ib_bdy)%nblenrim(3)   ! extent of rim 1
            bdyvmask(idx_bdy(ib_bdy)%nbi(ib,3), idx_bdy(ib_bdy)%nbj(ib,3)) = 0._wp
         END DO
      END DO

      CALL bdy_rim_treat( zumask, zvmask, zfmask, .false. )   ! compute flagu, flagv, ntreat on rim 1
      !
      ! Check which boundaries might need communication
      ALLOCATE( lsend_bdyint(nb_bdy,jpbgrd,8,0:1), lrecv_bdyint(nb_bdy,jpbgrd,8,0:1) )
      lsend_bdyint(:,:,:,:) = .false.
      lrecv_bdyint(:,:,:,:) = .false. 
      ALLOCATE( lsend_bdyext(nb_bdy,jpbgrd,8,0:1), lrecv_bdyext(nb_bdy,jpbgrd,8,0:1) )
      lsend_bdyext(:,:,:,:) = .false.
      lrecv_bdyext(:,:,:,:) = .false.
      !
      DO ib_bdy = 1, nb_bdy
         DO igrd = 1, jpbgrd
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
               IF( idx_bdy(ib_bdy)%ntreat(ib,igrd) == -1 ) CYCLE
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               ir = idx_bdy(ib_bdy)%nbr(ib,igrd)
               flagu = NINT(idx_bdy(ib_bdy)%flagu(ib,igrd))
               flagv = NINT(idx_bdy(ib_bdy)%flagv(ib,igrd))
               iibe = ii - flagu   ! neighbouring point towards the exterior of the computational domain
               ijbe = ij - flagv
               iibi = ii + flagu   ! neighbouring point towards the interior of the computational domain
               ijbi = ij + flagv
               CALL find_neib( ii, ij, idx_bdy(ib_bdy)%ntreat(ib,igrd), ii1, ij1, ii2, ij2, ii3, ij3 )   ! free ocean neighbours
               !
               !  take care of the 4 sides
               !
               DO icnt = 1, 4
                  SELECT CASE( icnt )
                     !                                           ... _____
                  CASE( 1 )   ! x: rim on rcvwe/sndea-side         o|  :
                     !          o: potential neighbour(s)          o|x :
                     !             outside of the MPI domain     ..o|__:__
                     iRnei    = jpwe             ;   iSnei    = jpea
                     iiRst    = 1                ;   ijRst    = Njs0            ! Rcv we-side starting point, excluding sw-corner
                     iiRnd    = nn_hls           ;   ijRnd    = Nje0            ! Rcv we-side   ending point, excluding nw-corner
                     iiSst    = Nie0-nn_hls+1    ;   ijSst    = Njs0            ! Snd ea-side starting point, excluding se-corner
                     iiSnd    = Nie0             ;   ijSnd    = Nje0            ! Snd ea-side   ending point, excluding ne-corner
                     iioutdir = -1               ;   ijoutdir = -999            ! outside MPI domain: westward
                     !                                           ______....
                  CASE( 2 )   ! x: rim on rcvea/sndwe-side          :  |o
                     !          o: potential neighbour(s)           : x|o
                     !             outside of the MPI domain     ___:__|o..
                     iRnei    = jpea             ;   iSnei    = jpwe
                     iiRst    = Nie0+1           ;   ijRst    =  Njs0            ! Rcv ea-side starting point, excluding se-corner
                     iiRnd    = jpi              ;   ijRnd    =  Nje0            ! Rcv ea-side   ending point, excluding ne-corner
                     iiSst    = Nis0             ;   ijSst    =  Njs0            ! Snd we-side starting point, excluding sw-corner
                     iiSnd    = Nis0+nn_hls-1    ;   ijSnd    =  Nje0            ! Snd we-side   ending point, excluding nw-corner
                     iioutdir = 1                ;   ijoutdir = -999             ! outside MPI domain: eastward
                     !
                  CASE( 3 )   ! x: rim on rcvso/sndno-side       |       |
                     !          o: potential neighbour(s)        |¨¨¨¨¨¨¨|
                     !             outside of the MPI domain     |___x___|
                     !                                           : o o o :
                     !                                           :       :
                     iRnei    = jpso             ;   iSnei    = jpno
                     iiRst    = Nis0             ;   ijRst    = 1                ! Rcv so-side starting point, excluding sw-corner
                     iiRnd    = Nie0             ;   ijRnd    = nn_hls           ! Rcv so-side   ending point, excluding se-corner
                     iiSst    = Nis0             ;   ijSst    = Nje0-nn_hls+1    ! Snd no-side starting point, excluding nw-corner
                     iiSnd    = Nie0             ;   ijSnd    = Nje0             ! Snd no-side   ending point, excluding ne-corner
                     iioutdir = -999             ;   ijoutdir = -1               ! outside MPI domain: southward
                     !                                           :       :
                  CASE( 4 )   ! x: rim on rcvno/sndso-side       :_o_o_o_:
                     !          o: potential neighbour(s)        |   x   |
                     !             outside of the MPI domain     |       |
                     !                                           |¨¨¨¨¨¨¨|
                     iRnei    = jpno             ;   iSnei    = jpso
                     iiRst    = Nis0             ;   ijRst    = Nje0+1           ! Rcv no-side starting point, excluding nw-corner
                     iiRnd    = Nie0             ;   ijRnd    = jpj              ! Rcv no-side   ending point, excluding ne-corner
                     iiSst    = Nis0             ;   ijSst    = Njs0             ! Snd so-side starting point, excluding sw-corner
                     iiSnd    = Nie0             ;   ijSnd    = Njs0+nn_hls-1    ! Snd so-side   ending point, excluding se-corner
                     iioutdir = -999             ;   ijoutdir = 1                ! outside MPI domain: northward
                  END SELECT
                  !
                  IF( ii >= iiRst .AND. ii <= iiRnd .AND. ij >= ijRst .AND. ij <= ijRnd ) THEN   ! rim point in recv side
                     iiout = ii+iioutdir ; ijout = ij+ijoutdir        ! in which direction do we go outside of the MPI domain?
                     ! take care of neighbourg(s) in the interior of the computational domain
                     IF(  iibi==iiout .OR. ii1==iiout .OR. ii2==iiout .OR. ii3==iiout .OR.   &   ! Neib outside of the MPI domain
                        & ijbi==ijout .OR. ij1==ijout .OR. ij2==ijout .OR. ij3==ijout ) THEN     ! -> I cannot compute it -> recv it
                        IF( mpiRnei(nn_hls,iRnei) > -1 )   lrecv_bdyint(ib_bdy,igrd,iRnei,ir) = .TRUE.
                     ENDIF
                     ! take care of neighbourg in the exterior of the computational domain
                     IF(  iibe==iiout .OR. ijbe==ijout ) THEN   ! Neib outside of the MPI domain -> I cannot compute it -> recv it
                        IF( mpiRnei(nn_hls,iRnei) > -1 )   lrecv_bdyext(ib_bdy,igrd,iRnei,ir) = .TRUE.
                     ENDIF
                  ENDIF
                  
                  IF( ii >= iiSst .AND. ii <= iiSnd .AND. ij >= ijSst .AND. ij <= ijSnd ) THEN   ! rim point in send side
                     iiout = ii+iioutdir ; ijout = ij+ijoutdir        ! in which direction do we go outside of the nei MPI domain?
                     ! take care of neighbourg(s) in the interior of the computational domain
                     IF(  iibi==iiout .OR. ii1==iiout .OR. ii2==iiout .OR. ii3==iiout .OR.   &   ! Neib outside of nei MPI domain
                        & ijbi==ijout .OR. ij1==ijout .OR. ij2==ijout .OR. ij3==ijout ) THEN     ! -> nei cannot compute it
                        IF( mpiSnei(nn_hls,iSnei) > -1 )   lsend_bdyint(ib_bdy,igrd,iSnei,ir) = .TRUE.   ! -> send to nei
                     ENDIF
                     ! take care of neighbourg in the exterior of the computational domain
                     IF( iibe == iiout .OR. ijbe == ijout ) THEN   ! Neib outside of the nei MPI domain -> nei cannot compute it
                        IF( mpiSnei(nn_hls,iSnei) > -1 )   lsend_bdyext(ib_bdy,igrd,iSnei,ir) = .TRUE.   ! -> send to nei
                     ENDIF
                  END IF

               END DO   ! 4 sides
               !
               ! specific treatment for the corners
               !
               DO icnt = 1, 4
                  SELECT CASE( icnt )
                     !                                       ...|....
                  CASE( 1 )   ! x: rim on sw-corner            o|   :
                     !          o: potential neighbour(s)      o|x__:__
                     !             outside of the MPI domain   o o o:
                     !                                              :
                     iRdiag    = jpsw            ;   iRsono    = jpso            ! Recv: for sw or so
                     iSdiag    = jpne            ;   iSsono    = jpno            ! Send: to ne or no
                     iiRst     = 1               ;   ijRst     = 1               ! Rcv sw-corner starting point
                     iiRnd     = nn_hls          ;   ijRnd     = nn_hls          ! Rcv sw-corner   ending point
                     iiSstdiag = Nie0-nn_hls+1   ;   ijSstdiag = Nje0-nn_hls+1   ! send to sw-corner of ne neighbourg
                     iiSnddiag = Nie0            ;   ijSnddiag = Nje0            ! send to sw-corner of ne neighbourg
                     iiSstsono = 1               ;   ijSstsono = Nje0-nn_hls+1   ! send to sw-corner of no neighbourg
                     iiSndsono = nn_hls          ;   ijSndsono = Nje0            ! send to sw-corner of no neighbourg
                     iioutdir  = -1              ;   ijoutdir  = -1              ! outside MPI domain: westward or southward
                     !                                          ....|...
                  CASE( 2 )   ! x: rim on se-corner             :   |o
                     !          o: potential neighbour(s)     __:__x|o
                     !             outside of the MPI domain    :o o o
                     !                                          :    
                     iRdiag    = jpse            ;   iRsono    = jpso            ! Recv: for se or so
                     iSdiag    = jpnw            ;   iSsono    = jpno            ! Send: to nw or no
                     iiRst     = Nie0+1          ;   ijRst     = 1               ! Rcv se-corner starting point
                     iiRnd     = jpi             ;   ijRnd     = nn_hls          ! Rcv se-corner   ending point
                     iiSstdiag = Nis0            ;   ijSstdiag = Nje0-nn_hls+1   ! send to se-corner of nw neighbourg
                     iiSnddiag = Nis0+nn_hls-1   ;   ijSnddiag = Nje0            ! send to se-corner of nw neighbourg
                     iiSstsono = Nie0+1          ;   ijSstsono = Nje0-nn_hls+1   ! send to se-corner of no neighbourg
                     iiSndsono = jpi             ;   ijSndsono = Nje0            ! send to se-corner of no neighbourg
                     iioutdir  = 1               ;   ijoutdir  = -1              ! outside MPI domain: eastward or southward
                     !                                              :       
                     !                                         o o_o:___
                  CASE( 3 )   ! x: rim on nw-corner            o|x  :
                     !          o: potential neighbour(s)    ..o|...:
                     !             outside of the MPI domain    |
                     iRdiag    = jpnw            ;   iRsono    = jpno            ! Recv: for nw or no
                     iSdiag    = jpse            ;   iSsono    = jpso            ! Send: to se or so
                     iiRst     = 1               ;   ijRst     = Nje0+1          ! Rcv nw-corner starting point
                     iiRnd     = nn_hls          ;   ijRnd     = jpj             ! Rcv nw-corner   ending point
                     iiSstdiag = Nie0-nn_hls+1   ;   ijSstdiag = Njs0            ! send to nw-corner of se neighbourg
                     iiSnddiag = Nie0            ;   ijSnddiag = Njs0+nn_hls-1   ! send to nw-corner of se neighbourg
                     iiSstsono = 1               ;   ijSstsono = Njs0            ! send to nw-corner of so neighbourg
                     iiSndsono = nn_hls          ;   ijSndsono = Njs0+nn_hls-1   ! send to nw-corner of so neighbourg
                     iioutdir  = -1              ;   ijoutdir  =  1              ! outside MPI domain: westward or northward
                     !                                          :       
                     !                                       ___:o_o o
                  CASE( 4 )   ! x: rim on ne-corner             :  x|o
                     !          o: potential neighbour(s)       :...|o...
                     !             outside of the MPI domain        |
                     iRdiag    = jpne            ;   iRsono    = jpno            ! Recv: for ne or no
                     iSdiag    = jpsw            ;   iSsono    = jpso            ! Send: to sw or so
                     iiRst     = Nie0+1          ;   ijRst     = Nje0+1          ! Rcv ne-corner starting point
                     iiRnd     = jpi             ;   ijRnd     = jpj             ! Rcv ne-corner   ending point
                     iiSstdiag = Nis0            ;   ijSstdiag = Njs0            ! send to ne-corner of sw neighbourg
                     iiSnddiag = Nis0+nn_hls-1   ;   ijSnddiag = Njs0+nn_hls-1   ! send to ne-corner of sw neighbourg
                     iiSstsono = Nie0+1          ;   ijSstsono = Njs0            ! send to ne-corner of so neighbourg
                     iiSndsono = jpi             ;   ijSndsono = Njs0+nn_hls-1   ! send to ne-corner of so neighbourg
                     iioutdir  = 1               ;   ijoutdir  = 1               ! outside MPI domain: eastward or southward
                  END SELECT
                  !
                  ! Check if we need to receive data for this rim point
                  IF( ii >= iiRst .AND. ii <= iiRnd .AND. ij >= ijRst .AND. ij <= ijRnd ) THEN   ! rim point on the corner
                     iiout = ii+iioutdir ; ijout = ij+ijoutdir        ! in which direction do we go outside of the MPI domain?
                     ! take care of neighbourg(s) in the interior of the computational domain
                     IF(  iibi==iiout .OR. ii1==iiout .OR. ii2==iiout .OR. ii3==iiout .OR.   &   ! Neib outside of the MPI domain
                        & ijbi==ijout .OR. ij1==ijout .OR. ij2==ijout .OR. ij3==ijout ) THEN     ! -> I cannot compute it -> recv it
                        IF( mpiRnei(nn_hls,iRdiag) > -1                    )   lrecv_bdyint(ib_bdy,igrd,iRdiag,ir) = .TRUE.   ! Receive directly from diagonal neighbourg
                        IF( mpiRnei(nn_hls,iRsono) > -1 .AND. nn_comm == 1 )   lrecv_bdyint(ib_bdy,igrd,iRsono,ir) = .TRUE.   ! Receive through the South/North neighbourg
                     ENDIF
                     ! take care of neighbourg in the exterior of the computational domain
                     IF(  iibe==iiout .OR. ijbe==ijout ) THEN   ! Neib outside of the MPI domain -> I cannot compute it -> recv it
                        IF( mpiRnei(nn_hls,iRdiag) > -1                    )   lrecv_bdyext(ib_bdy,igrd,iRdiag,ir) = .TRUE.   ! Receive directly from diagonal neighbourg
                        IF( mpiRnei(nn_hls,iRsono) > -1 .AND. nn_comm == 1 )   lrecv_bdyext(ib_bdy,igrd,iRsono,ir) = .TRUE.   ! Receive through the South/North neighbourg
                     ENDIF
                  ENDIF
                  !
                  ! Check if this rim point corresponds to the corner of one neighbourg. if yes, do we need to send data?
                  ! Direct send to diag: Is this rim point the corner point of a diag neighbour with which we communicate?
                  IF( ii >= iiSstdiag .AND. ii <= iiSnddiag .AND. ij >= ijSstdiag .AND. ij <= ijSnddiag   &
                     &                .AND. mpiSnei(nn_hls,iSdiag) > -1 ) THEN
                     iiout = ii+iioutdir ; ijout = ij+ijoutdir        ! in which direction do we go outside of the nei MPI domain?
                     ! take care of neighbourg(s) in the interior of the computational domain
                     IF(  iibi==iiout .OR. ii1==iiout .OR. ii2==iiout .OR. ii3==iiout .OR.   &   ! Neib outside of diag nei MPI 
                        & ijbi==ijout .OR. ij1==ijout .OR. ij2==ijout .OR. ij3==ijout )      &   ! domain -> nei cannot compute it
                        &    lsend_bdyint(ib_bdy,igrd,iSdiag,ir) = .TRUE.                        ! send rim point data to diag nei
                     ! take care of neighbourg in the exterior of the computational domain
                     IF(  iibe==iiout .OR. ijbe==ijout )   &                                 
                        &    lsend_bdyext(ib_bdy,igrd,iSdiag,ir) = .TRUE.
                  ENDIF
                  ! Indirect send to diag (through so/no): rim point is the corner point of a so/no nei with which we communicate
                  IF( ii >= iiSstsono .AND. ii <= iiSndsono .AND. ij >= ijSstsono .AND. ij <= ijSndsono   &
                     &                .AND. mpiSnei(nn_hls,iSsono) > -1 .AND. nn_comm == 1 ) THEN
                     iiout = ii+iioutdir ; ijout = ij+ijoutdir        ! in which direction do we go outside of the nei MPI domain?
                     ! take care of neighbourg(s) in the interior of the computational domain
                     IF(  iibi==iiout .OR. ii1==iiout .OR. ii2==iiout .OR. ii3==iiout .OR.   &   ! Neib outside of so/no nei MPI
                        & ijbi==ijout .OR. ij1==ijout .OR. ij2==ijout .OR. ij3==ijout )      &   ! domain -> nei cannot compute it
                        &    lsend_bdyint(ib_bdy,igrd,iSsono,ir) = .TRUE.                        ! send rim point data to so/no nei
                     ! take care of neighbourg in the exterior of the computational domain
                     IF(  iibe==iiout .OR. ijbe==ijout )   &
                        &    lsend_bdyext(ib_bdy,igrd,iSsono,ir) = .TRUE.
                  ENDIF
                  !
               END DO   ! 4 corners
            END DO   ! ib
         END DO   ! igrd

         ! Comment out for debug
!!$         DO ir = 0,1
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'T', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyint(ib_bdy,1,:,ir), lrecv = lrecv_bdyint(ib_bdy,1,:,ir) )
!!$            IF(lwp) WRITE(numout,*) ' bdy debug int T', ir ; CALL FLUSH(numout)
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'U', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyint(ib_bdy,2,:,ir), lrecv = lrecv_bdyint(ib_bdy,2,:,ir) )
!!$            IF(lwp) WRITE(numout,*) ' bdy debug int U', ir ; CALL FLUSH(numout)
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'V', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyint(ib_bdy,3,:,ir), lrecv = lrecv_bdyint(ib_bdy,3,:,ir) )    
!!$            IF(lwp) WRITE(numout,*) ' bdy debug int V', ir ; CALL FLUSH(numout)
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'T', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyext(ib_bdy,1,:,ir), lrecv = lrecv_bdyext(ib_bdy,1,:,ir) )
!!$            IF(lwp) WRITE(numout,*) ' bdy debug ext T', ir ; CALL FLUSH(numout)
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'U', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyext(ib_bdy,2,:,ir), lrecv = lrecv_bdyext(ib_bdy,2,:,ir) )
!!$            IF(lwp) WRITE(numout,*) ' bdy debug ext U', ir ; CALL FLUSH(numout)
!!$            zzbdy(:,:) = narea ; CALL lbc_lnk('bdy debug', zzbdy, 'V', 1._wp, kfillmode = jpfillnothing,   &
!!$               &                              lsend = lsend_bdyext(ib_bdy,3,:,ir), lrecv = lrecv_bdyext(ib_bdy,3,:,ir) )    
!!$            IF(lwp) WRITE(numout,*) ' bdy debug ext V', ir ; CALL FLUSH(numout)
!!$         END DO
         
      END DO   ! ib_bdy

      DO ib_bdy = 1,nb_bdy
         IF(  cn_dyn2d(ib_bdy) == 'orlanski' .OR. cn_dyn2d(ib_bdy) == 'orlanski_npo' .OR. &
            & cn_dyn3d(ib_bdy) == 'orlanski' .OR. cn_dyn3d(ib_bdy) == 'orlanski_npo' .OR. &
            & cn_tra(ib_bdy)   == 'orlanski' .OR. cn_tra(ib_bdy)   == 'orlanski_npo'      ) THEN
            DO igrd = 1, jpbgrd
               DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
                  ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
                  ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
                  IF(  mig0(ii) > 2 .AND. mig0(ii) < Ni0glo-2 .AND. mjg0(ij) > 2 .AND. mjg0(ij) < Nj0glo-2  ) THEN
                     WRITE(ctmp1,*) ' Orlanski is not safe when the open boundaries are on the interior of the computational domain'
                     CALL ctl_stop( ctmp1 )
                  END IF
               END DO
            END DO
         END IF
      END DO
      !
      DEALLOCATE( nbidta, nbjdta, nbrdta )
      !
   END SUBROUTINE bdy_def


   SUBROUTINE bdy_rim_treat( pumask, pvmask, pfmask, lrim0 )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_rim_treat  ***
      !!
      !! ** Purpose :   Initialize structures ( flagu, flagv, ntreat ) indicating how rim points
      !!                  are to be handled in the boundary condition treatment
      !!
      !! ** Method  :   - to handle rim 0 zmasks must indicate ocean points      (set at one on rim 0 and rim 1 and interior)
      !!                            and bdymasks must be set at 0 on rim 0       (set at one on rim 1 and interior)
      !!                    (as if rim 1 was free ocean)
      !!                - to handle rim 1 zmasks must be set at 0 on rim 0       (set at one on rim 1 and interior)
      !!                            and bdymasks must indicate free ocean points (set at one on interior)
      !!                    (as if rim 0 was land)
      !!                - we can then check in which direction the interior of the computational domain is with the difference
      !!                         mask array values on both sides to compute flagu and flagv
      !!                - and look at the ocean neighbours to compute ntreat
      !!----------------------------------------------------------------------
      REAL(wp), TARGET, DIMENSION(jpi,jpj), INTENT (in   ) :: pumask, pvmask   ! temporary u/v mask array
      REAL(wp), TARGET, DIMENSION(jpi,jpj), INTENT (in   ) :: pfmask           ! temporary fmask excluding coastal boundary condition (shlat)
      LOGICAL                             , INTENT (in   ) :: lrim0            ! .true. -> rim 0   .false. -> rim 1
      INTEGER  ::   ib_bdy, ii, ij, igrd, ib, icount       ! dummy loop indices
      INTEGER  ::   i_offset, j_offset, inn                ! local integer
      INTEGER  ::   ibeg, iend                             ! local integer
      LOGICAL  ::   llnon, llson, llean, llwen             ! local logicals indicating the presence of a ocean neighbour
      REAL(wp), POINTER, DIMENSION(:,:)       ::   zmask   ! pointer to 2D mask fields
      REAL(wp) ::   zefl, zwfl, znfl, zsfl                 ! local scalars
      CHARACTER(LEN=1), DIMENSION(jpbgrd)     ::   cgrid
      REAL(wp)        , DIMENSION(jpi,jpj)    ::   ztmp
      !!----------------------------------------------------------------------

      cgrid = (/'t','u','v'/)

      DO ib_bdy = 1, nb_bdy       ! Indices and directions of rim velocity components

         DO igrd = 1, jpbgrd
            
            IF( lrim0 ) THEN   ! extent of rim 0
               ibeg = 1                                     ;   iend = idx_bdy(ib_bdy)%nblenrim0(igrd)
            ELSE               ! extent of rim 1
               ibeg = idx_bdy(ib_bdy)%nblenrim0(igrd) + 1   ;   iend = idx_bdy(ib_bdy)%nblenrim(igrd)
            END IF

            ! Calculate relationship of U direction to the local orientation of the boundary
            ! flagu = -1 : u component is normal to the dynamical boundary and its direction is outward
            ! flagu =  0 : u is tangential
            ! flagu =  1 : u is normal to the boundary and is direction is inward
            SELECT CASE( igrd )
               CASE( 1 )   ;   zmask => pumask     ;   i_offset = 0   ! U(i-1)   T(i)   U(i  )
               CASE( 2 )   ;   zmask => bdytmask   ;   i_offset = 1   ! T(i  )   U(i)   T(i+1)
               CASE( 3 )   ;   zmask => pfmask     ;   i_offset = 0   ! F(i-1)   V(i)   F(i  )
            END SELECT 
            icount = 0
            ztmp(:,:) = -999._wp
            DO ib = ibeg, iend 
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               IF( ii < Nis0 .OR. ii > Nie0 .OR. ij < Njs0 .OR. ij > Nje0 )  CYCLE   ! call lbc_lnk -> no need to compute these pts
               zwfl = zmask(ii+i_offset-1,ij)
               zefl = zmask(ii+i_offset  ,ij)
               ! This error check only works if you are using the bdyXmask arrays (which are set to 0 on rims)
               IF( i_offset == 1 .and. zefl + zwfl == 2._wp ) THEN
                  icount = icount + 1
                  IF(lwp) WRITE(numout,*) 'Problem with igrd = ',igrd,' at (global) nbi, nbj : ',mig(ii),mjg(ij)
               ELSE
                  ztmp(ii,ij) = -zwfl + zefl
               ENDIF
            END DO
            IF( icount /= 0 ) THEN
               WRITE(ctmp1,*) 'Some ',cgrid(igrd),' grid points,',   &
                  ' are not boundary points (flagu calculation). Check nbi, nbj, indices for boundary set ',ib_bdy
               CALL ctl_stop( ctmp1 )
            ENDIF 
            SELECT CASE( igrd )
               CASE( 1 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'T', 1.0_wp )
               CASE( 2 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'U', 1.0_wp )
               CASE( 3 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'V', 1.0_wp )
            END SELECT 
            DO ib = ibeg, iend
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               idx_bdy(ib_bdy)%flagu(ib,igrd) = ztmp(ii,ij)
            END DO
            
            ! Calculate relationship of V direction to the local orientation of the boundary
            ! flagv = -1 : v component is normal to the dynamical boundary but its direction is outward
            ! flagv =  0 : v is tangential
            ! flagv =  1 : v is normal to the boundary and is direction is inward
            SELECT CASE( igrd )
               CASE( 1 )   ;   zmask => pvmask     ;   j_offset = 0
               CASE( 2 )   ;   zmask => pfmask     ;   j_offset = 0
               CASE( 3 )   ;   zmask => bdytmask   ;   j_offset = 1
            END SELECT 
            icount = 0
            ztmp(:,:) = -999._wp
            DO ib = ibeg, iend
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               IF( ii < Nis0 .OR. ii > Nie0 .OR. ij < Njs0 .OR. ij > Nje0 )  CYCLE   ! call lbc_lnk -> no need to compute these pts
               zsfl = zmask(ii,ij+j_offset-1)
               znfl = zmask(ii,ij+j_offset  )
               ! This error check only works if you are using the bdyXmask arrays (which are set to 0 on rims)
               IF( j_offset == 1 .and. znfl + zsfl == 2._wp ) THEN
                  IF(lwp) WRITE(numout,*) 'Problem with igrd = ',igrd,' at (global) nbi, nbj : ',mig(ii),mjg(ij)
                  icount = icount + 1
               ELSE
                  ztmp(ii,ij) = -zsfl + znfl
               END IF
            END DO
            IF( icount /= 0 ) THEN
               WRITE(ctmp1,*) 'Some ',cgrid(igrd),' grid points,',   &
                  ' are not boundary points (flagv calculation). Check nbi, nbj, indices for boundary set ',ib_bdy
               CALL ctl_stop( ctmp1 )
            ENDIF
            SELECT CASE( igrd )
               CASE( 1 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'T', 1.0_wp )
               CASE( 2 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'U', 1.0_wp )
               CASE( 3 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'V', 1.0_wp )
            END SELECT 
            DO ib = ibeg, iend
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               idx_bdy(ib_bdy)%flagv(ib,igrd) = ztmp(ii,ij)
            END DO
      
            ! Calculate ntreat
            SELECT CASE( igrd )
               CASE( 1 )   ;   zmask => bdytmask 
               CASE( 2 )   ;   zmask => bdyumask 
               CASE( 3 )   ;   zmask => bdyvmask 
            END SELECT
            ztmp(:,:) = -999._wp
            DO ib = ibeg, iend
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               IF( ii < Nis0 .OR. ii > Nie0 .OR. ij < Njs0 .OR. ij > Nje0 )  CYCLE   ! call lbc_lnk -> no need to compute these pts
               llnon = zmask(ii  ,ij+1) == 1._wp  
               llson = zmask(ii  ,ij-1) == 1._wp 
               llean = zmask(ii+1,ij  ) == 1._wp 
               llwen = zmask(ii-1,ij  ) == 1._wp 
               inn  = COUNT( (/ llnon, llson, llean, llwen /) )
               IF( inn == 0 ) THEN   ! no neighbours -> interior of a corner  or  cluster of rim points
                  !               !              !     _____     !     _____    !    __     __
                  !  1 |   o      !  2  o   |    !  3 | x        !  4     x |   !      |   |   -> error
                  !    |_x_ _     !    _ _x_|    !    |   o      !      o   |   !      |x_x|
                  IF(     zmask(ii+1,ij+1) == 1._wp ) THEN   ;   ztmp(ii,ij) = 1._wp
                  ELSEIF( zmask(ii-1,ij+1) == 1._wp ) THEN   ;   ztmp(ii,ij) = 2._wp
                  ELSEIF( zmask(ii+1,ij-1) == 1._wp ) THEN   ;   ztmp(ii,ij) = 3._wp
                  ELSEIF( zmask(ii-1,ij-1) == 1._wp ) THEN   ;   ztmp(ii,ij) = 4._wp
                  ELSE                                       ;   ztmp(ii,ij) = -1._wp
                     WRITE(ctmp1,*) 'Problem with  ',cgrid(igrd) ,' grid point', ii, ij,   &
                       ' on boundary set ', ib_bdy, ' has no free ocean neighbour'
                     IF( lrim0 ) THEN
                        WRITE(ctmp2,*) ' There seems to be a cluster of rim 0 points.'
                     ELSE
                        WRITE(ctmp2,*) ' There seems to be a cluster of rim 1 points.'
                     END IF
                     CALL ctl_warn( ctmp1, ctmp2 )
                  END IF
               END IF
               IF( inn == 1 ) THEN   ! middle of linear bdy  or incomplete corner  ! ___ o
                  !    |         !         |   !      o     !    ______            !    |x___
                  ! 5  | x o     ! 6   o x |   ! 7  __x__   ! 8    x
                  !    |         !         |   !            !      o
                  IF( llean )   ztmp(ii,ij) = 5._wp
                  IF( llwen )   ztmp(ii,ij) = 6._wp
                  IF( llnon )   ztmp(ii,ij) = 7._wp
                  IF( llson )   ztmp(ii,ij) = 8._wp
               END IF
               IF( inn == 2 ) THEN   ! exterior of a corner
                  !        o      !        o      !    _____|       !       |_____  
                  !  9 ____x o    ! 10   o x___   ! 11     x o      ! 12   o x      
                  !         |     !       |       !        o        !        o 
                  IF( llnon .AND. llean )   ztmp(ii,ij) =  9._wp
                  IF( llnon .AND. llwen )   ztmp(ii,ij) = 10._wp
                  IF( llson .AND. llean )   ztmp(ii,ij) = 11._wp
                  IF( llson .AND. llwen )   ztmp(ii,ij) = 12._wp
               END IF
               IF( inn == 3 ) THEN   ! 3 neighbours     __   __
                  !    |_  o      !        o  _|  !       |_|     !       o         
                  ! 13  _| x o    ! 14   o x |_   ! 15   o x o    ! 16  o x o       
                  !    |   o      !        o   |  !        o      !    __|¨|__    
                  IF( llnon .AND. llean .AND. llson )   ztmp(ii,ij) = 13._wp
                  IF( llnon .AND. llwen .AND. llson )   ztmp(ii,ij) = 14._wp
                  IF( llwen .AND. llson .AND. llean )   ztmp(ii,ij) = 15._wp
                  IF( llwen .AND. llnon .AND. llean )   ztmp(ii,ij) = 16._wp
               END IF
               IF( inn == 4 ) THEN
                  WRITE(ctmp1,*)  'Problem with  ',cgrid(igrd) ,' grid point', ii, ij,   &
                       ' on boundary set ', ib_bdy, ' have 4 neighbours'
                  CALL ctl_stop( ctmp1 )
               END IF
            END DO
            SELECT CASE( igrd )
               CASE( 1 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'T', 1.0_wp )
               CASE( 2 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'U', 1.0_wp )
               CASE( 3 )   ;   CALL lbc_lnk( 'bdyini', ztmp, 'V', 1.0_wp )
            END SELECT 
            DO ib = ibeg, iend
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               idx_bdy(ib_bdy)%ntreat(ib,igrd) = NINT(ztmp(ii,ij))
            END DO
            !
         END DO   ! jpbgrd
         !
      END DO   ! ib_bdy

    END SUBROUTINE bdy_rim_treat

   
    SUBROUTINE find_neib( ii, ij, itreat, ii1, ij1, ii2, ij2, ii3, ij3 )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE find_neib  ***
      !!
      !! ** Purpose :   get ii1, ij1, ii2, ij2, ii3, ij3, the indices of
      !!               the free ocean neighbours of (ii,ij) for bdy treatment
      !!
      !! ** Method  :  use itreat input to select a case
      !!               N.B. ntreat is defined for all bdy points in routine bdy_rim_treat
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   )      ::   ii, ij, itreat
      INTEGER, INTENT(  out)      ::   ii1, ij1, ii2, ij2, ii3, ij3
      !!----------------------------------------------------------------------
      SELECT CASE( itreat )   ! points that will be used by bdy routines, -1 will be discarded
         !               !               !     _____     !     _____     
         !  1 |   o      !  2  o   |     !  3 | x        !  4     x |    
         !    |_x_ _     !    _ _x_|     !    |   o      !      o   |
      CASE( 1 )    ;   ii1 = ii+1   ;   ij1 = ij+1   ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
      CASE( 2 )    ;   ii1 = ii-1   ;   ij1 = ij+1   ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
      CASE( 3 )    ;   ii1 = ii+1   ;   ij1 = ij-1   ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
      CASE( 4 )    ;   ii1 = ii-1   ;   ij1 = ij-1   ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
         !    |          !         |     !      o        !    ______                   ! or incomplete corner
         ! 5  | x o      ! 6   o x |     ! 7  __x__      ! 8    x                      !  7  ____ o
         !    |          !         |     !               !      o                      !         |x___
      CASE( 5 )    ;   ii1 = ii+1   ;   ij1 = ij     ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
      CASE( 6 )    ;   ii1 = ii-1   ;   ij1 = ij     ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
      CASE( 7 )    ;   ii1 = ii     ;   ij1 = ij+1   ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
      CASE( 8 )    ;   ii1 = ii     ;   ij1 = ij-1   ;   ii2 = -1     ;   ij2 = -1     ;   ii3 = -1     ;   ij3 = -1
         !        o      !        o      !    _____|     !       |_____  
         !  9 ____x o    ! 10   o x___   ! 11     x o    ! 12   o x      
         !         |     !       |       !        o      !        o      
      CASE( 9  )   ;   ii1 = ii     ;   ij1 = ij+1   ;   ii2 = ii+1   ;   ij2 = ij     ;   ii3 = -1     ;   ij3 = -1 
      CASE( 10 )   ;   ii1 = ii     ;   ij1 = ij+1   ;   ii2 = ii-1   ;   ij2 = ij     ;   ii3 = -1     ;   ij3 = -1
      CASE( 11 )   ;   ii1 = ii     ;   ij1 = ij-1   ;   ii2 = ii+1   ;   ij2 = ij     ;   ii3 = -1     ;   ij3 = -1
      CASE( 12 )   ;   ii1 = ii     ;   ij1 = ij-1   ;   ii2 = ii-1   ;   ij2 = ij     ;   ii3 = -1     ;   ij3 = -1
         !    |_  o      !        o  _|  !     ¨¨|_|¨¨   !       o         
         ! 13  _| x o    !  14  o x |_   !  15  o x o    ! 16  o x o       
         !    |   o      !        o   |  !        o      !    __|¨|__ 
      CASE( 13 )   ;   ii1 = ii     ;   ij1 = ij+1   ;   ii2 = ii+1   ;   ij2 = ij     ;   ii3 = ii     ;   ij3 = ij-1   
      CASE( 14 )   ;   ii1 = ii     ;   ij1 = ij+1   ;   ii2 = ii-1   ;   ij2 = ij     ;   ii3 = ii     ;   ij3 = ij-1 
      CASE( 15 )   ;   ii1 = ii-1   ;   ij1 = ij     ;   ii2 = ii     ;   ij2 = ij-1   ;   ii3 = ii+1   ;   ij3 = ij   
      CASE( 16 )   ;   ii1 = ii-1   ;   ij1 = ij     ;   ii2 = ii     ;   ij2 = ij+1   ;   ii3 = ii+1   ;   ij3 = ij
      END SELECT
   END SUBROUTINE find_neib
    

   SUBROUTINE bdy_read_seg( kb_bdy, knblendta ) 
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_read_seg  ***
      !!
      !! ** Purpose :  build bdy coordinates with segments defined in namelist
      !!
      !! ** Method  :  read namelist nambdy_index blocks
      !!
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT (in   ) ::   kb_bdy           ! bdy number
      INTEGER, DIMENSION(jpbgrd), INTENT (  out) ::   knblendta        ! length of index arrays 
      !!
      INTEGER          ::   ios                 ! Local integer output status for namelist read
      INTEGER          ::   nbdyind, nbdybeg, nbdyend
      INTEGER          ::   nbdy_count, nbdy_rdstart, nbdy_loc
      CHARACTER(LEN=1) ::   ctypebdy   !     -        - 
      CHARACTER(LEN=50)::   cerrmsg    !     -        - 
      NAMELIST/nambdy_index/ ctypebdy, nbdyind, nbdybeg, nbdyend
      !!----------------------------------------------------------------------
      ! Need to support possibility of reading more than one nambdy_index from
      ! the namelist_cfg internal file.
      ! Do this by finding the kb_bdy'th occurence of nambdy_index in the
      ! character buffer as the starting point.
      nbdy_rdstart = 1
      DO nbdy_count = 1, kb_bdy
       nbdy_loc = INDEX( numnam_cfg( nbdy_rdstart: ), 'nambdy_index' )
       IF( nbdy_loc .GT. 0 ) THEN
          nbdy_rdstart = nbdy_rdstart + nbdy_loc
       ELSE
          WRITE(cerrmsg,'(A,I4,A)') 'Error: entry number ',kb_bdy,' of nambdy_index not found'
          ios = -1
          CALL ctl_nam ( ios , cerrmsg )
       ENDIF
      END DO
      nbdy_rdstart = MAX( 1, nbdy_rdstart - 2 )
      READ  ( numnam_cfg( nbdy_rdstart: ), nambdy_index, IOSTAT = ios, ERR = 904)
904   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy_index in configuration namelist' )
      IF(lwm) WRITE ( numond, nambdy_index )
      
      SELECT CASE ( TRIM(ctypebdy) )
      CASE( 'N' )
         IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
            nbdyind  = Nj0glo - 2  ! set boundary to whole side of model domain.
            nbdybeg  = 2
            nbdyend  = Ni0glo - 1
         ENDIF
         nbdysegn = nbdysegn + 1
         npckgn(nbdysegn) = kb_bdy ! Save bdy package number
         jpjnob(nbdysegn) = nbdyind 
         jpindt(nbdysegn) = nbdybeg
         jpinft(nbdysegn) = nbdyend
         !
      CASE( 'S' )
         IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
            nbdyind  = 2           ! set boundary to whole side of model domain.
            nbdybeg  = 2
            nbdyend  = Ni0glo - 1
         ENDIF
         nbdysegs = nbdysegs + 1
         npckgs(nbdysegs) = kb_bdy ! Save bdy package number
         jpjsob(nbdysegs) = nbdyind
         jpisdt(nbdysegs) = nbdybeg
         jpisft(nbdysegs) = nbdyend
         !
      CASE( 'E' )
         IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
            nbdyind  = Ni0glo - 2  ! set boundary to whole side of model domain.
            nbdybeg  = 2
            nbdyend  = Nj0glo - 1
         ENDIF
         nbdysege = nbdysege + 1 
         npckge(nbdysege) = kb_bdy ! Save bdy package number
         jpieob(nbdysege) = nbdyind
         jpjedt(nbdysege) = nbdybeg
         jpjeft(nbdysege) = nbdyend
         !
      CASE( 'W' )
         IF( nbdyind == -1 ) THEN  ! Automatic boundary definition: if nbdysegX = -1
            nbdyind  = 2           ! set boundary to whole side of model domain.
            nbdybeg  = 2
            nbdyend  = Nj0glo - 1
         ENDIF
         nbdysegw = nbdysegw + 1
         npckgw(nbdysegw) = kb_bdy ! Save bdy package number
         jpiwob(nbdysegw) = nbdyind
         jpjwdt(nbdysegw) = nbdybeg
         jpjwft(nbdysegw) = nbdyend
         !
      CASE DEFAULT   ;   CALL ctl_stop( 'ctypebdy must be N, S, E or W' )
      END SELECT
      
      ! For simplicity we assume that in case of straight bdy, arrays have the same length
      ! (even if it is true that last tangential velocity points
      ! are useless). This simplifies a little bit boundary data format (and agrees with format
      ! used so far in obc package)
      
      knblendta(1:jpbgrd) =  (nbdyend - nbdybeg + 1) * nn_rimwidth(kb_bdy)
      
   END SUBROUTINE bdy_read_seg

   
   SUBROUTINE bdy_ctl_seg
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_ctl_seg  ***
      !!
      !! ** Purpose :   Check straight open boundary segments location
      !!
      !! ** Method  :   - Look for open boundary corners
      !!                - Check that segments start or end on land 
      !!----------------------------------------------------------------------
      INTEGER  ::   ib, ib1, ib2, ji ,jj, itest  
      INTEGER, DIMENSION(jp_nseg,2) :: icorne, icornw, icornn, icorns  
      REAL(wp), DIMENSION(2) ::   ztestmask
      !!----------------------------------------------------------------------
      !
      IF (lwp) WRITE(numout,*) ' '
      IF (lwp) WRITE(numout,*) 'bdy_ctl_seg: Check analytical segments'
      IF (lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      !
      IF(lwp) WRITE(numout,*) 'Number of east  segments     : ', nbdysege
      IF(lwp) WRITE(numout,*) 'Number of west  segments     : ', nbdysegw
      IF(lwp) WRITE(numout,*) 'Number of north segments     : ', nbdysegn
      IF(lwp) WRITE(numout,*) 'Number of south segments     : ', nbdysegs
      !
      ! 1. Check bounds
      !----------------
      DO ib = 1, nbdysegn
         IF (lwp) WRITE(numout,*) '**check north seg bounds pckg: ', npckgn(ib)
         IF ((jpjnob(ib).ge.Nj0glo-1).or.& 
            &(jpjnob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpindt(ib).ge.jpinft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpindt(ib).lt.1     )     CALL ctl_stop( 'Start index out of domain' )
         IF (jpinft(ib).gt.Ni0glo)     CALL ctl_stop( 'End index out of domain' )
      END DO
      !
      DO ib = 1, nbdysegs
         IF (lwp) WRITE(numout,*) '**check south seg bounds pckg: ', npckgs(ib)
         IF ((jpjsob(ib).ge.Nj0glo-1).or.& 
            &(jpjsob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpisdt(ib).ge.jpisft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpisdt(ib).lt.1     )     CALL ctl_stop( 'Start index out of domain' )
         IF (jpisft(ib).gt.Ni0glo)     CALL ctl_stop( 'End index out of domain' )
      END DO
      !
      DO ib = 1, nbdysege
         IF (lwp) WRITE(numout,*) '**check east  seg bounds pckg: ', npckge(ib)
         IF ((jpieob(ib).ge.Ni0glo-1).or.& 
            &(jpieob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpjedt(ib).ge.jpjeft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpjedt(ib).lt.1     )     CALL ctl_stop( 'Start index out of domain' )
         IF (jpjeft(ib).gt.Nj0glo)     CALL ctl_stop( 'End index out of domain' )
      END DO
      !
      DO ib = 1, nbdysegw
         IF (lwp) WRITE(numout,*) '**check west  seg bounds pckg: ', npckgw(ib)
         IF ((jpiwob(ib).ge.Ni0glo-1).or.& 
            &(jpiwob(ib).le.1))        CALL ctl_stop( 'nbdyind out of domain' )
         IF (jpjwdt(ib).ge.jpjwft(ib)) CALL ctl_stop( 'Bdy start index is greater than end index' )
         IF (jpjwdt(ib).lt.1     )     CALL ctl_stop( 'Start index out of domain' )
         IF (jpjwft(ib).gt.Nj0glo)     CALL ctl_stop( 'End index out of domain' )
      ENDDO
      !      
      ! 2. Look for segment crossings
      !------------------------------ 
      IF (lwp) WRITE(numout,*) '**Look for segments corners  :'
      !
      itest = 0 ! corner number
      !
      ! flag to detect if start or end of open boundary belongs to a corner
      ! if not (=0), it must be on land.
      ! if a corner is detected, save bdy package number for further tests
      icorne(:,:)=0. ; icornw(:,:)=0. ; icornn(:,:)=0. ; icorns(:,:)=0.
      ! South/West crossings
      IF ((nbdysegw > 0).AND.(nbdysegs > 0)) THEN
         DO ib1 = 1, nbdysegw        
            DO ib2 = 1, nbdysegs
               IF (( jpisdt(ib2)<=jpiwob(ib1)).AND. &
                &  ( jpisft(ib2)>=jpiwob(ib1)).AND. &
                &  ( jpjwdt(ib1)<=jpjsob(ib2)).AND. &
                &  ( jpjwft(ib1)>=jpjsob(ib2))) THEN
                  IF ((jpjwdt(ib1)==jpjsob(ib2)).AND.(jpisdt(ib2)==jpiwob(ib1))) THEN 
                     ! We have a possible South-West corner                      
!                     WRITE(numout,*) ' Found a South-West corner at (i,j): ', jpisdt(ib2), jpjwdt(ib1) 
!                     WRITE(numout,*) ' between segments: ', npckgw(ib1), npckgs(ib2)
                     icornw(ib1,1) = npckgs(ib2)
                     icorns(ib2,1) = npckgw(ib1)
                  ELSEIF ((jpisft(ib2)==jpiwob(ib1)).AND.(jpjwft(ib1)==jpjsob(ib2))) THEN
                     WRITE(ctmp1,*) ' Found an acute open boundary corner at point (i,j)= ', &
                        &                                     jpisft(ib2), jpjwft(ib1)
                     WRITE(ctmp2,*) ' Not allowed yet'
                     WRITE(ctmp3,*) ' Crossing problem with West segment: ',npckgw(ib1), & 
                        &                            ' and South segment: ',npckgs(ib2)
                     CALL ctl_stop( ctmp1, ctmp2, ctmp3 )
                  ELSE
                     WRITE(ctmp1,*) ' Check South and West Open boundary indices'
                     WRITE(ctmp2,*) ' Crossing problem with West segment: ',npckgw(ib1) , &
                        &                            ' and South segment: ',npckgs(ib2)
                     CALL ctl_stop( ctmp1, ctmp2 )
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! South/East crossings
      IF ((nbdysege > 0).AND.(nbdysegs > 0)) THEN
         DO ib1 = 1, nbdysege
            DO ib2 = 1, nbdysegs
               IF (( jpisdt(ib2)<=jpieob(ib1)+1).AND. &
                &  ( jpisft(ib2)>=jpieob(ib1)+1).AND. &
                &  ( jpjedt(ib1)<=jpjsob(ib2)  ).AND. &
                &  ( jpjeft(ib1)>=jpjsob(ib2)  )) THEN
                  IF ((jpjedt(ib1)==jpjsob(ib2)).AND.(jpisft(ib2)==jpieob(ib1)+1)) THEN
                     ! We have a possible South-East corner 
!                     WRITE(numout,*) ' Found a South-East corner at (i,j): ', jpisft(ib2), jpjedt(ib1) 
!                     WRITE(numout,*) ' between segments: ', npckge(ib1), npckgs(ib2)
                     icorne(ib1,1) = npckgs(ib2)
                     icorns(ib2,2) = npckge(ib1)
                  ELSEIF ((jpjeft(ib1)==jpjsob(ib2)).AND.(jpisdt(ib2)==jpieob(ib1)+1)) THEN
                     WRITE(ctmp1,*) ' Found an acute open boundary corner at point (i,j)= ', &
                        &                                     jpisdt(ib2), jpjeft(ib1)
                     WRITE(ctmp2,*) ' Not allowed yet'
                     WRITE(ctmp3,*) ' Crossing problem with East segment: ',npckge(ib1), &
                        &                            ' and South segment: ',npckgs(ib2)
                     CALL ctl_stop( ctmp1, ctmp2, ctmp3 )
                  ELSE
                     WRITE(ctmp1,*) ' Check South and East Open boundary indices'
                     WRITE(ctmp2,*) ' Crossing problem with East segment: ',npckge(ib1), &
                     &                               ' and South segment: ',npckgs(ib2)
                     CALL ctl_stop( ctmp1, ctmp2 )
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! North/West crossings
      IF ((nbdysegn > 0).AND.(nbdysegw > 0)) THEN
         DO ib1 = 1, nbdysegw        
            DO ib2 = 1, nbdysegn
               IF (( jpindt(ib2)<=jpiwob(ib1)  ).AND. &
                &  ( jpinft(ib2)>=jpiwob(ib1)  ).AND. &
                &  ( jpjwdt(ib1)<=jpjnob(ib2)+1).AND. &
                &  ( jpjwft(ib1)>=jpjnob(ib2)+1)) THEN
                  IF ((jpjwft(ib1)==jpjnob(ib2)+1).AND.(jpindt(ib2)==jpiwob(ib1))) THEN
                     ! We have a possible North-West corner 
!                     WRITE(numout,*) ' Found a North-West corner at (i,j): ', jpindt(ib2), jpjwft(ib1) 
!                     WRITE(numout,*) ' between segments: ', npckgw(ib1), npckgn(ib2)
                     icornw(ib1,2) = npckgn(ib2)
                     icornn(ib2,1) = npckgw(ib1)
                  ELSEIF ((jpjwdt(ib1)==jpjnob(ib2)+1).AND.(jpinft(ib2)==jpiwob(ib1))) THEN
                     WRITE(ctmp1,*) ' Found an acute open boundary corner at point (i,j)= ', &
                     &                                     jpinft(ib2), jpjwdt(ib1)
                     WRITE(ctmp2,*) ' Not allowed yet'
                     WRITE(ctmp3,*) ' Crossing problem with West segment: ',npckgw(ib1), &
                     &                               ' and North segment: ',npckgn(ib2)
                     CALL ctl_stop( ctmp1, ctmp2, ctmp3 )
                  ELSE
                     WRITE(ctmp1,*) ' Check North and West Open boundary indices'
                     WRITE(ctmp2,*) ' Crossing problem with West segment: ',npckgw(ib1), &
                     &                               ' and North segment: ',npckgn(ib2)
                     CALL ctl_stop( ctmp1, ctmp2 )
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! North/East crossings
      IF ((nbdysegn > 0).AND.(nbdysege > 0)) THEN
         DO ib1 = 1, nbdysege        
            DO ib2 = 1, nbdysegn
               IF (( jpindt(ib2)<=jpieob(ib1)+1).AND. &
                &  ( jpinft(ib2)>=jpieob(ib1)+1).AND. &
                &  ( jpjedt(ib1)<=jpjnob(ib2)+1).AND. &
                &  ( jpjeft(ib1)>=jpjnob(ib2)+1)) THEN
                  IF ((jpjeft(ib1)==jpjnob(ib2)+1).AND.(jpinft(ib2)==jpieob(ib1)+1)) THEN
                     ! We have a possible North-East corner 
!                     WRITE(numout,*) ' Found a North-East corner at (i,j): ', jpinft(ib2), jpjeft(ib1)
!                     WRITE(numout,*) ' between segments: ', npckge(ib1), npckgn(ib2)
                     icorne(ib1,2) = npckgn(ib2)
                     icornn(ib2,2) = npckge(ib1)
                  ELSEIF ((jpjedt(ib1)==jpjnob(ib2)+1).AND.(jpindt(ib2)==jpieob(ib1)+1)) THEN
                     WRITE(ctmp1,*) ' Found an acute open boundary corner at point (i,j)= ', &
                     &                                     jpindt(ib2), jpjedt(ib1)
                     WRITE(ctmp2,*) ' Not allowed yet'
                     WRITE(ctmp3,*) ' Crossing problem with East segment: ',npckge(ib1), &
                     &                               ' and North segment: ',npckgn(ib2)
                     CALL ctl_stop( ctmp1, ctmp2, ctmp3 )
                  ELSE
                     WRITE(ctmp1,*) ' Check North and East Open boundary indices'
                     WRITE(ctmp2,*) ' Crossing problem with East segment: ',npckge(ib1), &
                     &                               ' and North segment: ',npckgn(ib2)
                     CALL ctl_stop( ctmp1, ctmp2 )
                  END IF
               END IF
            END DO
         END DO
      END IF
      !
      ! 3. Check if segment extremities are on land
      !-------------------------------------------- 
      !
      ! West segments
      DO ib = 1, nbdysegw
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
              IF( mig0(ji) == jpiwob(ib) .AND. mjg0(jj) == jpjwdt(ib) )   ztestmask(1) = tmask(ji,jj,1)
              IF( mig0(ji) == jpiwob(ib) .AND. mjg0(jj) == jpjwft(ib) )   ztestmask(2) = tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF (ztestmask(1)==1) THEN 
            IF (icornw(ib,1)==0) THEN
               WRITE(ctmp1,*) ' Open boundary segment ', npckgw(ib)
               CALL ctl_stop( ctmp1, ' does not start on land or on a corner' )
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a South-West corner at (i,j): ', jpiwob(ib), jpjwdt(ib)
               CALL bdy_ctl_corn(npckgw(ib), icornw(ib,1))
               itest=itest+1
            ENDIF
         ENDIF
         IF (ztestmask(2)==1) THEN
            IF (icornw(ib,2)==0) THEN
               WRITE(ctmp1,*) ' Open boundary segment ', npckgw(ib)
               CALL ctl_stop( ' ', ctmp1, ' does not end on land or on a corner' )
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a North-West corner at (i,j): ', jpiwob(ib), jpjwft(ib)
               CALL bdy_ctl_corn(npckgw(ib), icornw(ib,2))
               itest=itest+1
            ENDIF
         ENDIF
      END DO
      !
      ! East segments
      DO ib = 1, nbdysege
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
              IF( mig0(ji) == jpieob(ib)+1 .AND. mjg0(jj) == jpjedt(ib) )   ztestmask(1) = tmask(ji,jj,1)
              IF( mig0(ji) == jpieob(ib)+1 .AND. mjg0(jj) == jpjeft(ib) )   ztestmask(2) = tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF (ztestmask(1)==1) THEN
            IF (icorne(ib,1)==0) THEN
               WRITE(ctmp1,*) ' Open boundary segment ', npckge(ib)
               CALL ctl_stop( ctmp1, ' does not start on land or on a corner' )
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a South-East corner at (i,j): ', jpieob(ib)+1, jpjedt(ib)
               CALL bdy_ctl_corn(npckge(ib), icorne(ib,1))
               itest=itest+1
            ENDIF
         ENDIF
         IF (ztestmask(2)==1) THEN
            IF (icorne(ib,2)==0) THEN
               WRITE(ctmp1,*) ' Open boundary segment ', npckge(ib)
               CALL ctl_stop( ctmp1, ' does not end on land or on a corner' )
            ELSE
               ! This is a corner
               IF(lwp) WRITE(numout,*) 'Found a North-East corner at (i,j): ', jpieob(ib)+1, jpjeft(ib)
               CALL bdy_ctl_corn(npckge(ib), icorne(ib,2))
               itest=itest+1
            ENDIF
         ENDIF
      END DO
      !
      ! South segments
      DO ib = 1, nbdysegs
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
              IF( mjg0(jj) == jpjsob(ib) .AND. mig0(ji) == jpisdt(ib) )   ztestmask(1) = tmask(ji,jj,1)
              IF( mjg0(jj) == jpjsob(ib) .AND. mig0(ji) == jpisft(ib) )   ztestmask(2) = tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF ((ztestmask(1)==1).AND.(icorns(ib,1)==0)) THEN
            WRITE(ctmp1,*) ' Open boundary segment ', npckgs(ib)
            CALL ctl_stop( ctmp1, ' does not start on land or on a corner' )
         ENDIF
         IF ((ztestmask(2)==1).AND.(icorns(ib,2)==0)) THEN
            WRITE(ctmp1,*) ' Open boundary segment ', npckgs(ib)
            CALL ctl_stop( ctmp1, ' does not end on land or on a corner' )
         ENDIF
      END DO
      !
      ! North segments
      DO ib = 1, nbdysegn
         ! get mask at boundary extremities:
         ztestmask(1:2)=0.
         DO ji = 1, jpi
            DO jj = 1, jpj             
               IF( mjg0(jj) == jpjnob(ib)+1 .AND. mig0(ji) == jpindt(ib) )   ztestmask(1) = tmask(ji,jj,1)
               IF( mjg0(jj) == jpjnob(ib)+1 .AND. mig0(ji) == jpinft(ib) )   ztestmask(2) = tmask(ji,jj,1)  
            END DO
         END DO
         CALL mpp_sum( 'bdyini', ztestmask, 2 )   ! sum over the global domain

         IF ((ztestmask(1)==1).AND.(icornn(ib,1)==0)) THEN
            WRITE(ctmp1,*) ' Open boundary segment ', npckgn(ib)
            CALL ctl_stop( ctmp1, ' does not start on land' )
         ENDIF
         IF ((ztestmask(2)==1).AND.(icornn(ib,2)==0)) THEN
            WRITE(ctmp1,*) ' Open boundary segment ', npckgn(ib)
            CALL ctl_stop( ctmp1, ' does not end on land' )
         ENDIF
      END DO
      !
      IF ((itest==0).AND.(lwp)) WRITE(numout,*) 'NO open boundary corner found'
      !
      ! Other tests TBD: 
      ! segments completly on land
      ! optimized open boundary array length according to landmask
      ! Nudging layers that overlap with interior domain
      !
   END SUBROUTINE bdy_ctl_seg

   
   SUBROUTINE bdy_coords_seg( nbidta, nbjdta, nbrdta ) 
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_coords_seg  ***
      !!
      !! ** Purpose :  build nbidta, nbidta, nbrdta for bdy built with segments
      !!
      !! ** Method  :  
      !!
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:,:), intent(  out)  ::   nbidta, nbjdta, nbrdta   ! Index arrays: i and j indices of bdy dta
      !!
      INTEGER  ::   ii, ij, ir, iseg
      INTEGER  ::   igrd         ! grid type (t=1, u=2, v=3)
      INTEGER  ::   icount       ! 
      INTEGER  ::   ib_bdy       ! bdy number
      !!----------------------------------------------------------------------

      ! East
      !-----
      DO iseg = 1, nbdysege
         ib_bdy = npckge(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjedt(iseg), jpjeft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 2 - ir + nn_hls
               nbjdta(icount, igrd, ib_bdy) = ij + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjedt(iseg), jpjeft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 1 - ir + nn_hls
               nbjdta(icount, igrd, ib_bdy) = ij + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            !            DO ij = jpjedt(iseg), jpjeft(iseg) - 1
            DO ij = jpjedt(iseg), jpjeft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpieob(iseg) + 2 - ir + nn_hls
               nbjdta(icount, igrd, ib_bdy) = ij + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
      ENDDO
      !
      ! West
      !-----
      DO iseg = 1, nbdysegw
         ib_bdy = npckgw(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjwdt(iseg), jpjwft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1 + nn_hls
               nbjdta(icount, igrd, ib_bdy) = ij + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ij = jpjwdt(iseg), jpjwft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1 + nn_hls
               nbjdta(icount, igrd, ib_bdy) = ij + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            !            DO ij = jpjwdt(iseg), jpjwft(iseg) - 1
            DO ij = jpjwdt(iseg), jpjwft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = jpiwob(iseg) + ir - 1 + nn_hls
               nbjdta(icount, igrd, ib_bdy) = ij + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
      ENDDO
      !
      ! North
      !-----
      DO iseg = 1, nbdysegn
         ib_bdy = npckgn(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpindt(iseg), jpinft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii + nn_hls
               nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 2 - ir + nn_hls 
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            !            DO ii = jpindt(iseg), jpinft(iseg) - 1
            DO ii = jpindt(iseg), jpinft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii + nn_hls
               nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 2 - ir + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpindt(iseg), jpinft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii + nn_hls
               nbjdta(icount, igrd, ib_bdy) = jpjnob(iseg) + 1 - ir + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
      ENDDO
      !
      ! South
      !-----
      DO iseg = 1, nbdysegs
         ib_bdy = npckgs(iseg)
         !
         ! ------------ T points -------------
         igrd=1
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpisdt(iseg), jpisft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii + nn_hls
               nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1 + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
         !
         ! ------------ U points -------------
         igrd=2
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            !            DO ii = jpisdt(iseg), jpisft(iseg) - 1
            DO ii = jpisdt(iseg), jpisft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii + nn_hls
               nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1 + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
            nbidta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
            nbjdta(icount, igrd, ib_bdy) = -ib_bdy ! Discount this point
         ENDDO
         !
         ! ------------ V points -------------
         igrd=3
         icount=0
         DO ir = 1, nn_rimwidth(ib_bdy)
            DO ii = jpisdt(iseg), jpisft(iseg)
               icount = icount + 1
               nbidta(icount, igrd, ib_bdy) = ii + nn_hls
               nbjdta(icount, igrd, ib_bdy) = jpjsob(iseg) + ir - 1 + nn_hls
               nbrdta(icount, igrd, ib_bdy) = ir
            ENDDO
         ENDDO
      ENDDO

      
   END SUBROUTINE bdy_coords_seg
   
   
   SUBROUTINE bdy_ctl_corn( ib1, ib2 )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_ctl_corn  ***
      !!
      !! ** Purpose :   Check numerical schemes consistency between
      !!                segments having a common corner
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  ::   ib1, ib2
      INTEGER :: itest
      !!----------------------------------------------------------------------
      itest = 0

      IF( cn_dyn2d(ib1) /= cn_dyn2d(ib2) )   itest = itest + 1
      IF( cn_dyn3d(ib1) /= cn_dyn3d(ib2) )   itest = itest + 1
      IF( cn_tra  (ib1) /= cn_tra  (ib2) )   itest = itest + 1
      !
      IF( nn_dyn2d_dta(ib1) /= nn_dyn2d_dta(ib2) )   itest = itest + 1
      IF( nn_dyn3d_dta(ib1) /= nn_dyn3d_dta(ib2) )   itest = itest + 1
      IF( nn_tra_dta  (ib1) /= nn_tra_dta  (ib2) )   itest = itest + 1
      !
      IF( nn_rimwidth(ib1) /= nn_rimwidth(ib2) )   itest = itest + 1   
      !
      IF( itest>0 ) THEN
         WRITE(ctmp1,*) ' Segments ', ib1, 'and ', ib2
         CALL ctl_stop( ctmp1, ' have different open bdy schemes' )
      ENDIF
      !
   END SUBROUTINE bdy_ctl_corn


   SUBROUTINE bdy_meshwri()
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_meshwri  ***
      !!         
      !! ** Purpose :   write netcdf file with nbr, flagu, flagv, ntreat for T, U 
      !!                and V points in 2D arrays for easier visualisation/control
      !!
      !! ** Method  :   use iom_rstput as in domwri.F
      !!----------------------------------------------------------------------      
      INTEGER  ::   ib_bdy, ii, ij, igrd, ib     ! dummy loop indices
      INTEGER  ::   inum                                   !   -       -
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zmask                   ! pointer to 2D mask fields
      REAL(wp)         , DIMENSION(jpi,jpj) ::   ztmp
      CHARACTER(LEN=1) , DIMENSION(jpbgrd)  ::   cgrid
      !!----------------------------------------------------------------------      
      cgrid = (/'t','u','v'/)
      CALL iom_open( 'bdy_mesh', inum, ldwrt = .TRUE. )
      DO igrd = 1, jpbgrd
         SELECT CASE( igrd )
         CASE( 1 )   ;   zmask => tmask(:,:,1)
         CASE( 2 )   ;   zmask => umask(:,:,1)
         CASE( 3 )   ;   zmask => vmask(:,:,1)
         END SELECT
         ztmp(:,:) = zmask(:,:)
         DO ib_bdy = 1, nb_bdy
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)      ! nbr deined for all rims
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               ztmp(ii,ij) = REAL(idx_bdy(ib_bdy)%nbr(ib,igrd), wp) + 10.
               IF( zmask(ii,ij) == 0. ) ztmp(ii,ij) = - ztmp(ii,ij)
            END DO
         END DO
         CALL iom_rstput( 0, 0, inum, 'bdy_nbr_'//cgrid(igrd), ztmp(:,:), ktype = jp_i4 )
         ztmp(:,:) = zmask(:,:)
         DO ib_bdy = 1, nb_bdy
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)   ! flagu defined only for rims 0 and 1
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               ztmp(ii,ij) = REAL(idx_bdy(ib_bdy)%flagu(ib,igrd), wp) + 10.
               IF( zmask(ii,ij) == 0. ) ztmp(ii,ij) = - ztmp(ii,ij)
            END DO
         END DO
         CALL iom_rstput( 0, 0, inum, 'flagu_'//cgrid(igrd), ztmp(:,:), ktype = jp_i4 )
         ztmp(:,:) = zmask(:,:)
         DO ib_bdy = 1, nb_bdy
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)   ! flagv defined only for rims 0 and 1
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               ztmp(ii,ij) = REAL(idx_bdy(ib_bdy)%flagv(ib,igrd), wp) + 10.
               IF( zmask(ii,ij) == 0. ) ztmp(ii,ij) = - ztmp(ii,ij)
            END DO
         END DO
         CALL iom_rstput( 0, 0, inum, 'flagv_'//cgrid(igrd), ztmp(:,:), ktype = jp_i4 )
         ztmp(:,:) = zmask(:,:)
         DO ib_bdy = 1, nb_bdy
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)   ! ntreat defined only for rims 0 and 1
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               ztmp(ii,ij) = REAL(idx_bdy(ib_bdy)%ntreat(ib,igrd), wp) + 10.
               IF( zmask(ii,ij) == 0. ) ztmp(ii,ij) = - ztmp(ii,ij)
            END DO
         END DO
         CALL iom_rstput( 0, 0, inum, 'ntreat_'//cgrid(igrd), ztmp(:,:), ktype = jp_i4 )
      END DO
      CALL iom_close( inum )

   END SUBROUTINE bdy_meshwri
   
   !!=================================================================================
END MODULE bdyini
