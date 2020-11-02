MODULE sbcwave
   !!======================================================================
   !!                       ***  MODULE  sbcwave  ***
   !! Wave module 
   !!======================================================================
   !! History :  3.3  !  2011-09  (M. Adani)  Original code: Drag Coefficient 
   !!         :  3.4  !  2012-10  (M. Adani)  Stokes Drift 
   !!            3.6  !  2014-09  (E. Clementi,P. Oddo) New Stokes Drift Computation
   !!             -   !  2016-12  (G. Madec, E. Clementi) update Stoke drift computation
   !!                                                    + add sbc_wave_ini routine
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_stokes    : calculate 3D Stokes-drift velocities
   !!   sbc_wave      : wave data from wave model in netcdf files 
   !!   sbc_wave_init : initialisation fo surface waves 
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants 
   USE oce            ! ocean variables
   USE sbc_oce	       ! Surface boundary condition: ocean fields
   USE zdf_oce,  ONLY : ln_zdfswm, ln_zdfst ! NB add second boolean
   USE bdy_oce        ! open boundary condition variables
   USE domvvl         ! domain: variable volume layers
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE fldread	       ! read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_stokes      ! routine called in sbccpl
   PUBLIC   sbc_wstress     ! routine called in sbcmod 
   PUBLIC   sbc_wave        ! routine called in sbcmod
   PUBLIC   sbc_wave_init   ! routine called in sbcmod
   
   ! Variables checking if the wave parameters are coupled (if not, they are read from file)
   LOGICAL, PUBLIC ::   cpl_hsig   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_phioc  = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrftx = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrfty = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wper   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wfreq  = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wnum   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_tauwoc = .FALSE.
   LOGICAL, PUBLIC ::   cpl_tauw   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wdrag  = .FALSE.
   !--- NB -------
   LOGICAL, PUBLIC ::   cpl_swell  = .FALSE.
   !--- END NB ---

   INTEGER ::   jpfld    ! number of files to read for stokes drift
   INTEGER ::   jp_usd   ! index of stokes drift  (i-component) (m/s)    at T-point
   INTEGER ::   jp_vsd   ! index of stokes drift  (j-component) (m/s)    at T-point
   INTEGER ::   jp_hsw   ! index of significant wave hight      (m)      at T-point
   INTEGER ::   jp_wmp   ! index of mean wave period            (s)      at T-point
   INTEGER ::   jp_wfr   ! index of wave peak frequency         (1/s)    at T-point
   !NB
   INTEGER ::   jp_wmd   ! index of mean wave direction         (s)      at T-point

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_cd      ! structure of input fields (file informations, fields read) Drag Coefficient
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sd      ! structure of input fields (file informations, fields read) Stokes Drift
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_wn      ! structure of input fields (file informations, fields read) wave number for Qiao
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tauwoc  ! structure of input fields (file informations, fields read) normalized wave stress into the ocean
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tauw    ! structure of input fields (file informations, fields read) ocean stress components from wave model
   !--- NB -------
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_swell   ! structure of input fields (file informations, fields read) swell characteristics from wave model   
   !--- END NB ---

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   cdn_wave            !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hsw, wmp, wnum      !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wfreq               !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wave          !:  
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauw_x, tauw_y      !:  
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tsd2d               !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   div_sd              !: barotropic stokes drift divergence
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   ut0sd, vt0sd        !: surface Stokes drift velocities at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   usd  , vsd  , wsd   !: Stokes drift velocities at u-, v- & w-points, resp.
   !--- NB -------
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   usd_b, vsd_b, wsd_b   !: Stokes drift velocities at u-, v- & w-points, resp. before
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wmd
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   swell_swh, swell_wmp, swell_wmd
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   swell_usd, swell_vsd
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   swell_dut_dz, swell_dvt_dz
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   wnd_dut_dz, wnd_dvt_dz
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   sw_tsd2d               !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   swell_ut, swell_vt          ! Surface Stokes drift (Mc Williams)

   !--- END NB ---

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcwave.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_stokes( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_stokes  ***
      !!
      !! ** Purpose :   compute the 3d Stokes Drift according to Breivik et al.,
      !!                2014 (DOI: 10.1175/JPO-D-14-0020.1)
      !!
      !! ** Method  : - Calculate Stokes transport speed 
      !!              - Calculate horizontal divergence 
      !!              - Integrate the horizontal divergenze from the bottom 
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      INTEGER  ::   jj, ji, jk   ! dummy loop argument
      INTEGER  ::   ik           ! local integer 
      REAL(wp) ::  ztransp, zfac, zsp0
      REAL(wp) ::  zdepth, zsqrt_depth,  zexp_depth, z_two_thirds, zsqrtpi !sqrt of pi
      REAL(wp) ::  zbot_u, zbot_v, zkb_u, zkb_v, zke3_u, zke3_v, zda_u, zda_v
      REAL(wp) ::  zstokes_psi_u_bot, zstokes_psi_v_bot
      REAL(wp) ::  zdep_u, zdep_v, zkh_u, zkh_v
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   zk_t, zk_u, zk_v, zu0_sd, zv0_sd     ! 2D workspace
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   zstokes_psi_u_top, zstokes_psi_v_top ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ze3divh                              ! 3D workspace
      !--- NB -------
      REAL(wp) ::  mydummy, fct_all, factor, zkh
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   sw_zk_t, sw_zk_u, sw_zk_v, sw_zu0_sd, sw_zv0_sd
      !--- END NB ---
      !!---------------------------------------------------------------------
      !
      ALLOCATE( ze3divh(jpi,jpj,jpk) )
      ALLOCATE( zk_t(jpi,jpj), zk_u(jpi,jpj), zk_v(jpi,jpj), zu0_sd(jpi,jpj), zv0_sd(jpi,jpj) )
      !--- NB -------
      wnd_dut_dz   = 0.; wnd_dvt_dz   = 0.; 
      ut0sd = 0._wp; vt0sd = 0._wp; 
      usd = 0._wp; vsd = 0._wp; wsd = 0._wp
      IF( ln_st_swell ) THEN
          ALLOCATE( sw_zk_t(jpi,jpj), sw_zk_u(jpi,jpj), sw_zk_v(jpi,jpj), sw_zu0_sd(jpi,jpj), sw_zv0_sd(jpi,jpj) )
          swell_dut_dz = 0.; swell_dvt_dz = 0.; 
          swell_ut = 0._wp; swell_vt = 0._wp
          swell_usd = 0._wp; swell_vsd = 0._wp;
      ENDIF
      !
      ! Compute surface stokes drift if not provided - next version would be coupled to WW3
      ! at T point as waves known at centre of elements

      !WRITE(numout,'(A,6F10.2)'),  'WW', MINVAL( hsw ), MAXVAL( hsw ), MINVAL( wmd ), MAXVAL( wmd ), MINVAL( wmp ), MAXVAL( wmp )
      !WRITE(numout,'(A,6F10.2)'),  'TS', MINVAL( swell_swh), MAXVAL( swell_swh), MINVAL( swell_wmd ), MAXVAL( swell_wmd ), MINVAL( swell_wmp ), MAXVAL( swell_wmp )
      IF (ln_compute_st) THEN 
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF ( gdept_n(ji,jj,mbkt(ji,jj)) < 10._wp ) hsw(ji,jj) = 0._wp
               CALL ESTIMATE_SURF_STK_VEC( hsw(ji,jj), wmp(ji,jj), wmd(ji,jj), ut0sd(ji,jj), vt0sd(ji,jj) )
               ut0sd(ji,jj) = ut0sd(ji,jj) * tmask(ji,jj,1)
               vt0sd(ji,jj) = vt0sd(ji,jj) * tmask(ji,jj,1) 
               IF ( gdept_n(ji,jj,mbkt(ji,jj)) < 10._wp ) ut0sd(ji,jj) = 0._wp
               IF ( gdept_n(ji,jj,mbkt(ji,jj)) < 10._wp ) vt0sd(ji,jj) = 0._wp

               ! COMPUTE STOKES SWELL
               IF( ln_st_swell ) THEN
                   IF ( gdept_n(ji,jj,mbkt(ji,jj)) < 10._wp ) swell_swh(ji,jj) = 0._wp
                   CALL ESTIMATE_SURF_STK_VEC( swell_swh(ji,jj), swell_wmp(ji,jj), swell_wmd(ji,jj), swell_ut(ji,jj), swell_vt(ji,jj) )
                   swell_ut(ji,jj) = swell_ut(ji,jj) * tmask(ji,jj,1)
                   swell_vt(ji,jj) = swell_vt(ji,jj) * tmask(ji,jj,1) 
                   IF ( gdept_n(ji,jj,mbkt(ji,jj)) < 10._wp )  swell_ut(ji,jj)= 0._wp
                   IF ( gdept_n(ji,jj,mbkt(ji,jj)) < 10._wp )  swell_vt(ji,jj)= 0._wp
               ENDIF

            ENDDO
         ENDDO
      ENDIF
      !--- END NB ---
      !
      ! select parameterization for the calculation of vertical Stokes drift
      ! exp. wave number at t-point
      IF( ll_st_bv_li ) THEN   ! (Eq. (19) in Breivik et al. (2014) )
         zfac = 2.0_wp * rpi / 16.0_wp
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Stokes drift velocity estimated from Hs and Tmean
               ztransp = zfac * hsw(ji,jj)*hsw(ji,jj) / MAX( wmp(ji,jj), 0.0000001_wp )
               ! Stokes surface speed
               tsd2d(ji,jj) = SQRT( ut0sd(ji,jj)*ut0sd(ji,jj) + vt0sd(ji,jj)*vt0sd(ji,jj))
               ! Wavenumber scale
               zk_t (ji,jj) = ABS( tsd2d(ji,jj) ) / MAX( ABS( 5.97_wp*ztransp ), 0.0000001_wp )

               ! COMPUTE STOKES SWELL
               IF( ln_st_swell ) THEN
                   ! Stokes drift velocity estimated from Hs and Tmean
                   ztransp = zfac * swell_swh(ji,jj)*swell_swh(ji,jj) / MAX( swell_wmp(ji,jj), 0.0000001_wp )
                   ! Stokes surface speed
                   sw_tsd2d(ji,jj) = SQRT( swell_ut(ji,jj)*swell_ut(ji,jj) + swell_vt(ji,jj)*swell_vt(ji,jj))
                   ! Wavenumber scale
                   sw_zk_t (ji,jj) = ABS( sw_tsd2d(ji,jj) ) / MAX( ABS( 5.97_wp*ztransp ), 0.0000001_wp )
               ENDIF

            END DO
         END DO
         DO jj = 1, jpjm1              ! exp. wave number & Stokes drift velocity at u- & v-points
            DO ji = 1, jpim1
               zk_u(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji+1,jj) )
               zk_v(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji,jj+1) )
               !
               zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
               zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )

               ! COMPUTE STOKES SWELL
               IF( ln_st_swell ) THEN
                  sw_zk_u(ji,jj) = 0.5_wp * ( sw_zk_t(ji,jj) + sw_zk_t(ji+1,jj) )
                  sw_zk_v(ji,jj) = 0.5_wp * ( sw_zk_t(ji,jj) + sw_zk_t(ji,jj+1) )
               !
                  sw_zu0_sd(ji,jj) = 0.5_wp * ( swell_ut(ji,jj) + swell_ut(ji+1,jj) )
                  sw_zv0_sd(ji,jj) = 0.5_wp * ( swell_vt(ji,jj) + swell_vt(ji,jj+1) )
               ENDIF

            END DO
         END DO
      ELSE IF( ll_st_peakfr ) THEN    ! peak wave number calculated from the peak frequency received by the wave model
         DO jj = 1, jpj
            DO ji = 1, jpi
               zk_t(ji,jj) = ( 2.0_wp * rpi * wfreq(ji,jj) ) * ( 2.0_wp * rpi * wfreq(ji,jj) ) / grav
            END DO
         END DO
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zk_u(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji+1,jj) )
               zk_v(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji,jj+1) )
               !
               zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
               zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )
            END DO
         END DO
      ENDIF
      !
      !                       !==  horizontal Stokes Drift 3D velocity  ==!
      IF( ll_st_bv2014 ) THEN
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zdep_u = 0.5_wp * ( gdept_n(ji,jj,jk) + gdept_n(ji+1,jj,jk) )
                  zdep_v = 0.5_wp * ( gdept_n(ji,jj,jk) + gdept_n(ji,jj+1,jk) )
                  !                          
                  zkh_u = zk_u(ji,jj) * zdep_u     ! k * depth
                  zkh_v = zk_v(ji,jj) * zdep_v
                  !                                ! Depth attenuation
                  zda_u = EXP( -2.0_wp*zkh_u ) / ( 1.0_wp + 8.0_wp*zkh_u )
                  zda_v = EXP( -2.0_wp*zkh_v ) / ( 1.0_wp + 8.0_wp*zkh_v )
                  !
                  usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
                  vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)

                  ! COMPUTE STOKES SWELL
                  IF( ln_st_swell ) THEN
                      zkh_u = sw_zk_u(ji,jj) * zdep_u     ! k * depth
                      zkh_v = sw_zk_v(ji,jj) * zdep_v
                  !
                      zda_u = EXP( -2.0_wp*zkh_u ) / ( 1.0_wp + 8.0_wp*zkh_u )
                      zda_v = EXP( -2.0_wp*zkh_v ) / ( 1.0_wp + 8.0_wp*zkh_v )
                  !
                      swell_usd(ji,jj,jk) = zda_u * sw_zu0_sd(ji,jj) * umask(ji,jj,jk)
                      swell_vsd(ji,jj,jk) = zda_v * sw_zv0_sd(ji,jj) * vmask(ji,jj,jk)
                  ENDIF

               END DO
            END DO
         END DO
      ELSE IF( ll_st_li2017 .OR. ll_st_peakfr ) THEN
         ALLOCATE( zstokes_psi_u_top(jpi,jpj), zstokes_psi_v_top(jpi,jpj) )
         DO jj = 1, jpjm1              ! exp. wave number & Stokes drift velocity at u- & v-points
            DO ji = 1, jpim1
               zstokes_psi_u_top(ji,jj) = 0._wp
               zstokes_psi_v_top(ji,jj) = 0._wp
            END DO
         END DO
         zsqrtpi = SQRT(rpi)
         z_two_thirds = 2.0_wp / 3.0_wp
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zbot_u = ( gdepw_n(ji,jj,jk+1) + gdepw_n(ji+1,jj,jk+1) )  ! 2 * bottom depth
                  zbot_v = ( gdepw_n(ji,jj,jk+1) + gdepw_n(ji,jj+1,jk+1) )  ! 2 * bottom depth
                  zkb_u  = zk_u(ji,jj) * zbot_u                             ! 2 * k * bottom depth
                  zkb_v  = zk_v(ji,jj) * zbot_v                             ! 2 * k * bottom depth
                  !
                  zke3_u = MAX(1.e-8_wp, 2.0_wp * zk_u(ji,jj) * e3u_n(ji,jj,jk))     ! 2k * thickness
                  zke3_v = MAX(1.e-8_wp, 2.0_wp * zk_v(ji,jj) * e3v_n(ji,jj,jk))     ! 2k * thickness

                  ! Depth attenuation .... do u component first..
                  zdepth      = zkb_u
                  zsqrt_depth = SQRT(zdepth)
                  zexp_depth  = EXP(-zdepth)
                  zstokes_psi_u_bot = 1.0_wp - zexp_depth  &
                       &              - z_two_thirds * ( zsqrtpi*zsqrt_depth*zdepth*ERFC(zsqrt_depth) &
                       &              + 1.0_wp - (1.0_wp + zdepth)*zexp_depth )
                  zda_u                    = ( zstokes_psi_u_bot - zstokes_psi_u_top(ji,jj) ) / zke3_u
                  zstokes_psi_u_top(ji,jj) =   zstokes_psi_u_bot

                  !         ... and then v component
                  zdepth      =zkb_v
                  zsqrt_depth = SQRT(zdepth)
                  zexp_depth  = EXP(-zdepth)
                  zstokes_psi_v_bot = 1.0_wp - zexp_depth  &
                       &              - z_two_thirds * ( zsqrtpi*zsqrt_depth*zdepth*ERFC(zsqrt_depth) &
                       &              + 1.0_wp - (1.0_wp + zdepth)*zexp_depth )
                  zda_v                    = ( zstokes_psi_v_bot - zstokes_psi_v_top(ji,jj) ) / zke3_v
                  zstokes_psi_v_top(ji,jj) =   zstokes_psi_v_bot
                  !
                  usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
                  vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)
                  ! NB - dont consider waves if not enough water depth
                  !IF ( gdept_n(ji,jj,mbkt(ji,jj)) <  swell_swh(ji,jj)/ 0.73_wp ) THEN
                  !     usd(ji,jj,jk)= 0._wp
                  !     vsd(ji,jj,jk)= 0._wp
                  !ENDIF
                  ! END NB
               END DO
            END DO
         END DO
         DEALLOCATE( zstokes_psi_u_top, zstokes_psi_v_top )
      ENDIF
      !--- NB -------
      ! Compute stokes drift vertical derivative using analytical formulation
      ! following Breivik formulation - not developed otherwise
      ! for GLS turbulence scheme at W points
      IF( ln_zdfst ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  IF (jk.eq.1) THEN
                     zkh = 0.0_wp
                  ELSE
                     zkh = - zk_t(ji,jj) * 0.5_wp *                        &
                         &   (gdept_n(ji,jj,jk-1)+gdept_n(ji,jj,jk) )          ! compute kz at t points (multiply by -1 to swith why axis)
                  ENDIF 

                  factor = 2.0_wp * ( 1.0_wp - 4.0_wp - 8.0_wp * zkh )     &
                         &        / ( 1.0_wp - 8.0_wp * zkh )**2.0_wp          ! compute factor : analytical derivative of Breivik
                  
                  fct_all = zk_t(ji,jj) * EXP( 2.0_wp * zkh) * factor
                  wnd_dut_dz(ji,jj,jk) = fct_all * zk_t(ji,jj) * ut0sd(ji,jj) * wmask(ji,jj,jk)    
                  wnd_dvt_dz(ji,jj,jk) = fct_all * zk_t(ji,jj) * vt0sd(ji,jj) * wmask(ji,jj,jk)

                  ! COMPUTE STOKES SWELL
                  IF( ln_st_swell ) THEN
                      swell_dut_dz(ji,jj,jk) = fct_all * sw_zk_t(ji,jj) * swell_ut(ji,jj) * wmask(ji,jj,jk)
                      swell_dvt_dz(ji,jj,jk) = fct_all * sw_zk_t(ji,jj) * swell_vt(ji,jj) * wmask(ji,jj,jk)
                  ENDIF

               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'sbcwave',   wnd_dut_dz, 'U', -1.,   wnd_dvt_dz, 'V', -1. )
         CALL lbc_lnk_multi( 'sbcwave', swell_dut_dz, 'U', -1., swell_dvt_dz, 'V', -1. )
      ENDIF
      !--- END NB ---
      !
      CALL lbc_lnk_multi( 'sbcwave',       usd, 'U', -1.,       vsd, 'V', -1. )
      IF( ln_st_swell ) CALL lbc_lnk_multi( 'sbcwave', swell_usd, 'U', -1., swell_vsd, 'V', -1. )
      !
      ! Save 3D horizontal wind wave part (what is computed originally)
      ! Save new swell component
      IF( ln_st_swell ) THEN
         CALL iom_put( "ustokes_wndsea",  usd  )
         CALL iom_put( "vstokes_wndsea",  vsd  )

         CALL iom_put( "ustokes_swell",  swell_usd  )
         CALL iom_put( "vstokes_swell",  swell_vsd  )
         !
         ! Combine them before computing vertical component
         ! to not copy twice the code
         ! if using 1 to jpm, then code crrashes with infinity
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  usd(ji,jj,jk) = (usd(ji,jj,jk) + swell_usd(ji,jj,jk)) !* umask(ji,jj,jk)
                  vsd(ji,jj,jk) = (vsd(ji,jj,jk) + swell_vsd(ji,jj,jk)) !* vmask(ji,jj,jk)
               END DO
            END DO
         END DO 
      !
         IF(lwp) WRITE(numout,*), "NB - STOKES DRIFT WW + SWELL"
      ENDIF
      !
      CALL lbc_lnk_multi( 'sbcwave',       usd, 'U', -1.,       vsd, 'V', -1. )
      !
      !--- END NB -----
      !
      !                       !==  vertical Stokes Drift 3D velocity  ==!
      !
      DO jk = 1, jpkm1               ! Horizontal e3*divergence
         DO jj = 2, jpj
            DO ji = fs_2, jpi
               ze3divh(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * usd(ji  ,jj,jk)    &
                  &                 - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * usd(ji-1,jj,jk)    &
                  &                 + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * vsd(ji,jj  ,jk)    &
                  &                 - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * vsd(ji,jj-1,jk)  ) * r1_e1e2t(ji,jj)
            END DO
         END DO
      END DO
      !
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
         IF( nbondi == -1 .OR. nbondi == 2 )   ze3divh( 2:nbghostcells+1,:        ,:) = 0._wp      ! west
         IF( nbondi ==  1 .OR. nbondi == 2 )   ze3divh( nlci-nbghostcells:nlci-1,:,:) = 0._wp      ! east
         IF( nbondj == -1 .OR. nbondj == 2 )   ze3divh( :,2:nbghostcells+1        ,:) = 0._wp      ! south
         IF( nbondj ==  1 .OR. nbondj == 2 )   ze3divh( :,nlcj-nbghostcells:nlcj-1,:) = 0._wp      ! north
      ENDIF
#endif
      !
      CALL lbc_lnk( 'sbcwave', ze3divh, 'T', 1. )
      !
      IF( ln_linssh ) THEN   ;   ik = 1   ! none zero velocity through the sea surface
      ELSE                   ;   ik = 2   ! w=0 at the surface (set one for all in sbc_wave_init)
      ENDIF
      DO jk = jpkm1, ik, -1          ! integrate from the bottom the hor. divergence (NB: at k=jpk w is always zero)
         wsd(:,:,jk) = wsd(:,:,jk+1) - ze3divh(:,:,jk)
      END DO
      !
      IF( ln_bdy ) THEN
         DO jk = 1, jpkm1
            wsd(:,:,jk) = wsd(:,:,jk) * bdytmask(:,:)
         END DO
      ENDIF
      !                       !==  Horizontal divergence of barotropic Stokes transport  ==!
      div_sd(:,:) = 0._wp
      DO jk = 1, jpkm1                                 ! 
        div_sd(:,:) = div_sd(:,:) + ze3divh(:,:,jk)
      END DO
      !
      CALL iom_put( "ustokes",  usd  )
      CALL iom_put( "vstokes",  vsd  )
      CALL iom_put( "wstokes",  wsd  )
      !
      DEALLOCATE( ze3divh )
      DEALLOCATE( zk_t, zk_u, zk_v, zu0_sd, zv0_sd )
      !
      !--- NB -------
      IF( ln_st_swell ) THEN
          DEALLOCATE( sw_zk_t, sw_zk_u, sw_zk_v, sw_zu0_sd, sw_zv0_sd )
      ENDIF
      !--- END NB ---
   END SUBROUTINE sbc_stokes


   SUBROUTINE sbc_wstress( )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wstress  ***
      !!
      !! ** Purpose :   Updates the ocean momentum modified by waves
      !!
      !! ** Method  : - Calculate u,v components of stress depending on stress
      !!                model 
      !!              - Calculate the stress module
      !!              - The wind module is not modified by waves 
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER  ::   jj, ji   ! dummy loop argument
      !
      IF( ln_tauwoc ) THEN
         utau(:,:) = utau(:,:)*tauoc_wave(:,:)
         vtau(:,:) = vtau(:,:)*tauoc_wave(:,:)
         taum(:,:) = taum(:,:)*tauoc_wave(:,:)
      ENDIF
      !
      IF( ln_tauw ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               ! Stress components at u- & v-points
               utau(ji,jj) = 0.5_wp * ( tauw_x(ji,jj) + tauw_x(ji+1,jj) )
               vtau(ji,jj) = 0.5_wp * ( tauw_y(ji,jj) + tauw_y(ji,jj+1) )
               !
               ! Stress module at t points
               taum(ji,jj) = SQRT( tauw_x(ji,jj)*tauw_x(ji,jj) + tauw_y(ji,jj)*tauw_y(ji,jj) )
            END DO
         END DO
         CALL lbc_lnk_multi( 'sbcwave', utau(:,:), 'U', -1. , vtau(:,:), 'V', -1. , taum(:,:) , 'T', -1. )
      ENDIF
      !
   END SUBROUTINE sbc_wstress


   SUBROUTINE sbc_wave( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave  ***
      !!
      !! ** Purpose :   read wave parameters from wave model  in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - Read Cd_n10 fields in netcdf files 
      !!              - Read stokes drift 2d in netcdf files 
      !!              - Read wave number in netcdf files 
      !!              - Compute 3d stokes drift using Breivik et al.,2014
      !!                formulation
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      !
      IF( ln_cdgw .AND. .NOT. cpl_wdrag ) THEN     !==  Neutral drag coefficient  ==!
         CALL fld_read( kt, nn_fsbc, sf_cd )             ! read from external forcing
         cdn_wave(:,:) = sf_cd(1)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      IF( ln_tauwoc .AND. .NOT. cpl_tauwoc ) THEN  !==  Wave induced stress  ==!
         CALL fld_read( kt, nn_fsbc, sf_tauwoc )         ! read wave norm stress from external forcing
         tauoc_wave(:,:) = sf_tauwoc(1)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      IF( ln_tauw .AND. .NOT. cpl_tauw ) THEN      !==  Wave induced stress  ==!
         CALL fld_read( kt, nn_fsbc, sf_tauw )           ! read ocean stress components from external forcing (T grid)
         tauw_x(:,:) = sf_tauw(1)%fnow(:,:,1) * tmask(:,:,1)
         tauw_y(:,:) = sf_tauw(2)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      !--- NB -------
      IF( ln_st_swell .AND. .NOT. cpl_swell ) THEN
         CALL fld_read( kt, nn_fsbc, sf_swell )          ! read swell characteristics
         swell_swh(:,:) = sf_swell(1)%fnow(:,:,1) * tmask(:,:,1)   ! significant wave height
         swell_wmp(:,:) = sf_swell(2)%fnow(:,:,1) * tmask(:,:,1)   ! wave mean period
         swell_wmd(:,:) = sf_swell(3)%fnow(:,:,1) * tmask(:,:,1)   ! wave mean direction
      ENDIF
      !--- END NB ---

      IF( ln_sdw )  THEN                           !==  Computation of the 3d Stokes Drift  ==! 
         !
         IF( jpfld > 0 ) THEN                            ! Read from file only if the field is not coupled
            CALL fld_read( kt, nn_fsbc, sf_sd )          ! read wave parameters from external forcing
            IF( jp_hsw > 0 )   hsw  (:,:) = sf_sd(jp_hsw)%fnow(:,:,1) * tmask(:,:,1)  ! significant wave height
            IF( jp_wmp > 0 )   wmp  (:,:) = sf_sd(jp_wmp)%fnow(:,:,1) * tmask(:,:,1)  ! wave mean period
            IF( jp_wfr > 0 )   wfreq(:,:) = sf_sd(jp_wfr)%fnow(:,:,1) * tmask(:,:,1)  ! Peak wave frequency
            IF( jp_usd > 0 )   ut0sd(:,:) = sf_sd(jp_usd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D zonal Stokes Drift at T point
            IF( jp_vsd > 0 )   vt0sd(:,:) = sf_sd(jp_vsd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D meridional Stokes Drift at T point
            IF( jp_wmd > 0 )   wmd  (:,:) = sf_sd(jp_wmd)%fnow(:,:,1) * tmask(:,:,1)  ! wave mean direction
         ENDIF
         !
         ! Read also wave number if needed, so that it is available in coupling routines
         IF( ln_zdfswm .AND. .NOT.cpl_wnum ) THEN
            CALL fld_read( kt, nn_fsbc, sf_wn )          ! read wave parameters from external forcing
            wnum(:,:) = sf_wn(1)%fnow(:,:,1) * tmask(:,:,1)
         ENDIF
           
         ! Calculate only if required fields have been read
         ! In coupled wave model-NEMO case the call is done after coupling
         !
         IF( ( ll_st_bv_li   .AND. jp_hsw>0 .AND. jp_wmp>0 .AND. jp_usd>0 .AND. jp_vsd>0 ) .OR. &
         !--- NB -------
           & ( ll_st_bv_li   .AND. ln_compute_st .AND. jp_hsw>0 .AND. jp_wmp>0 .AND. jp_wmd>0) .OR. &
         !--- END NB ---
           & ( ll_st_peakfr  .AND. jp_wfr>0 .AND. jp_usd>0 .AND. jp_vsd>0                ) ) CALL sbc_stokes( kt )
         !
      ENDIF
      !
   END SUBROUTINE sbc_wave


   SUBROUTINE sbc_wave_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave_init  ***
      !!
      !! ** Purpose :   read wave parameters from wave model  in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - Read Cd_n10 fields in netcdf files 
      !!              - Read stokes drift 2d in netcdf files 
      !!              - Read wave number in netcdf files 
      !!              - Compute 3d stokes drift using Breivik et al.,2014
      !!                formulation
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER ::   ierror, ios   ! local integer
      INTEGER ::   ifpr
      !!
      CHARACTER(len=100)     ::  cn_dir                            ! Root directory for location of drag coefficient files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   slf_i, slf_j     ! array of namelist informations on the fields to read
      TYPE(FLD_N)            ::  sn_cdg, sn_usd, sn_vsd,  &
                             &   sn_hsw, sn_wmp, sn_wfr, sn_wnum, &
                             &   sn_tauwoc, sn_tauwx, sn_tauwy     ! informations about the fields to be read
      !--- NB -------
      TYPE(FLD_N)            ::  sn_wmd
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   slf_s
      TYPE(FLD_N)            ::  sn_swell_swh, sn_swell_wmp,      &
                             &   sn_swell_wmd                      ! swell fields to be read
      !--- END NB ---
      !
      NAMELIST/namsbc_wave/  sn_cdg, cn_dir, sn_usd, sn_vsd, sn_hsw, sn_wmp, sn_wfr, &
                             sn_wnum, sn_tauwoc, sn_tauwx, sn_tauwy,                 &
      !--- NB -------
                             sn_swell_swh, sn_swell_wmp, sn_swell_wmd, sn_wmd
      !--- END NB ---
      !!---------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namsbc_wave in reference namelist : File for drag coeff. from wave model
      READ  ( numnam_ref, namsbc_wave, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_wave in reference namelist' )
         
      REWIND( numnam_cfg )              ! Namelist namsbc_wave in configuration namelist : File for drag coeff. from wave model
      READ  ( numnam_cfg, namsbc_wave, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_wave in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_wave )
      !
      IF( ln_cdgw ) THEN
         IF( .NOT. cpl_wdrag ) THEN
            ALLOCATE( sf_cd(1), STAT=ierror )               !* allocate and fill sf_wave with sn_cdg
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
                                   ALLOCATE( sf_cd(1)%fnow(jpi,jpj,1)   )
            IF( sn_cdg%ln_tint )   ALLOCATE( sf_cd(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_cd, (/ sn_cdg /), cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
         ENDIF
         ALLOCATE( cdn_wave(jpi,jpj) )
      ENDIF

      IF( ln_tauwoc ) THEN
         IF( .NOT. cpl_tauwoc ) THEN
            ALLOCATE( sf_tauwoc(1), STAT=ierror )           !* allocate and fill sf_wave with sn_tauwoc
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
                                     ALLOCATE( sf_tauwoc(1)%fnow(jpi,jpj,1)   )
            IF( sn_tauwoc%ln_tint )  ALLOCATE( sf_tauwoc(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_tauwoc, (/ sn_tauwoc /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( tauoc_wave(jpi,jpj) )
      ENDIF

      IF( ln_tauw ) THEN
         IF( .NOT. cpl_tauw ) THEN
            ALLOCATE( sf_tauw(2), STAT=ierror )           !* allocate and fill sf_wave with sn_tauwx/y
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_tauw structure' )
            !
            ALLOCATE( slf_j(2) )
            slf_j(1) = sn_tauwx
            slf_j(2) = sn_tauwy
                                    ALLOCATE( sf_tauw(1)%fnow(jpi,jpj,1)   )
                                    ALLOCATE( sf_tauw(2)%fnow(jpi,jpj,1)   )
            IF( slf_j(1)%ln_tint )  ALLOCATE( sf_tauw(1)%fdta(jpi,jpj,1,2) )
            IF( slf_j(2)%ln_tint )  ALLOCATE( sf_tauw(2)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_tauw, (/ slf_j /), cn_dir, 'sbc_wave_init', 'read wave input', 'namsbc_wave' )
         ENDIF
         ALLOCATE( tauw_x(jpi,jpj) )
         ALLOCATE( tauw_y(jpi,jpj) )
      ENDIF
 
      !--- NB -------
      ! option to read a second set of wave parameters for SWELL 
      ! to compute its stokes drift according to McWilliams (2014)
      IF( ln_st_swell ) THEN
         IF( .NOT. cpl_swell ) THEN
            ALLOCATE( sf_swell(3), STAT=ierror )           !* allocate and fill sf_wave with sn_swell
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_swell structure' )
            !
            ALLOCATE( slf_s(3) )
            slf_s(1) = sn_swell_swh
            slf_s(2) = sn_swell_wmp
            slf_s(3) = sn_swell_wmd
                                    ALLOCATE( sf_swell(1)%fnow(jpi,jpj,1)   )
                                    ALLOCATE( sf_swell(2)%fnow(jpi,jpj,1)   )
                                    ALLOCATE( sf_swell(3)%fnow(jpi,jpj,1)   )
            IF( slf_s(1)%ln_tint )  ALLOCATE( sf_swell(1)%fdta(jpi,jpj,1,2) )
            IF( slf_s(2)%ln_tint )  ALLOCATE( sf_swell(2)%fdta(jpi,jpj,1,2) )
            IF( slf_s(3)%ln_tint )  ALLOCATE( sf_swell(3)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_swell, (/ slf_s /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( swell_swh(jpi,jpj), swell_wmp(jpi,jpj), swell_wmd(jpi,jpj) )
         ALLOCATE( swell_usd(jpi,jpj,jpk), swell_vsd(jpi,jpj,jpk))
         ALLOCATE( swell_dut_dz(jpi,jpj,jpk), swell_dvt_dz(jpi,jpj,jpk) )
         ALLOCATE( sw_tsd2d (jpi,jpj) )  !NB
         ALLOCATE( swell_ut(jpi,jpj), swell_vt(jpi,jpj) )
         
         swell_swh(:,:) = 0._wp
         swell_wmd(:,:) = 0._wp
         swell_wmp(:,:) = 0._wp

         swell_usd(:,:,:) = 0._wp
         swell_vsd(:,:,:) = 0._wp

         swell_dut_dz(:,:,:) = 0._wp
         swell_dvt_dz(:,:,:) = 0._wp

         swell_ut(:,:) = 0._wp
         swell_vt(:,:) = 0._wp
         sw_tsd2d(:,:) = 0._wp
      ENDIF 
      !--- END NB ---

      IF( ln_sdw ) THEN   ! Find out how many fields have to be read from file if not coupled
         jpfld=0
         jp_usd=0   ;   jp_vsd=0   ;   jp_hsw=0   ;   jp_wmp=0   ;   jp_wfr=0 
         jp_wmd = 0 ! NB
         ! NB IF( .NOT. cpl_sdrftx ) THEN
         IF( .NOT. cpl_sdrftx .AND. .NOT. ln_compute_st) THEN
            jpfld  = jpfld + 1
            jp_usd = jpfld
         ENDIF
         ! NB IF( .NOT. cpl_sdrfty ) THEN
         IF( .NOT. cpl_sdrfty .AND. .NOT. ln_compute_st) THEN
            jpfld  = jpfld + 1
            jp_vsd = jpfld
         ENDIF
         IF( .NOT. cpl_hsig  .AND. ll_st_bv_li  ) THEN
            jpfld  = jpfld + 1
            jp_hsw = jpfld
         ENDIF
         IF( .NOT. cpl_wper  .AND. ll_st_bv_li  ) THEN
            jpfld  = jpfld + 1
            jp_wmp = jpfld
         ENDIF
         IF( .NOT. cpl_wfreq .AND. ll_st_peakfr ) THEN
            jpfld  = jpfld + 1
            jp_wfr = jpfld
         ENDIF
         !--- NB -------
         IF ( ln_compute_st ) THEN
            jpfld  = jpfld + 1
            jp_wmd = jpfld
         ENDIF
         !--- END NB ---

         ! Read from file only the non-coupled fields 
         IF( jpfld > 0 ) THEN
            ALLOCATE( slf_i(jpfld) )
            IF( jp_usd > 0 )   slf_i(jp_usd) = sn_usd
            IF( jp_vsd > 0 )   slf_i(jp_vsd) = sn_vsd
            IF( jp_hsw > 0 )   slf_i(jp_hsw) = sn_hsw
            IF( jp_wmp > 0 )   slf_i(jp_wmp) = sn_wmp
            IF( jp_wfr > 0 )   slf_i(jp_wfr) = sn_wfr
            IF( jp_wmd > 0 )   slf_i(jp_wmd) = sn_wmd   !--- NB add direction

            ALLOCATE( sf_sd(jpfld), STAT=ierror )   !* allocate and fill sf_sd with stokes drift
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
            DO ifpr= 1, jpfld
               ALLOCATE( sf_sd(ifpr)%fnow(jpi,jpj,1) )
               IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf_sd(ifpr)%fdta(jpi,jpj,1,2) )
            END DO
            !
            CALL fld_fill( sf_sd, slf_i, cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
         ENDIF
         ALLOCATE( usd  (jpi,jpj,jpk), vsd  (jpi,jpj,jpk), wsd(jpi,jpj,jpk) )
         ALLOCATE( hsw  (jpi,jpj)    , wmp  (jpi,jpj)     )
         ALLOCATE( wmd  (jpi,jpj) )  !--- NB - add direction
         ALLOCATE( wnd_dut_dz(jpi,jpj,jpk), wnd_dvt_dz(jpi,jpj,jpk) )

         ALLOCATE( wfreq(jpi,jpj) )
         ALLOCATE( ut0sd(jpi,jpj)    , vt0sd(jpi,jpj)     )
         ALLOCATE( div_sd(jpi,jpj) )
         ALLOCATE( tsd2d (jpi,jpj) )
         ALLOCATE( usd_b(jpi,jpj,jpk), vsd_b(jpi,jpj,jpk), wsd_b(jpi,jpj,jpk) )  !--- NB - keep previous stokes drift

         ut0sd(:,:) = 0._wp
         vt0sd(:,:) = 0._wp
         hsw  (:,:) = 0._wp
         wmp  (:,:) = 0._wp
         wmd  (:,:) = 0._wp !--- NB

         usd(:,:,:) = 0._wp
         vsd(:,:,:) = 0._wp
         wsd(:,:,:) = 0._wp

         !--- NB - keep previous stokes drift
         usd_b(:,:,:) = 0._wp
         vsd_b(:,:,:) = 0._wp
         wsd_b(:,:,:) = 0._wp

         wnd_dut_dz(:,:,:) = 0._wp
         wnd_dvt_dz(:,:,:) = 0._wp

         ! Wave number needed only if ln_zdfswm=T
         IF( .NOT. cpl_wnum ) THEN
            ALLOCATE( sf_wn(1), STAT=ierror )           !* allocate and fill sf_wave with sn_wnum
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable toallocate sf_wave structure' )
                                   ALLOCATE( sf_wn(1)%fnow(jpi,jpj,1)   )
            IF( sn_wnum%ln_tint )  ALLOCATE( sf_wn(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_wn, (/ sn_wnum /), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( wnum(jpi,jpj) )
      ENDIF
      !
   END SUBROUTINE sbc_wave_init

   !!======================================================================

   !--- NB -------
   ! Compute surface stokes drift as monochromatic waves
   ! Equation 14 of Webb and Fox-Kemper (2011)
   ! we assume Hs == Hm0 and monochromatic Tp = T0 = T1, etc
   FUNCTION ESTIMATE_SURF_STK( Hm0, Tm ) RESULT( STK )
      REAL(wp) :: Hm0, Tm , STK
      REAL(wp) :: Tm3, H2, pi3
      pi3 = rpi * rpi * rpi
      H2  = Hm0 * Hm0
      Tm3 = Tm * Tm * Tm 
      STK = pi3 * H2 / grav / Tm3
      IF (Hm0 .LE. 0.0_wp) STK = 0.0000001_wp
   END FUNCTION ESTIMATE_SURF_STK

   SUBROUTINE ESTIMATE_SURF_STK_VEC( Hm0, Tm, Dm, st_u, st_v )
      REAL(wp), INTENT(INOUT) :: Tm , Dm
      REAL(wp), INTENT(INOUT) :: Hm0
      REAL(wp), INTENT(INOUT) :: st_u, st_v
      REAL(wp) :: STK
      IF (Hm0 .LE. 0.0_wp .OR. Tm .LE. 0.0_wp)  THEN
         Hm0 = 0.0_wp
         Tm  = 0.0000001_wp
      ENDIF
      IF ( Tm .LE. 2.5_wp ) Tm = 2.5_wp 
      STK  = ESTIMATE_SURF_STK( Hm0, Tm )
      st_u = STK * COS( (270.0_wp-Dm) * rad )
      st_v = STK * SIN( (270.0_wp-Dm) * rad )
   END SUBROUTINE ESTIMATE_SURF_STK_VEC
   !--- END NB ---


END MODULE sbcwave
