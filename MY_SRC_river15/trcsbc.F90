MODULE trcsbc
   !!==============================================================================
   !!                       ***  MODULE  trcsbc  ***
   !! Ocean passive tracers:  surface boundary condition
   !!======================================================================
   !! History :  8.2  !  1998-10  (G. Madec, G. Roullet, M. Imbard)  Original code
   !!            8.2  !  2001-02  (D. Ludicone)  sea ice and free surface
   !!            8.5  !  2002-06  (G. Madec)  F90: Free form and module
   !!            9.0  !  2004-03  (C. Ethe)  adapted for passive tracers
   !!                 !  2006-08  (C. Deltel) Diagnose ML trends for passive tracers
   !!==============================================================================
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_sbc      : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean  passive tracers variables
   USE prtctl_trc      ! Print control for debbuging
   USE iom
   USE trd_oce
   USE trdtra

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sbc   ! routine called by step.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsbc.F90 10788 2019-03-21 11:15:14Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sbc ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sbc  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (concentration/dilution effect) and add it to the general 
      !!       trend of tracer equations.
      !!
      !! ** Method :
      !!      * concentration/dilution effect:
      !!            The surface freshwater flux modify the ocean volume
      !!         and thus the concentration of a tracer as :
      !!            tra = tra + emp * trn / e3t   for k=1
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation ) given in kg/m2/s is divided
      !!         by 1035 kg/m3 (density of ocean water) to obtain m/s.
      !!
      !! ** Action  : - Update the 1st level of tra with the trend associated
      !!                with the tracer surface boundary condition 
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jn                      ! dummy loop indices
      REAL(wp) ::   zse3t, zrtrn, zfact     ! local scalars
      REAL(wp) ::   zftra, zdtra, ztfx, ztra   !   -      -
      CHARACTER (len=22) :: charout
      REAL(wp), DIMENSION(jpi,jpj)   ::   zsfx
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrtrd
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sbc')
      !
      ! Allocate temporary workspace
      IF( l_trdtrc )  ALLOCATE( ztrtrd(jpi,jpj,jpk) )
      !
      zrtrn = 1.e-15_wp

      !-- NB -- INITIALISE TABLES THAT WERE NOT
      sbc_trc_b(:,:,:) = 0._wp
      sbc_trc  (:,:,:) = 0._wp
      !-- END NB ------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_sbc : Passive tracers surface boundary condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         !
         IF( ln_rsttr .AND. .NOT.ln_top_euler .AND.   &                     ! Restart: read in restart  file
            iom_varid( numrtr, 'sbc_'//TRIM(ctrcnm(1))//'_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nittrc000-nn_dttrc surface tracer content forcing fields red in the restart file'
            zfact = 0.5_wp
            DO jn = 1, jptra
               CALL iom_get( numrtr, jpdom_autoglo, 'sbc_'//TRIM(ctrcnm(jn))//'_b', sbc_trc_b(:,:,jn) )   ! before tracer content sbc
            END DO
         ELSE                                         ! No restart or restart not found: Euler forward time stepping
           zfact = 1._wp
           sbc_trc_b(:,:,:) = 0._wp
         ENDIF
      ELSE                                         ! Swap of forcing fields
         IF( ln_top_euler ) THEN
            zfact = 1._wp
            sbc_trc_b(:,:,:) = 0._wp
         ELSE
            zfact = 0.5_wp
            sbc_trc_b(:,:,:) = sbc_trc(:,:,:)
         ENDIF
         !
      ENDIF

      ! Coupling online : river runoff is added to the horizontal divergence (hdivn) in the subroutine sbc_rnf_div 
      ! one only consider the concentration/dilution effect due to evaporation minus precipitation + freezing/melting of sea-ice
      ! Coupling offline : runoff are in emp which contains E-P-R
      !
      IF( .NOT.ln_linssh ) THEN  ! online coupling with vvl
         zsfx(:,:) = 0._wp
      ELSE                                      ! online coupling free surface or offline with free surface
         zsfx(:,:) = emp(:,:)
      ENDIF

      ! 0. initialization
      SELECT CASE ( nn_ice_tr )

      CASE ( -1 ) ! No tracers in sea ice (null concentration in sea ice)
         !
         DO jn = 1, jptra
            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  sbc_trc(ji,jj,jn) = zsfx(ji,jj) * r1_rau0 * trn(ji,jj,1,jn)
               END DO
            END DO
         END DO
         !
      CASE ( 0 )  ! Same concentration in sea ice and in the ocean
         !
         DO jn = 1, jptra
            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  sbc_trc(ji,jj,jn) = ( zsfx(ji,jj) + fmmflx(ji,jj) ) * r1_rau0 * trn(ji,jj,1,jn)
               END DO
            END DO
         END DO
         !
      CASE ( 1 )  ! Specific treatment of sea ice fluxes with an imposed concentration in sea ice 
         !
         DO jn = 1, jptra
            DO jj = 2, jpj
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zse3t = 1. / e3t_n(ji,jj,1)
                  ! tracer flux at the ice/ocean interface (tracer/m2/s)
                  zftra = - trc_i(ji,jj,jn) * fmmflx(ji,jj) ! uptake of tracer in the sea ice
                  !                                         ! only used in the levitating sea ice case
                  ! tracer flux only       : add concentration dilution term in net tracer flux, no F-M in volume flux
                  ! tracer and mass fluxes : no concentration dilution term in net tracer flux, F-M term in volume flux
                  ztfx  = zftra                        ! net tracer flux
                  !
                  zdtra = r1_rau0 * ( ztfx + ( zsfx(ji,jj) + fmmflx(ji,jj) ) * trn(ji,jj,1,jn) ) 
                  IF ( zdtra < 0. ) THEN
                     zdtra  = MAX(zdtra, -trn(ji,jj,1,jn) * e3t_n(ji,jj,1) / r2dttrc )   ! avoid negative concentrations to arise
                  ENDIF
                  sbc_trc(ji,jj,jn) =  zdtra 
               END DO
            END DO
         END DO
      END SELECT
      !
      CALL lbc_lnk( 'trcsbc', sbc_trc(:,:,:), 'T', 1. )
      !                                       Concentration dilution effect on tracers due to evaporation & precipitation 
      DO jn = 1, jptra
         !
         IF( l_trdtrc )   ztrtrd(:,:,:) = tra(:,:,:,jn)  ! save trends
         !
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zse3t = zfact / e3t_n(ji,jj,1)
               tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ( sbc_trc_b(ji,jj,jn) + sbc_trc(ji,jj,jn) ) * zse3t
            END DO
         END DO
         !
         IF( l_trdtrc ) THEN
            ztrtrd(:,:,:) = tra(:,:,:,jn) - ztrtrd(:,:,:)
            CALL trd_tra( kt, 'TRC', jn, jptra_nsr, ztrtrd )
         END IF
         !                                                       ! ===========
      END DO                                                     ! tracer loop
      !                                                          ! ===========
      !
      !                                           Write in the tracer restar  file
      !                                          *******************************
      IF( lrst_trc .AND. .NOT.ln_top_euler ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : ocean surface tracer content forcing fields written in tracer restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         DO jn = 1, jptra
            CALL iom_rstput( kt, nitrst, numrtw, 'sbc_'//TRIM(ctrcnm(jn))//'_b', sbc_trc(:,:,jn) )
         END DO
      ENDIF
      !
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('sbc ')") ;  CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      ENDIF
      IF( l_trdtrc )  DEALLOCATE( ztrtrd )
      !
      IF( ln_timing )   CALL timing_stop('trc_sbc')
      !
   END SUBROUTINE trc_sbc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sbc (kt)              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_sbc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sbc
#endif
   
   !!======================================================================
END MODULE trcsbc
