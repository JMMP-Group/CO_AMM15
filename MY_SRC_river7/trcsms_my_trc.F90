MODULE trcsms_my_trc
   !!======================================================================
   !!                         ***  MODULE trcsms_my_trc  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_sms_my_trc       : MY_TRC model main routine
   !! trc_sms_my_trc_alloc : allocate arrays specific to MY_TRC sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc
   USE trcbc, only : trc_bc
!-- NB load open boundary variables
   USE bdy_oce
!-- END NB ------------------------

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_my_trc       ! called by trcsms.F90 module
   PUBLIC   trc_sms_my_trc_alloc ! called by trcini_my_trc.F90 module

   ! Defined HERE the arrays specific to MY_TRC sms and ALLOCATE them in trc_sms_my_trc_alloc
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: tra_bkp ! keep in memory the initial conditions

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_my_trc.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_my_trc( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_my_trc  ***
      !!
      !! ** Purpose :   main routine of MY_TRC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::   jn   ! dummy loop index
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrmyt

      !-- NB
      INTEGER :: jgrd, jb_bdy,jb,jj,ji,jk
      REAL(wp):: uflg, vflg
      !-- END NB

      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_my_trc')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_my_trc:  MY_TRC model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      IF( l_trdtrc )  ALLOCATE( ztrmyt(jpi,jpj,jpk) )

      CALL trc_bc ( kt )       ! tracers: surface and lateral Boundary Conditions

      ! add here the call to BGC model
      !-- NB ------------------------------------------------- 
      ! Add river mass flux at each river mouth
      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jk = 1,nk_rnf(ji,jj)
               WHERE ( flagRIV(ji,jj,jk,1:jp_bgc ) .GT. 0. ) trn (ji,jj,jk,1:jp_bgc) = trn (ji,jj,jk,1:jp_bgc) + rnf(ji,jj) / h_rnf(ji,jj)
            ENDDO
         ENDDO
      ENDDO

      ! if other point, add mass flux based on velocity
      ! for baltic
      jgrd   = 2                        ! Everything is at U-points here
      jb_bdy = 2                        ! baltic sea bnd
      DO jn  = 2, 2                     ! baltic tracer
         DO jb = 1, idx_bdy(jb_bdy)%nblenrim(jgrd)
            ji   = idx_bdy(jb_bdy)%nbi  (jb,jgrd)
            jj   = idx_bdy(jb_bdy)%nbj  (jb,jgrd)
            uflg = idx_bdy(jb_bdy)%flagu(jb,jgrd)
            vflg = idx_bdy(jb_bdy)%flagv(jb,jgrd)
            DO jk = 1, jpkm1
               trn (ji,jj,jk,jn) = trn(ji,jj,jk,jn)     &
                          &      + uflg * (un(ji,jj,jk)+un(ji-1,jj,jk))*0.5 * rau0 / e2t(ji,jj) &
                          &      + vflg * (vn(ji,jj,jk)+vn(ji,jj-1,jk))*0.5 * rau0 / e1t(ji,jj)
               if ( trn (ji,jj,jk,jn) .lt. 0. ) trn (ji,jj,jk,jn) = 0.
            ENDDO
           
         END DO
      END DO

      jj=1

      !-- END NB ----------------------------------------------


      ! Save the trends in the mixed layer
      IF( l_trdtrc ) THEN
          DO jn = jp_myt0, jp_myt1
            ztrmyt(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrmyt, jn, jptra_sms, kt )   ! save trends
          END DO
          DEALLOCATE( ztrmyt )
      END IF
      !
      IF( ln_timing )   CALL timing_stop('trc_sms_my_trc')
      !
   END SUBROUTINE trc_sms_my_trc


   INTEGER FUNCTION trc_sms_my_trc_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_my_trc_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to MY_TRC
      ! ALLOCATE( tab(...) , STAT=trc_sms_my_trc_alloc )
      !-- NB initialise array
!      ALLOCATE( tra_bkp(jpi,jpj,jpk,jp_myt1-jp_myt0+1), STAT=trc_sms_my_trc_alloc )
!      IF( trc_sms_my_trc_alloc /= 0 ) CALL ctl_warn( 'trc_sms_my_trc_alloc : failed to allocate arrays tra 1' )
!      tra_bkp(:,:,:,:) = -999.0
      !-- END NB ------------
      !
      trc_sms_my_trc_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_my_trc_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_my_trc_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_my_trc_alloc

   !!======================================================================
END MODULE trcsms_my_trc
