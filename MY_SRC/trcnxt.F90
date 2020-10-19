MODULE trcnxt
   !!======================================================================
   !!                       ***  MODULE  trcnxt  ***
   !! Ocean passive tracers:  time stepping on passives tracers
   !!======================================================================
   !! History :  7.0  !  1991-11  (G. Madec)  Original code
   !!                 !  1993-03  (M. Guyon)  symetrical conditions
   !!                 !  1995-02  (M. Levy)   passive tracers
   !!                 !  1996-02  (G. Madec & M. Imbard)  opa release 8.0
   !!            8.0  !  1996-04  (A. Weaver)  Euler forward step
   !!            8.2  !  1999-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!                 !  2002-08  (G. Madec)  F90: Free form and module
   !!                 !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!                 !  2004-03  (C. Ethe) passive tracers
   !!                 !  2007-02  (C. Deltel) Diagnose ML trends for passive tracers
   !!            2.0  !  2006-02  (L. Debreu, C. Mazauric) Agrif implementation
   !!            3.0  !  2008-06  (G. Madec)  time stepping always done in trazdf
   !!            3.1  !  2009-02  (G. Madec, R. Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_nxt     : time stepping on passive tracers
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables
   USE trd_oce
   USE trdtra
   USE tranxt
   USE bdy_oce   , ONLY: ln_bdy
   USE trcbdy          ! BDY open boundaries
# if defined key_agrif
   USE agrif_top_interp
# endif
   !
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nxt   ! routine called by step.F90

   REAL(wp) ::   rfact1, rfact2

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnxt.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_nxt( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trcnxt  ***
      !!
      !! ** Purpose :   Compute the passive tracers fields at the 
      !!      next time-step from their temporal trends and swap the fields.
      !! 
      !! ** Method  :   Apply lateral boundary conditions on (ua,va) through 
      !!      call to lbc_lnk routine
      !!   default:
      !!      arrays swap
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!         (trb) = (trn) 
      !!
      !!   For Arakawa or TVD Scheme : 
      !!      A Asselin time filter applied on now tracers (trn) to avoid
      !!      the divergence of two consecutive time-steps and tr arrays
      !!      to prepare the next time_step:
      !!         (trb) = (trn) + atfp [ (trb) + (tra) - 2 (trn) ]
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!
      !!
      !! ** Action  : - update trb, trn
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
      !
      INTEGER  ::   jk, jn   ! dummy loop indices
      REAL(wp) ::   zfact            ! temporary scalar
      CHARACTER (len=22) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ztrdt    ! 4D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_nxt')
      !
      IF( kt == nittrc000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nxt : time stepping on passive tracers'
      ENDIF
      !
#if defined key_agrif
      CALL Agrif_trc                   ! AGRIF zoom boundaries
#endif
      ! Update after tracer on domain lateral boundaries
      CALL lbc_lnk( 'trcnxt', tra(:,:,:,:), 'T', 1. )   

!      IF( ln_bdy )  CALL trc_bdy( kt )
!--- NB -------
! hardcoded 0 tracer at boundaries
!
      igrd = 1                          ! Everything is at T-points here
      DO ib_bdy = 1, nb_bdy               ! Loop over open boundary only
         DO jn = 1, jptra
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               DO ik = 1, jpkm1
                  tra(ii,ij,ik,jn) = 0.  !1. * tmask(ii,ij,ik)
               ENDDO
            END DO
         END DO
      END DO
!
!--- END NB ---


      IF( l_trdtrc )  THEN             ! trends: store now fields before the Asselin filter application
         ALLOCATE( ztrdt(jpi,jpj,jpk,jptra) )
         ztrdt(:,:,:,:)  = 0._wp
         IF( ln_traldf_iso ) THEN                       ! diagnose the "pure" Kz diffusive trend 
            DO jn = 1, jptra
               CALL trd_tra( kt, 'TRC', jn, jptra_zdfp, ztrdt(:,:,:,jn) )
            ENDDO
         ENDIF

         ! total trend for the non-time-filtered variables. 
         zfact = 1.0 / rdttrc
         ! G Nurser 23 Mar 2017. Recalculate trend as Delta(e3t*T)/e3tn; e3tn cancel from tsn terms
         IF( ln_linssh ) THEN       ! linear sea surface height only
            DO jn = 1, jptra
               DO jk = 1, jpkm1
                  ztrdt(:,:,jk,jn) = ( tra(:,:,jk,jn)*e3t_a(:,:,jk) / e3t_n(:,:,jk) - trn(:,:,jk,jn)) * zfact
               END DO
            END DO
         ELSE
            DO jn = 1, jptra
               DO jk = 1, jpkm1
                  ztrdt(:,:,jk,jn) = ( tra(:,:,jk,jn) - trn(:,:,jk,jn) ) * zfact
               END DO
            END DO
         ENDIF
         !
         DO jn = 1, jptra
            CALL trd_tra( kt, 'TRC', jn, jptra_tot, ztrdt(:,:,:,jn) )
         ENDDO
         !
         IF( ln_linssh ) THEN       ! linear sea surface height only
            ! Store now fields before applying the Asselin filter 
            ! in order to calculate Asselin filter trend later.
            ztrdt(:,:,:,:) = trn(:,:,:,:) 
         ENDIF

      ENDIF
      !                                ! Leap-Frog + Asselin filter time stepping
      IF( (neuler == 0 .AND. kt == nittrc000) .OR. ln_top_euler ) THEN    ! Euler time-stepping (only swap)
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               trn(:,:,jk,jn) = tra(:,:,jk,jn)
               trb(:,:,jk,jn) = trn(:,:,jk,jn)  
            END DO
         END DO
         IF (l_trdtrc .AND. .NOT. ln_linssh ) THEN   ! Zero Asselin filter contribution must be explicitly written out since for vvl
            !                                        ! Asselin filter is output by tra_nxt_vvl that is not called on this time step
            ztrdt(:,:,:,:) = 0._wp            
            DO jn = 1, jptra
               CALL trd_tra( kt, 'TRC', jn, jptra_atf, ztrdt(:,:,:,jn) )
            ENDDO
         END IF
         !
      ELSE     
         IF( .NOT. l_offline ) THEN ! Leap-Frog + Asselin filter time stepping
            IF( ln_linssh ) THEN   ;   CALL tra_nxt_fix( kt, nittrc000,         'TRC', trb, trn, tra, jptra )  !     linear ssh
            ELSE                   ;   CALL tra_nxt_vvl( kt, nittrc000, rdttrc, 'TRC', trb, trn, tra,      &
              &                                                                   sbc_trc, sbc_trc_b, jptra )  ! non-linear ssh
            ENDIF
         ELSE
                                       CALL trc_nxt_off( kt )       ! offline 
         ENDIF
         !
         CALL lbc_lnk_multi( 'trcnxt', trb(:,:,:,:), 'T', 1._wp, trn(:,:,:,:), 'T', 1._wp, tra(:,:,:,:), 'T', 1._wp )
      ENDIF
      !
      IF( l_trdtrc .AND. ln_linssh ) THEN      ! trend of the Asselin filter (tb filtered - tb)/dt )
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               zfact = 1._wp / r2dttrc  
               ztrdt(:,:,jk,jn) = ( trb(:,:,jk,jn) - ztrdt(:,:,jk,jn) ) * zfact 
            END DO
            CALL trd_tra( kt, 'TRC', jn, jptra_atf, ztrdt(:,:,:,jn) )
         END DO
      END IF
      IF( l_trdtrc ) DEALLOCATE( ztrdt ) 
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nxt')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_nxt')
      !
   END SUBROUTINE trc_nxt


   SUBROUTINE trc_nxt_off( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_vvl  ***
      !!
      !! ** Purpose :   Time varying volume: apply the Asselin time filter  
      !!                and swap the tracer fields.
      !! 
      !! ** Method  : - Apply a thickness weighted Asselin time filter on now fields.
      !!              - save in (ta,sa) a thickness weighted average over the three 
      !!             time levels which will be used to compute rdn and thus the semi-
      !!             implicit hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = ( e3t_n*tn + rbcp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )   ln_dynhpg_imp = T
      !!                  /( e3t_n    + rbcp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )   
      !!             ztm = 0                                                       otherwise
      !!             tb  = ( e3t_n*tn + atfp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )
      !!                  /( e3t_n    + atfp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )
      !!             tn  = ta 
      !!             ta  = zt        (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )   ::  kt       ! ocean time-step index
      !!     
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices
      REAL(wp) ::   ztc_a , ztc_n , ztc_b , ztc_f , ztc_d    ! local scalar
      REAL(wp) ::   ze3t_b, ze3t_n, ze3t_a, ze3t_f, ze3t_d   !   -      -
      !!----------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_nxt_off : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF( .NOT. ln_linssh ) THEN
            rfact1 = atfp * rdttrc
            rfact2 = rfact1 / rau0
         ENDIF
        !  
      ENDIF
      !
      DO jn = 1, jptra      
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ze3t_b = e3t_b(ji,jj,jk)
                  ze3t_n = e3t_n(ji,jj,jk)
                  ze3t_a = e3t_a(ji,jj,jk)
                  !                                         ! tracer content at Before, now and after
                  ztc_b  = trb(ji,jj,jk,jn) * ze3t_b
                  ztc_n  = trn(ji,jj,jk,jn) * ze3t_n
                  ztc_a  = tra(ji,jj,jk,jn) * ze3t_a
                  !
                  ze3t_d = ze3t_a - 2. * ze3t_n + ze3t_b
                  ztc_d  = ztc_a  - 2. * ztc_n  + ztc_b
                  !
                  ze3t_f = ze3t_n + atfp * ze3t_d
                  ztc_f  = ztc_n  + atfp * ztc_d
                  !
                  IF( .NOT. ln_linssh .AND. jk == mikt(ji,jj) ) THEN           ! first level 
                     ze3t_f = ze3t_f - rfact2 * ( emp_b(ji,jj)      - emp(ji,jj)   ) 
                     ztc_f  = ztc_f  - rfact1 * ( sbc_trc(ji,jj,jn) - sbc_trc_b(ji,jj,jn) )
                  ENDIF

                  ze3t_f = 1.e0 / ze3t_f
                  trb(ji,jj,jk,jn) = ztc_f * ze3t_f       ! ptb <-- ptn filtered
                  trn(ji,jj,jk,jn) = tra(ji,jj,jk,jn)     ! ptn <-- pta
                  !
               END DO
            END DO
         END DO
         ! 
      END DO
      !
   END SUBROUTINE trc_nxt_off

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nxt( kt )  
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_nxt: You should not have seen this print! error?', kt
   END SUBROUTINE trc_nxt
#endif
   !!======================================================================
END MODULE trcnxt
