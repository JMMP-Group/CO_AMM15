MODULE icethd_pnd 
   !!======================================================================
   !!                     ***  MODULE  icethd_pnd   ***
   !!   sea-ice: Melt ponds on top of sea ice
   !!======================================================================
   !! history :       !  2012     (O. Lecomte)       Adaptation from Flocco and Turner
   !!                 !  2017     (M. Vancoppenolle, O. Lecomte, C. Rousset) Implementation
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3' :                                     SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_thd_pnd_init : some initialization and namelist read
   !!   ice_thd_pnd      : main calling routine
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icetab         ! sea-ice: 1D <==> 2D transformation
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_pnd_init    ! routine called by icestp.F90
   PUBLIC   ice_thd_pnd         ! routine called by icestp.F90

   INTEGER ::              nice_pnd    ! choice of the type of pond scheme
   !                                   ! associated indices:
   INTEGER, PARAMETER ::   np_pndNO  = 0   ! No pond scheme
   INTEGER, PARAMETER ::   np_pndCST = 1   ! Constant ice pond scheme
   INTEGER, PARAMETER ::   np_pndLEV = 2   ! Level ice pond scheme

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_pnd
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_thd_pnd   ***
      !!               
      !! ** Purpose :   change melt pond fraction and thickness
      !!                
      !!-------------------------------------------------------------------
      !
      SELECT CASE ( nice_pnd )
      !
      CASE (np_pndCST)   ;   CALL pnd_CST    !==  Constant melt ponds  ==!
         !
      CASE (np_pndLEV)   ;   CALL pnd_LEV    !==  Level ice melt ponds  ==!
         !
      END SELECT
      !
   END SUBROUTINE ice_thd_pnd 


   SUBROUTINE pnd_CST 
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE pnd_CST  ***
      !!
      !! ** Purpose :   Compute melt pond evolution
      !!
      !! ** Method  :   Melt pond fraction and thickness are prescribed 
      !!                to non-zero values when t_su = 0C
      !!
      !! ** Tunable parameters : pond fraction (rn_apnd), pond depth (rn_hpnd)
      !!                
      !! ** Note   : Coupling with such melt ponds is only radiative
      !!             Advection, ridging, rafting... are bypassed
      !!
      !! ** References : Bush, G.W., and Trump, D.J. (2017)
      !!-------------------------------------------------------------------
      INTEGER  ::   ji        ! loop indices
      !!-------------------------------------------------------------------
      DO ji = 1, npti
         !
         IF( a_i_1d(ji) > 0._wp .AND. t_su_1d(ji) >= rt0 ) THEN
            h_ip_1d(ji)      = rn_hpnd    
            a_ip_1d(ji)      = rn_apnd * a_i_1d(ji)
            h_il_1d(ji)      = 0._wp    ! no pond lids whatsoever
         ELSE
            h_ip_1d(ji)      = 0._wp    
            a_ip_1d(ji)      = 0._wp
            h_il_1d(ji)      = 0._wp
         ENDIF
         !
      END DO
      !
   END SUBROUTINE pnd_CST


   SUBROUTINE pnd_LEV
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE pnd_LEV  ***
      !!
      !! ** Purpose : Compute melt pond evolution
      !!
      !! ** Method  : A fraction of meltwater is accumulated in ponds and sent to ocean when surface is freezing
      !!              We  work with volumes and then redistribute changes into thickness and concentration
      !!              assuming linear relationship between the two. 
      !!
      !! ** Action  : - pond growth:      Vp = Vp + dVmelt                                          --- from Holland et al 2012 ---
      !!                                     dVmelt = (1-r)/rhow * ( rhoi*dh_i + rhos*dh_s ) * a_i
      !!                                        dh_i  = meltwater from ice surface melt
      !!                                        dh_s  = meltwater from snow melt
      !!                                        (1-r) = fraction of melt water that is not flushed
      !!
      !!              - limtations:       a_ip must not exceed (1-r)*a_i
      !!                                  h_ip must not exceed 0.5*h_i
      !!
      !!              - pond shrinking:
      !!                       if lids:   Vp = Vp -dH * a_ip
      !!                                     dH = lid thickness change. Retrieved from this eq.:    --- from Flocco et al 2010 ---
      !!
      !!                                                                   rhoi * Lf * dH/dt = ki * MAX(Tp-Tsu,0) / H 
      !!                                                                      H = lid thickness
      !!                                                                      Lf = latent heat of fusion
      !!                                                                      Tp = -2C
      !!
      !!                                                                And solved implicitely as:
      !!                                                                   H(t+dt)**2 -H(t) * H(t+dt) -ki * (Tp-Tsu) * dt / (rhoi*Lf) = 0
      !!
      !!                    if no lids:   Vp = Vp * exp(0.01*MAX(Tp-Tsu,0)/Tp)                      --- from Holland et al 2012 ---
      !!
      !!              - Flushing:         w = -perm/visc * rho_oce * grav * Hp / Hi                 --- from Flocco et al 2007 ---
      !!                                     perm = permability of sea-ice
      !!                                     visc = water viscosity
      !!                                     Hp   = height of top of the pond above sea-level
      !!                                     Hi   = ice thickness thru which there is flushing
      !!
      !!              - Corrections:      remove melt ponds when lid thickness is 10 times the pond thickness
      !!
      !!              - pond thickness and area is retrieved from pond volume assuming a linear relationship between h_ip and a_ip:
      !!                                  a_ip/a_i = a_ip_frac = h_ip / zaspect
      !!
      !! ** Tunable parameters : rn_apnd_max, rn_apnd_min, rn_pnd_flush
      !! 
      !! ** Note       :   mostly stolen from CICE
      !!
      !! ** References :   Flocco and Feltham (JGR, 2007)
      !!                   Flocco et al       (JGR, 2010)
      !!                   Holland et al      (J. Clim, 2012)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nlay_i) ::   ztmp           ! temporary array
      !!
      REAL(wp), PARAMETER ::   zaspect =  0.8_wp      ! pond aspect ratio
      REAL(wp), PARAMETER ::   zTp     = -2._wp       ! reference temperature
      REAL(wp), PARAMETER ::   zvisc   =  1.79e-3_wp  ! water viscosity
      !!
      REAL(wp) ::   zfr_mlt, zdv_mlt                  ! fraction and volume of available meltwater retained for melt ponding
      REAL(wp) ::   zdv_frz, zdv_flush                ! Amount of melt pond that freezes, flushes
      REAL(wp) ::   zhp                               ! heigh of top of pond lid wrt ssh
      REAL(wp) ::   zv_ip_max                         ! max pond volume allowed
      REAL(wp) ::   zdT                               ! zTp-t_su
      REAL(wp) ::   zsbr                              ! Brine salinity
      REAL(wp) ::   zperm                             ! permeability of sea ice
      REAL(wp) ::   zfac, zdum                        ! temporary arrays
      REAL(wp) ::   z1_rhow, z1_aspect, z1_Tp         ! inverse
      !!
      INTEGER  ::   ji, jk                            ! loop indices
      !!-------------------------------------------------------------------
      z1_rhow   = 1._wp / rhow 
      z1_aspect = 1._wp / zaspect
      z1_Tp     = 1._wp / zTp 

      DO ji = 1, npti
         !                                                            !----------------------------------------------------!
         IF( h_i_1d(ji) < rn_himin .OR. a_i_1d(ji) < epsi10 ) THEN    ! Case ice thickness < rn_himin or tiny ice fraction !
            !                                                         !----------------------------------------------------!
            !--- Remove ponds on thin ice or tiny ice fractions
            a_ip_1d(ji)      = 0._wp
            h_ip_1d(ji)      = 0._wp
            h_il_1d(ji)      = 0._wp
            !                                                         !--------------------------------!
         ELSE                                                         ! Case ice thickness >= rn_himin !
            !                                                         !--------------------------------!
            v_ip_1d(ji) = h_ip_1d(ji) * a_ip_1d(ji)   ! retrieve volume from thickness
            v_il_1d(ji) = h_il_1d(ji) * a_ip_1d(ji)
            !
            !------------------!
            ! case ice melting !
            !------------------!
            !
            !--- available meltwater for melt ponding ---!
            zdum    = -( dh_i_sum(ji)*rhoi + dh_s_mlt(ji)*rhos ) * z1_rhow * a_i_1d(ji)
            zfr_mlt = rn_apnd_min + ( rn_apnd_max - rn_apnd_min ) * at_i_1d(ji) !  = ( 1 - r ) = fraction of melt water that is not flushed
            zdv_mlt = MAX( 0._wp, zfr_mlt * zdum ) ! max for roundoff errors? 
            !
            !--- overflow ---!
            ! If pond area exceeds zfr_mlt * a_i_1d(ji) then reduce the pond volume
            !    a_ip_max = zfr_mlt * a_i
            !    => from zaspect = h_ip / (a_ip / a_i), set v_ip_max as: 
            zv_ip_max = zfr_mlt**2 * a_i_1d(ji) * zaspect
            zdv_mlt = MAX( 0._wp, MIN( zdv_mlt, zv_ip_max - v_ip_1d(ji) ) )

            ! If pond depth exceeds half the ice thickness then reduce the pond volume
            !    h_ip_max = 0.5 * h_i
            !    => from zaspect = h_ip / (a_ip / a_i), set v_ip_max as: 
            zv_ip_max = z1_aspect * a_i_1d(ji) * 0.25 * h_i_1d(ji) * h_i_1d(ji)
            zdv_mlt = MAX( 0._wp, MIN( zdv_mlt, zv_ip_max - v_ip_1d(ji) ) )
            
            !--- Pond growing ---!
            v_ip_1d(ji) = v_ip_1d(ji) + zdv_mlt
            !
            !--- Lid melting ---!
            IF( ln_pnd_lids )   v_il_1d(ji) = MAX( 0._wp, v_il_1d(ji) - zdv_mlt ) ! must be bounded by 0
            !
            !--- mass flux ---!
            IF( zdv_mlt > 0._wp ) THEN
               zfac = zdv_mlt * rhow * r1_rdtice                        ! melt pond mass flux < 0 [kg.m-2.s-1]
               wfx_pnd_1d(ji) = wfx_pnd_1d(ji) - zfac
               !
               zdum = zfac / ( wfx_snw_sum_1d(ji) + wfx_sum_1d(ji) )    ! adjust ice/snow melting flux > 0 to balance melt pond flux
               wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) * (1._wp + zdum)
               wfx_sum_1d(ji)     = wfx_sum_1d(ji)     * (1._wp + zdum)
            ENDIF

            !-------------------!
            ! case ice freezing ! i.e. t_su_1d(ji) < (zTp+rt0)
            !-------------------!
            !
            zdT = MAX( zTp+rt0 - t_su_1d(ji), 0._wp )
            !
            !--- Pond contraction (due to refreezing) ---!
            IF( ln_pnd_lids ) THEN
               !
               !--- Lid growing and subsequent pond shrinking ---! 
               zdv_frz = 0.5_wp * MAX( 0._wp, -v_il_1d(ji) + & ! Flocco 2010 (eq. 5) solved implicitly as aH**2 + bH + c = 0
                  &                    SQRT( v_il_1d(ji)**2 + a_ip_1d(ji)**2 * 4._wp * rcnd_i * zdT * rdt_ice / (rLfus * rhow) ) ) ! max for roundoff errors
               
               ! Lid growing
               v_il_1d(ji) = MAX( 0._wp, v_il_1d(ji) + zdv_frz )
               
               ! Pond shrinking
               v_ip_1d(ji) = MAX( 0._wp, v_ip_1d(ji) - zdv_frz )

            ELSE
               ! Pond shrinking
               v_ip_1d(ji) = v_ip_1d(ji) * EXP( 0.01_wp * zdT * z1_Tp ) ! Holland 2012 (eq. 6)
            ENDIF
            !
            !--- Set new pond area and depth ---! assuming linear relation between h_ip and a_ip_frac
            ! v_ip     = h_ip * a_ip
            ! a_ip/a_i = a_ip_frac = h_ip / zaspect (cf Holland 2012, fitting SHEBA so that knowing v_ip we can distribute it to a_ip and h_ip)
            a_ip_1d(ji)      = MIN( a_i_1d(ji), SQRT( v_ip_1d(ji) * z1_aspect * a_i_1d(ji) ) ) ! make sure a_ip < a_i
            h_ip_1d(ji)      = zaspect * a_ip_1d(ji) / a_i_1d(ji)

            !---------------!            
            ! Pond flushing !
            !---------------!
            ! height of top of the pond above sea-level
            zhp = ( h_i_1d(ji) * ( rau0 - rhoi ) + h_ip_1d(ji) * ( rau0 - rhow * a_ip_1d(ji) / a_i_1d(ji) ) ) * r1_rau0
            
            ! Calculate the permeability of the ice (Assur 1958, see Flocco 2010)
            DO jk = 1, nlay_i
               zsbr = - 1.2_wp                                  &
                  &   - 21.8_wp    * ( t_i_1d(ji,jk) - rt0 )    &
                  &   - 0.919_wp   * ( t_i_1d(ji,jk) - rt0 )**2 &
                  &   - 0.0178_wp  * ( t_i_1d(ji,jk) - rt0 )**3
               ztmp(jk) = sz_i_1d(ji,jk) / zsbr
            END DO
            zperm = MAX( 0._wp, 3.e-08_wp * MINVAL(ztmp)**3 )
            
            ! Do the drainage using Darcy's law
            zdv_flush   = -zperm * rau0 * grav * zhp * rdt_ice / (zvisc * h_i_1d(ji)) * a_ip_1d(ji) * rn_pnd_flush ! tunable rn_pnd_flush from hunke et al. (2013)
            zdv_flush   = MAX( zdv_flush, -v_ip_1d(ji) )
            v_ip_1d(ji) = v_ip_1d(ji) + zdv_flush
            
            !--- Set new pond area and depth ---! assuming linear relation between h_ip and a_ip_frac
            a_ip_1d(ji)      = MIN( a_i_1d(ji), SQRT( v_ip_1d(ji) * z1_aspect * a_i_1d(ji) ) ) ! make sure a_ip < a_i
            h_ip_1d(ji)      = zaspect * a_ip_1d(ji) / a_i_1d(ji)

            !--- Corrections and lid thickness ---!
            IF( ln_pnd_lids ) THEN
               !--- retrieve lid thickness from volume ---!
               IF( a_ip_1d(ji) > epsi10 ) THEN   ;   h_il_1d(ji) = v_il_1d(ji) / a_ip_1d(ji)
               ELSE                              ;   h_il_1d(ji) = 0._wp
               ENDIF
               !--- remove ponds if lids are much larger than ponds ---!
               IF ( h_il_1d(ji) > h_ip_1d(ji) * 10._wp ) THEN
                  a_ip_1d(ji)      = 0._wp
                  h_ip_1d(ji)      = 0._wp
                  h_il_1d(ji)      = 0._wp
               ENDIF
            ENDIF
            !
         ENDIF
         
      END DO
      !
   END SUBROUTINE pnd_LEV


   SUBROUTINE ice_thd_pnd_init 
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_thd_pnd_init   ***
      !!
      !! ** Purpose : Physical constants and parameters linked to melt ponds
      !!              over sea ice
      !!
      !! ** Method  :  Read the namthd_pnd  namelist and check the melt pond  
      !!               parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_pnd  
      !!-------------------------------------------------------------------
      INTEGER  ::   ios, ioptio   ! Local integer
      !!
      NAMELIST/namthd_pnd/  ln_pnd, ln_pnd_LEV , rn_apnd_min, rn_apnd_max, rn_pnd_flush, &
         &                          ln_pnd_CST , rn_apnd, rn_hpnd,         &
         &                          ln_pnd_lids, ln_pnd_alb
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )              ! Namelist namthd_pnd  in reference namelist : Melt Ponds  
      READ  ( numnam_ice_ref, namthd_pnd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namthd_pnd  in reference namelist' )
      REWIND( numnam_ice_cfg )              ! Namelist namthd_pnd  in configuration namelist : Melt Ponds
      READ  ( numnam_ice_cfg, namthd_pnd, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namthd_pnd in configuration namelist' )
      IF(lwm) WRITE ( numoni, namthd_pnd )
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_pnd_init: ice parameters for melt ponds'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namicethd_pnd:'
         WRITE(numout,*) '      Melt ponds activated or not                                 ln_pnd       = ', ln_pnd
         WRITE(numout,*) '         Level ice melt pond scheme                               ln_pnd_LEV   = ', ln_pnd_LEV
         WRITE(numout,*) '            Minimum ice fraction that contributes to melt ponds   rn_apnd_min  = ', rn_apnd_min
         WRITE(numout,*) '            Maximum ice fraction that contributes to melt ponds   rn_apnd_max  = ', rn_apnd_max
         WRITE(numout,*) '            Pond flushing efficiency                              rn_pnd_flush = ', rn_pnd_flush
         WRITE(numout,*) '         Constant ice melt pond scheme                            ln_pnd_CST   = ', ln_pnd_CST
         WRITE(numout,*) '            Prescribed pond fraction                              rn_apnd      = ', rn_apnd
         WRITE(numout,*) '            Prescribed pond depth                                 rn_hpnd      = ', rn_hpnd
         WRITE(numout,*) '         Frozen lids on top of melt ponds                         ln_pnd_lids  = ', ln_pnd_lids
         WRITE(numout,*) '         Melt ponds affect albedo or not                          ln_pnd_alb   = ', ln_pnd_alb
      ENDIF
      !
      !                             !== set the choice of ice pond scheme ==!
      ioptio = 0
      IF( .NOT.ln_pnd ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndNO     ;   ENDIF
      IF( ln_pnd_CST  ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndCST    ;   ENDIF
      IF( ln_pnd_LEV  ) THEN   ;   ioptio = ioptio + 1   ;   nice_pnd = np_pndLEV    ;   ENDIF
      IF( ioptio /= 1 )   &
         & CALL ctl_stop( 'ice_thd_pnd_init: choose either none (ln_pnd=F) or only one pond scheme (ln_pnd_LEV or ln_pnd_CST)' )
      !
      SELECT CASE( nice_pnd )
      CASE( np_pndNO )         
         IF( ln_pnd_alb  ) THEN ; ln_pnd_alb  = .FALSE. ; CALL ctl_warn( 'ln_pnd_alb=false when no ponds' )  ; ENDIF
         IF( ln_pnd_lids ) THEN ; ln_pnd_lids = .FALSE. ; CALL ctl_warn( 'ln_pnd_lids=false when no ponds' ) ; ENDIF
      CASE( np_pndCST )         
         IF( ln_pnd_lids ) THEN ; ln_pnd_lids = .FALSE. ; CALL ctl_warn( 'ln_pnd_lids=false when constant ponds' ) ; ENDIF
      END SELECT
      !
   END SUBROUTINE ice_thd_pnd_init
   
#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif 

   !!======================================================================
END MODULE icethd_pnd 
