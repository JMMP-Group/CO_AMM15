MODULE dom_oce
   !!======================================================================
   !!                       ***  MODULE dom_oce  ***
   !!       
   !! ** Purpose :   Define in memory all the ocean space domain variables
   !!======================================================================
   !! History :  1.0  ! 2005-10  (A. Beckmann, G. Madec)  reactivate s-coordinate 
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.5  ! 2012     (S. Mocavero, I. Epicoco) Add arrays associated
   !!                             to the optimization of BDY communications
   !!            3.7  ! 2015-11  (G. Madec) introduce surface and scale factor ratio
   !!             -   ! 2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   Agrif_Root    : dummy function used when lk_agrif=F
   !!   Agrif_CFixed  : dummy function used when lk_agrif=F
   !!   dom_oce_alloc : dynamical allocation of dom_oce arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters

   IMPLICIT NONE
   PUBLIC             ! allows the acces to par_oce when dom_oce is used (exception to coding rules)

   PUBLIC dom_oce_alloc  ! Called from nemogcm.F90

   !!----------------------------------------------------------------------
   !! time & space domain namelist
   !! ----------------------------
   !                                    !!* Namelist namdom : time & space domain *
   INTEGER , PUBLIC ::   nn_bathy        !: = 0/1 ,compute/read the bathymetry file
   REAL(wp), PUBLIC ::   rn_wd_ref_depth !: Reference depth for Wet and Dry
   REAL(wp), PUBLIC ::   rn_bathy        !: depth of flat bottom (active if nn_bathy=0; if =0 depth=jpkm1)
   REAL(wp), PUBLIC ::   rn_hmin         !: minimum ocean depth (>0) or minimum number of ocean levels (<0)
   REAL(wp), PUBLIC ::   rn_isfhmin      !: threshold to discriminate grounded ice to floating ice
   REAL(wp), PUBLIC ::   rn_e3zps_min    !: miminum thickness for partial steps (meters)
   REAL(wp), PUBLIC ::   rn_e3zps_rat    !: minimum thickness ration for partial steps
   INTEGER , PUBLIC ::   nn_msh          !: = 1 create a mesh-mask file
   REAL(wp), PUBLIC ::   rn_atfp         !: asselin time filter parameter
   REAL(wp), PUBLIC ::   rn_rdt          !: time step for the dynamics and tracer
   INTEGER , PUBLIC ::   nn_closea       !: =0 suppress closed sea/lake from the ORCA domain or not (=1)
   INTEGER , PUBLIC ::   nn_euler        !: =0 start with forward time step or not (=1)
   LOGICAL , PUBLIC ::   ln_iscpl       !: coupling with ice sheet
   LOGICAL , PUBLIC ::   ln_crs          !: Apply grid coarsening to dynamical model output or online passive tracers

   !! Free surface parameters
   !! =======================
   LOGICAL , PUBLIC :: ln_dynspg_exp     !: Explicit free surface flag
   LOGICAL , PUBLIC :: ln_dynspg_ts      !: Split-Explicit free surface flag

   !! Time splitting parameters
   !! =========================
   LOGICAL,  PUBLIC :: ln_bt_fw          !: Forward integration of barotropic sub-stepping
   LOGICAL,  PUBLIC :: ln_bt_av          !: Time averaging of barotropic variables
   LOGICAL,  PUBLIC :: ln_bt_auto        !: Set number of barotropic iterations automatically
   INTEGER,  PUBLIC :: nn_bt_flt         !: Filter choice
   INTEGER,  PUBLIC :: nn_baro           !: Number of barotropic iterations during one baroclinic step (rdt)
   REAL(wp), PUBLIC :: rn_bt_cmax        !: Maximum allowed courant number (used if ln_bt_auto=T)

   !! Horizontal grid parameters for domhgr
   !! =====================================
   INTEGER       ::   jphgr_msh          !: type of horizontal mesh
   !                                       !  = 0 curvilinear coordinate on the sphere read in coordinate.nc
   !                                       !  = 1 geographical mesh on the sphere with regular grid-spacing
   !                                       !  = 2 f-plane with regular grid-spacing
   !                                       !  = 3 beta-plane with regular grid-spacing
   !                                       !  = 4 Mercator grid with T/U point at the equator

   REAL(wp)      ::   ppglam0            !: longitude of first raw and column T-point (jphgr_msh = 1)
   REAL(wp)      ::   ppgphi0            !: latitude  of first raw and column T-point (jphgr_msh = 1)
   !                                                        !  used for Coriolis & Beta parameters (jphgr_msh = 2 or 3)
   REAL(wp)      ::   ppe1_deg           !: zonal      grid-spacing (degrees)
   REAL(wp)      ::   ppe2_deg           !: meridional grid-spacing (degrees)
   REAL(wp)      ::   ppe1_m             !: zonal      grid-spacing (degrees)
   REAL(wp)      ::   ppe2_m             !: meridional grid-spacing (degrees)

   !! Vertical grid parameter for domzgr
   !! ==================================
   REAL(wp)      ::   ppsur              !: ORCA r4, r2 and r05 coefficients
   REAL(wp)      ::   ppa0               !: (default coefficients)
   REAL(wp)      ::   ppa1               !:
   REAL(wp)      ::   ppkth              !:
   REAL(wp)      ::   ppacr              !:
   !
   !  If both ppa0 ppa1 and ppsur are specified to 0, then
   !  they are computed from ppdzmin, pphmax , ppkth, ppacr in dom_zgr
   REAL(wp)      ::   ppdzmin            !: Minimum vertical spacing
   REAL(wp)      ::   pphmax             !: Maximum depth
   !
   LOGICAL       ::   ldbletanh          !: Use/do not use double tanf function for vertical coordinates
   REAL(wp)      ::   ppa2               !: Double tanh function parameters
   REAL(wp)      ::   ppkth2             !:
   REAL(wp)      ::   ppacr2             !:

   !                                    !! old non-DOCTOR names still used in the model
   INTEGER , PUBLIC ::   ntopo           !: = 0/1 ,compute/read the bathymetry file
   REAL(wp), PUBLIC ::   e3zps_min       !: miminum thickness for partial steps (meters)
   REAL(wp), PUBLIC ::   e3zps_rat       !: minimum thickness ration for partial steps
   INTEGER , PUBLIC ::   nmsh            !: = 1 create a mesh-mask file
   REAL(wp), PUBLIC ::   atfp            !: asselin time filter parameter
   REAL(wp), PUBLIC ::   rdt             !: time step for the dynamics and tracer

   !                                                  !!! associated variables
   INTEGER , PUBLIC                 ::   neuler        !: restart euler forward option (0=Euler)
   REAL(wp), PUBLIC                 ::   atfp1         !: asselin time filter coeff. (atfp1= 1-2*atfp)
   REAL(wp), PUBLIC                 ::   r2dt          !: = 2*rdt except at nit000 (=rdt) if neuler=0

   !!----------------------------------------------------------------------
   !! space domain parameters
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC ::   lzoom      =  .FALSE.   !: zoom flag
   LOGICAL, PUBLIC ::   lzoom_e    =  .FALSE.   !: East  zoom type flag
   LOGICAL, PUBLIC ::   lzoom_w    =  .FALSE.   !: West  zoom type flag
   LOGICAL, PUBLIC ::   lzoom_s    =  .FALSE.   !: South zoom type flag
   LOGICAL, PUBLIC ::   lzoom_n    =  .FALSE.   !: North zoom type flag

   !                                     !!! domain parameters linked to mpp
   INTEGER, PUBLIC ::   nperio            !: type of lateral boundary condition
   INTEGER, PUBLIC ::   nimpp, njmpp      !: i- & j-indexes for mpp-subdomain left bottom
   INTEGER, PUBLIC ::   nreci, nrecj      !: overlap region in i and j
   INTEGER, PUBLIC ::   nproc             !: number for local processor
   INTEGER, PUBLIC ::   narea             !: number for local area
   INTEGER, PUBLIC ::   nbondi, nbondj    !: mark of i- and j-direction local boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondi_bdy(:)    !: mark i-direction local boundaries for BDY open boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondj_bdy(:)    !: mark j-direction local boundaries for BDY open boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondi_bdy_b(:)  !: mark i-direction of neighbours local boundaries for BDY open boundaries  
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondj_bdy_b(:)  !: mark j-direction of neighbours local boundaries for BDY open boundaries  

   INTEGER, PUBLIC ::   npolj             !: north fold mark (0, 3 or 4)
   INTEGER, PUBLIC ::   nlci, nldi, nlei  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER, PUBLIC ::   nlcj, nldj, nlej  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER, PUBLIC ::   noea, nowe        !: index of the local neighboring processors in
   INTEGER, PUBLIC ::   noso, nono        !: east, west, south and north directions
   INTEGER, PUBLIC ::   npne, npnw        !: index of north east and north west processor
   INTEGER, PUBLIC ::   npse, npsw        !: index of south east and south west processor
   INTEGER, PUBLIC ::   nbne, nbnw        !: logical of north east & north west processor
   INTEGER, PUBLIC ::   nbse, nbsw        !: logical of south east & south west processor
   INTEGER, PUBLIC ::   nidom             !: ???

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mig        !: local  ==> global domain i-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mjg        !: local  ==> global domain j-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mi0, mi1   !: global ==> local  domain i-index    !!bug ==> other solution?
   !                                                  ! (mi0=1 and mi1=0 if the global index is not in the local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mj0, mj1   !: global ==> local  domain j-index     !!bug ==> other solution?
   !                                                  ! (mi0=1 and mi1=0 if the global index is not in the local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nimppt, njmppt   !: i-, j-indexes for each processor
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ibonit, ibonjt   !: i-, j- processor neighbour existence
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nlcit , nlcjt    !: dimensions of every subdomain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nldit , nldjt    !: first, last indoor index for each i-domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nleit , nlejt    !: first, last indoor index for each j-domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: nfiimpp, nfipproc, nfilcit

   !!----------------------------------------------------------------------
   !! horizontal curvilinear coordinate and scale factors
   !! ---------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   glamt , glamu, glamv , glamf    !: longitude at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   gphit , gphiu, gphiv , gphif    !: latitude  at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1t   , e2t  , r1_e1t, r1_e2t   !: t-point horizontal scale factors    [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1u   , e2u  , r1_e1u, r1_e2u   !: horizontal scale factors at u-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1v   , e2v  , r1_e1v, r1_e2v   !: horizontal scale factors at v-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1f   , e2f  , r1_e1f, r1_e2f   !: horizontal scale factors at f-point [m]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e1e2t , r1_e1e2t                !: associated metrics at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e1e2u , r1_e1e2u , e2_e1u       !: associated metrics at u-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e1e2v , r1_e1e2v , e1_e2v       !: associated metrics at v-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e1e2f , r1_e1e2f                !: associated metrics at f-point
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ff_f, ff_t                      !: coriolis factor                   [1/s]

   !!----------------------------------------------------------------------
   !! vertical coordinate and scale factors
   !! ---------------------------------------------------------------------
   !                                !!* Namelist namzgr : vertical coordinate *
   LOGICAL, PUBLIC ::   ln_zco       !: z-coordinate - full step
   LOGICAL, PUBLIC ::   ln_zps       !: z-coordinate - partial step
   LOGICAL, PUBLIC ::   ln_sco       !: s-coordinate or hybrid z-s coordinate
   LOGICAL, PUBLIC ::   ln_mes       !: Multi Enveloped s-coordinate (MEs)
   LOGICAL, PUBLIC ::   ln_isfcav    !: presence of ISF 
   LOGICAL, PUBLIC ::   ln_linssh    !: variable grid flag

   !                                                        !  ref.   ! before  !   now   ! after  !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3t_0 ,   e3t_b ,   e3t_n ,  e3t_a   !: t- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3u_0 ,   e3u_b ,   e3u_n ,  e3u_a   !: u- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3v_0 ,   e3v_b ,   e3v_n ,  e3v_a   !: v- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3f_0           ,   e3f_n            !: f- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3w_0 ,   e3w_b ,   e3w_n            !: w- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::    e3uw_0 ,  e3uw_b ,  e3uw_n            !: uw-vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::    e3vw_0 ,  e3vw_b ,  e3vw_n            !: vw-vert. scale factor [m]

   !                                                        !  ref.   ! before  !   now   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gdept_0 , gdept_b , gdept_n   !: t- depth              [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gdepw_0 , gdepw_b , gdepw_n   !: w- depth              [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gde3w_0           , gde3w_n   !: w- depth (sum of e3w) [m]
   
   !                                                      !  ref. ! before  !   now   !  after  !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ht_0            ,    ht_n             !: t-depth              [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hu_0  ,    hu_b ,    hu_n ,    hu_a   !: u-depth              [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hv_0  ,    hv_b ,    hv_n ,    hv_a   !: u-depth              [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::           r1_hu_b , r1_hu_n , r1_hu_a   !: inverse of u-depth [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::           r1_hv_b , r1_hv_n , r1_hv_a   !: inverse of v-depth [1/m]


   INTEGER, PUBLIC ::   nla10              !: deepest    W level Above  ~10m (nlb10 - 1)
   INTEGER, PUBLIC ::   nlb10              !: shallowest W level Bellow ~10m (nla10 + 1) 

   !! 1D reference  vertical coordinate
   !! =-----------------====------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   gdept_1d, gdepw_1d !: reference depth of t- and w-points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   e3t_1d  , e3w_1d   !: reference vertical scale factors at T- and W-pts (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e3tp    , e3wp     !: ocean bottom level thickness at T and W points

!!gm  This should be removed from here....  ==>>> only used in domzgr at initialization phase
   !! s-coordinate and hybrid z-s-coordinate
   !! =----------------======---------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   gsigt, gsigw       !: model level depth coefficient at t-, w-levels (analytic)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   gsi3w              !: model level depth coefficient at w-level (sum of gsigw)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   esigt, esigw       !: vertical scale factor coef. at t-, w-levels

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hbatv , hbatf      !: ocean depth at the vertical of  v--f
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hbatt , hbatu      !:                                 t--u points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   scosrf, scobot     !: ocean surface and bottom topographies 
   !                                                                           !  (if deviating from coordinate surfaces in HYBRID)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hifv  , hiff       !: interface depth between stretching at v--f
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hift  , hifu       !: and quasi-uniform spacing             t--u points (m)
!!gm end

   !!----------------------------------------------------------------------
   !! masks, bathymetry
   !! ---------------------------------------------------------------------
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbathy             !: number of ocean level (=0, 1, ... , jpk-1)
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbkt               !: vertical index of the bottom last T- ocean level
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbku, mbkv         !: vertical index of the bottom last U- and W- ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   bathy              !: ocean depth (meters)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tmask_i            !: interior domain T-point mask
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tmask_h            !: internal domain T-point mask (Figure 8.5 NEMO book)

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   misfdep                 !: top first ocean level                (ISF)
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mikt, miku, mikv, mikf  !: first wet T-, U-, V-, F- ocean level (ISF)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   risfdep                 !: Iceshelf draft                       (ISF)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssmask, ssumask, ssvmask, ssfmask    !: surface mask at T-,U-, V- and F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:), TARGET :: tmask, umask, vmask, fmask   !: land/ocean mask at T-, U-, V- and F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:), TARGET :: wmask, wumask, wvmask        !: land/ocean mask at WT-, WU- and WV-pts

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   tpol, fpol          !: north fold mask (jperio= 3 or 4)

   !!----------------------------------------------------------------------
   !! calendar variables
   !! ---------------------------------------------------------------------
   INTEGER , PUBLIC ::   nyear         !: current year
   INTEGER , PUBLIC ::   nmonth        !: current month
   INTEGER , PUBLIC ::   nday          !: current day of the month
   INTEGER , PUBLIC ::   nhour         !: current hour
   INTEGER , PUBLIC ::   nminute       !: current minute
   INTEGER , PUBLIC ::   ndastp        !: time step date in yyyymmdd format
   INTEGER , PUBLIC ::   nday_year     !: current day counted from jan 1st of the current year
   INTEGER , PUBLIC ::   nsec_year     !: current time step counted in second since 00h jan 1st of the current year
   INTEGER , PUBLIC ::   nsec_month    !: current time step counted in second since 00h 1st day of the current month
   INTEGER , PUBLIC ::   nsec_week     !: current time step counted in second since 00h of last monday
   INTEGER , PUBLIC ::   nsec_day      !: current time step counted in second since 00h of the current day
   REAL(wp), PUBLIC ::   fjulday       !: current julian day 
   REAL(wp), PUBLIC ::   fjulstartyear !: first day of the current year in julian days
   REAL(wp), PUBLIC ::   adatrj        !: number of elapsed days since the begining of the whole simulation
   !                                   !: (cumulative duration of previous runs that may have used different time-step size)
   INTEGER , PUBLIC, DIMENSION(0: 2) ::   nyear_len     !: length in days of the previous/current/next year
   INTEGER , PUBLIC, DIMENSION(0:13) ::   nmonth_len    !: length in days of the months of the current year
   INTEGER , PUBLIC, DIMENSION(0:13) ::   nmonth_half   !: second since Jan 1st 0h of the current year and the half of the months
   INTEGER , PUBLIC, DIMENSION(0:13) ::   nmonth_end    !: second since Jan 1st 0h of the current year and the end of the months
   INTEGER , PUBLIC                  ::   nsec1jan000   !: second since Jan 1st 0h of nit000 year and Jan 1st 0h the current year

   !!----------------------------------------------------------------------
   !! mpp reproducibility
   !!----------------------------------------------------------------------



   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp_rep = .FALSE.   !: agrif flag


   !!----------------------------------------------------------------------
   !! agrif domain
   !!----------------------------------------------------------------------



   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .FALSE.   !: agrif flag


   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: dom_oce.F90 6140 2015-12-21 11:35:23Z timgraham $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   !!----------------------------------------------------------------------
   !! NOT 'key_agrif'      dummy function                     No AGRIF zoom
   !!----------------------------------------------------------------------
   LOGICAL FUNCTION Agrif_Root()
      Agrif_Root = .TRUE.
   END FUNCTION Agrif_Root

   CHARACTER(len=3) FUNCTION Agrif_CFixed()
      Agrif_CFixed = '0' 
   END FUNCTION Agrif_CFixed


   INTEGER FUNCTION dom_oce_alloc()
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(13) :: ierr
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( mig(jpi), mjg(jpj), nfiimpp(jpni,jpnj),  &
         &      nfipproc(jpni,jpnj), nfilcit(jpni,jpnj), STAT=ierr(1) )
         !
      ALLOCATE( nimppt(jpnij) , ibonit(jpnij) , nlcit(jpnij) , nlcjt(jpnij) ,     &
         &      njmppt(jpnij) , ibonjt(jpnij) , nldit(jpnij) , nldjt(jpnij) ,     &
         &                                      nleit(jpnij) , nlejt(jpnij) ,     &
         &      mi0(jpidta)   , mi1 (jpidta),  mj0(jpjdta)   , mj1 (jpjdta),      &
         &      tpol(jpiglo)  , fpol(jpiglo)                               , STAT=ierr(2) )
         !
      ALLOCATE( glamt(jpi,jpj) ,    glamu(jpi,jpj) ,  glamv(jpi,jpj) ,  glamf(jpi,jpj) ,     &
         &      gphit(jpi,jpj) ,    gphiu(jpi,jpj) ,  gphiv(jpi,jpj) ,  gphif(jpi,jpj) ,     &
         &       e1t (jpi,jpj) ,     e2t (jpi,jpj) , r1_e1t(jpi,jpj) , r1_e2t(jpi,jpj) ,     &
         &       e1u (jpi,jpj) ,     e2u (jpi,jpj) , r1_e1u(jpi,jpj) , r1_e2u(jpi,jpj) ,     &
         &       e1v (jpi,jpj) ,     e2v (jpi,jpj) , r1_e1v(jpi,jpj) , r1_e2v(jpi,jpj) ,     &
         &       e1f (jpi,jpj) ,     e2f (jpi,jpj) , r1_e1f(jpi,jpj) , r1_e2f(jpi,jpj) ,     &
         &      e1e2t(jpi,jpj) , r1_e1e2t(jpi,jpj)                                     ,     &
         &      e1e2u(jpi,jpj) , r1_e1e2u(jpi,jpj) , e2_e1u(jpi,jpj)                   ,     &
         &      e1e2v(jpi,jpj) , r1_e1e2v(jpi,jpj) , e1_e2v(jpi,jpj)                   ,     &
         &      e1e2f(jpi,jpj) , r1_e1e2f(jpi,jpj)                                     ,     &
         &       ff_f(jpi,jpj) ,     ff_t(jpi,jpj)                                     , STAT=ierr(3) )
         !
      ALLOCATE( gdept_0(jpi,jpj,jpk) , gdepw_0(jpi,jpj,jpk) , gde3w_0(jpi,jpj,jpk) ,     &
         &      gdept_b(jpi,jpj,jpk) , gdepw_b(jpi,jpj,jpk) ,                             &
         &      gdept_n(jpi,jpj,jpk) , gdepw_n(jpi,jpj,jpk) , gde3w_n(jpi,jpj,jpk) , STAT=ierr(4) )
         !
      ALLOCATE( e3t_0(jpi,jpj,jpk) , e3u_0(jpi,jpj,jpk) , e3v_0(jpi,jpj,jpk) , e3f_0(jpi,jpj,jpk) , e3w_0(jpi,jpj,jpk) ,   &
         &      e3t_b(jpi,jpj,jpk) , e3u_b(jpi,jpj,jpk) , e3v_b(jpi,jpj,jpk) ,                      e3w_b(jpi,jpj,jpk) ,   & 
         &      e3t_n(jpi,jpj,jpk) , e3u_n(jpi,jpj,jpk) , e3v_n(jpi,jpj,jpk) , e3f_n(jpi,jpj,jpk) , e3w_n(jpi,jpj,jpk) ,   & 
         &      e3t_a(jpi,jpj,jpk) , e3u_a(jpi,jpj,jpk) , e3v_a(jpi,jpj,jpk) ,                                             &
         !                                                          !
         &      e3uw_0(jpi,jpj,jpk) , e3vw_0(jpi,jpj,jpk) ,         &
         &      e3uw_b(jpi,jpj,jpk) , e3vw_b(jpi,jpj,jpk) ,         &               
         &      e3uw_n(jpi,jpj,jpk) , e3vw_n(jpi,jpj,jpk) ,     STAT=ierr(5) )                       
         !
      ALLOCATE( ht_0(jpi,jpj) , hu_0(jpi,jpj) , hv_0(jpi,jpj) ,                                           &
         &                      hu_b(jpi,jpj) , hv_b(jpi,jpj) , r1_hu_b(jpi,jpj) , r1_hv_b(jpi,jpj) ,     &
         &      ht_n(jpi,jpj) , hu_n(jpi,jpj) , hv_n(jpi,jpj) , r1_hu_n(jpi,jpj) , r1_hv_n(jpi,jpj) ,     &
         &                      hu_a(jpi,jpj) , hv_a(jpi,jpj) , r1_hu_a(jpi,jpj) , r1_hv_a(jpi,jpj) , STAT=ierr(6)  )
         !
         !
      ALLOCATE( gdept_1d(jpk) , gdepw_1d(jpk) ,                                     &
         &      e3t_1d  (jpk) , e3w_1d  (jpk) , e3tp (jpi,jpj), e3wp(jpi,jpj) ,     &
         &      gsigt   (jpk) , gsigw   (jpk) , gsi3w(jpk)    ,                     &
         &      esigt   (jpk) , esigw   (jpk)                                 , STAT=ierr(7) )
         !
      ALLOCATE( hbatv (jpi,jpj) , hbatf (jpi,jpj) ,     &
         &      hbatt (jpi,jpj) , hbatu (jpi,jpj) ,     &
         &      scosrf(jpi,jpj) , scobot(jpi,jpj) ,     &
         &      hifv  (jpi,jpj) , hiff  (jpi,jpj) ,     &
         &      hift  (jpi,jpj) , hifu  (jpi,jpj) , STAT=ierr(8) )

      ALLOCATE( mbathy(jpi,jpj) , bathy  (jpi,jpj) ,                                       &
         &     tmask_i(jpi,jpj) , tmask_h(jpi,jpj) ,                                       & 
         &     ssmask (jpi,jpj) , ssumask(jpi,jpj) , ssvmask(jpi,jpj) , ssfmask(jpi,jpj) , &
         &     mbkt   (jpi,jpj) , mbku   (jpi,jpj) , mbkv   (jpi,jpj) , STAT=ierr(9) )

! (ISF) Allocation of basic array   
      ALLOCATE( misfdep(jpi,jpj) , risfdep(jpi,jpj),     &
         &     mikt(jpi,jpj), miku(jpi,jpj), mikv(jpi,jpj) ,           &
         &     mikf(jpi,jpj), STAT=ierr(10) )

      ALLOCATE( tmask(jpi,jpj,jpk) , umask(jpi,jpj,jpk),     & 
         &      vmask(jpi,jpj,jpk) , fmask(jpi,jpj,jpk), STAT=ierr(11) )

      ALLOCATE( wmask(jpi,jpj,jpk) , wumask(jpi,jpj,jpk), wvmask(jpi,jpj,jpk) , STAT=ierr(12) )
      !
      dom_oce_alloc = MAXVAL(ierr)
      !
   END FUNCTION dom_oce_alloc

   !!======================================================================
END MODULE dom_oce
