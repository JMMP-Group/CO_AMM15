MODULE trcini_age
   !!======================================================================
   !!                         ***  MODULE trcini_age  ***
   !! TOP :   initialisation of the AGE tracer
   !!======================================================================
   !! History :   2.0  !  2007-12  (G. Nurser, G. Madec, C. Ethe ) Original code
   !!----------------------------------------------------------------------
   !! trc_ini_age   : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trc
   USE trcnam_age
   USE trcsms_age

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_age   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_age.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_age
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_age  ***  
      !!
      !! ** Purpose :   initialization for AGE model
      !!
      !!----------------------------------------------------------------------
      INTEGER    ::  jn
      CHARACTER(len = 20)  ::  cltra
      !!----------------------------------------------------------------------
      !
      CALL trc_nam_age
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_age: passive tracer age'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)

      rryear  = 1._wp / ( nyear_len(1) * rday )    ! recip number of seconds in one year
      rryear  = 1._wp / rday

      IF( .NOT. ln_rsttr ) trn(:,:,:,jp_bgc+1:jptra) = 0.
      !
   END SUBROUTINE trc_ini_age

   !!======================================================================
END MODULE trcini_age
