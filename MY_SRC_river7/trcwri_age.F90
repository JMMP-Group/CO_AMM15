MODULE trcwri_age
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    age :   Output of age tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top &&  defined key_iomput
   !!----------------------------------------------------------------------
   !! trc_wri_age   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE par_age     
   USE trc         
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_age 

CONTAINS

   SUBROUTINE trc_wri_age
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn
      !!---------------------------------------------------------------------

      ! write the tracer concentrations in the file
      DO jn = 1, jp_age
         !dum_age(:,:,:) = trn(:,:,:,jn) 
         !where( dum_age(:,:,:) .GE. 10.**-9 ) dum_age(:,:,:) = trn(:,:,:,jn+jp_bgc) / trn(:,:,:,jn)
         cltra = TRIM( ctrcnm(jn+jp_bgc) )                  ! short title for tracer
         !CALL iom_put( cltra, dum_age )
         CALL iom_put( cltra, trn(:,:,:,jn+jp_bgc) )
      END DO

      !
   END SUBROUTINE trc_wri_age

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_age
CONTAINS
   SUBROUTINE trc_wri_age                     ! Empty routine  
   END SUBROUTINE trc_wri_age
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_age.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE trcwri_age
