MODULE par_age
   !!======================================================================
   !!                        ***  par_age  ***
   !! TOP :   set the AGE parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_age.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER, PUBLIC  :: jp_age     !: number of passive tracers in PISCES
   ! NB: We want to store both the the number of total tracer and just the age ones
   INTEGER, PUBLIC  :: jp_age_all !: number of passive tracers in PISCES

   !!======================================================================
END MODULE par_age
