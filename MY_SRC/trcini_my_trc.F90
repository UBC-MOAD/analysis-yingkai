MODULE trcini_my_trc
   !!======================================================================
   !!                         ***  MODULE trcini_my_trc  ***
   !! TOP :   initialisation of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ini_my_trc   : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE trcsms_my_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_my_trc   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_my_trc.F90 2787 2011-06-27 09:54:00Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_my_trc
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_my_trc  ***  
      !!
      !! ** Purpose :   initialization for MY_TRC model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------

      !                       ! Allocate MY_TRC arrays
      IF( trc_sms_my_trc_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_my_trc: unable to allocate MY_TRC arrays' )
      !ALLOCATE(con_dil(jpi, jpj))
      ! Assign structure to Ba file
      
      CALL fld_fill( sf_Ba_ini       , (/ sn_Ba_ini /)       , cn_dir, 'trc_ini_my_trc', 'Ini-condition', 'trc_ini'     )
      CALL fld_fill( sf_Ba_river     , (/ sn_Ba_river /)     , cn_dir, 'trc_ini_my_trc', 'read Ba data' , 'trc_source'  )
      CALL fld_fill( sf_Ba_boundary  , (/ sn_Ba_boundary /)  , cn_dir, 'trc_ini_my_trc', 'Read Boundary', 'trc_boundary')

      CALL fld_fill( sf_d18O_ini     , (/ sn_d18O_ini /)     , cn_dir, 'trc_ini_my_trc', 'Ini-condition', 'trc_ini'     )
      CALL fld_fill( sf_d18O_river   , (/ sn_d18O_river /)   , cn_dir, 'trc_ini_my_trc', 'read Ba data' , 'trc_source'  )
      CALL fld_fill( sf_d18O_boundary, (/ sn_d18O_boundary /), cn_dir, 'trc_ini_my_trc', 'Read Boundary', 'trc_boundary')

      CALL fld_fill( sf_ice          , (/ sn_ice /)          , cn_dir, 'trc_ini_my_trc', 'Read Ice'     , 'trc_ice'     )
      CALL fld_fill( sf_rnf          , (/ sn_rnf /)          , cn_dir, 'trc_ini_my_trc', 'Read RNF'     , 'trc_rnf'     )
      !
   END SUBROUTINE trc_ini_my_trc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_my_trc             ! Empty routine
   END SUBROUTINE trc_ini_my_trc
#endif

   !!======================================================================
END MODULE trcini_my_trc
