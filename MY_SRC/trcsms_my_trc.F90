MODULE trcsms_my_trc
   !!======================================================================
   !!                         ***  MODULE trcsms_my_trc  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_sms_my_trc       : MY_TRC model main routine
   !! trc_sms_my_trc_alloc : allocate arrays specific to MY_TRC sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE fldread         ! time interpolation and file reading
   USE trdmod_oce
   USE trdmod_trc
   !USE ice_2

   IMPLICIT NONE
   PUBLIC

   PUBLIC   trc_sms_my_trc       ! called by trcsms.F90 module
   PUBLIC   trc_sms_my_trc_alloc ! called by trcini_my_trc.F90 module

   CHARACTER(len=100), PUBLIC :: cn_dir = './'                                                                                                  ! Root directory
   TYPE(FLD_N) :: sn_ice, sn_rnf, sn_Ba_boundary, sn_Ba_ini, sn_Ba_river, sn_d18O_boundary, sn_d18O_ini, sn_d18O_river                          ! file information array
   REAL(wp) , ALLOCATABLE, DIMENSION(:, :) :: precip, Ba, sornf, d18O                                                                            ! 2-D global array
   REAL(wp) , ALLOCATABLE, DIMENSION(:, :, :) :: Ba_ini, Ba_boundary, d18O_ini, d18O_boundary                                                   ! 3-D global array
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ice, sf_rnf, sf_Ba_ini, sf_Ba_boundary, sf_Ba_river, sf_d18O_ini, sf_d18O_boundary, sf_d18O_river ! file structure array (for allocation)

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_my_trc.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
      INTEGER, INTENT(in) ::   kt   ! ocean e-step index
      INTEGER :: i, j, k
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 ) CALL timing_start('trc_sms_my_trc')

      CALL fld_read (kt, 1, sf_Ba_river)
      CALL fld_read (kt, 1, sf_Ba_boundary)
      CALL fld_read (kt, 1, sf_Ba_ini)

      CALL fld_read (kt, 1, sf_d18O_river)
      CALL fld_read (kt, 1, sf_d18O_boundary)
      CALL fld_read (kt, 1, sf_d18O_ini)

      CALL fld_read (kt, 1, sf_ice)
      CALL fld_read (kt, 1, sf_rnf)
   
      Ba(:, :)             = sf_Ba_river(1)%fnow(:, :, 1)
      Ba_boundary(:, :, :) = sf_Ba_boundary(1)%fnow(:, :, :)
      Ba_ini(:, :, :)      = sf_Ba_ini(1)%fnow(:, :, :)

      d18O(:, :)             = sf_d18O_river(1)%fnow(:, :, 1)
      d18O_boundary(:, :, :) = sf_d18O_boundary(1)%fnow(:, :, :)
      d18O_ini(:, :, :)      = sf_d18O_ini(1)%fnow(:, :, :)

      precip(:, :)          = sf_ice(1)%fnow(:, :, 1)*1000.0
      sornf(:, :)           = sf_rnf(1)%fnow(:, :, 1)
      ! Initialization
      IF( kt < 5 ) THEN
         WRITE(*, *) '~~~~~~~~ Initialization ~~~~~~~~'
         trn(:, :, :, jpmyt1) = Ba_ini(:, :, :)
         trn(:, :, :, jpmyt2) = d18O_ini(:, :, :)
      ENDIF
      
      tra(:, :, 1, jpmyt1) = tra(:, :, 1, jpmyt1) + Ba*sornf/e3t(:, :, 1)/1000.0
      precip = precip * (1-qsr) ! <---- Impact of sea-ice cover on the received precipitation on the sea surface (qsr is sea-ice cover [0, 1])
      tra(:, :, 1, jpmyt2) = tra(:, :, 1, jpmyt2) + (d18O*sornf*0.9 + (-1.5)*(emps+precip+sornf) + (-18.0)*precip)/e3t(:, :, 1)/1000.0
      
      !WRITE(*, *) 'trn Lena:', trn(600, 450, 1, jpmyt1)
      !WRITE(*, *) 'trn diff:', trn(600, 450, 1, jpmyt1) - trb(600, 450, 1, jpmyt1)
      WRITE(*, *) 'Time step: ', kt


      IF(nn_timing == 1)  CALL timing_stop('trc_sms_my_trc')
   END SUBROUTINE trc_sms_my_trc


   INTEGER FUNCTION trc_sms_my_trc_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_my_trc_alloc  ***
      !!----------------------------------------------------------------------
      !
      INTEGER           ::   ierror ! error flag
      !ALLOCATE(con_dil(jpi, jpj))
      ALLOCATE(Ba(jpi,jpj), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_Ba_river(1), STAT=ierror)
      
      ALLOCATE(Ba_ini(jpi, jpj, jpk), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_Ba_ini(1), STAT=ierror)

      ALLOCATE(Ba_boundary(jpi, jpj, jpk), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_Ba_boundary(1), STAT=ierror)

      ALLOCATE(d18O(jpi,jpj), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_d18O_river(1), STAT=ierror)

      ALLOCATE(d18O_ini(jpi, jpj, jpk), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_d18O_ini(1), STAT=ierror)

      ALLOCATE(d18O_boundary(jpi, jpj, jpk), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_d18O_boundary(1), STAT=ierror)

      ALLOCATE(precip(jpi, jpj), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_ice(1), STAT=ierror)

      ALLOCATE(sornf(jpi, jpj), STAT=trc_sms_my_trc_alloc)
      ALLOCATE(sf_rnf(1), STAT=ierror)
    
      ALLOCATE(sf_Ba_river(1)%fnow(jpi, jpj, 1))
      ALLOCATE(sf_Ba_ini(1)%fnow(jpi, jpj, jpk))
      ALLOCATE(sf_Ba_boundary(1)%fnow(jpi, jpj, jpk))

      ALLOCATE(sf_d18O_river(1)%fnow(jpi, jpj, 1))
      ALLOCATE(sf_d18O_ini(1)%fnow(jpi, jpj, jpk))
      ALLOCATE(sf_d18O_boundary(1)%fnow(jpi, jpj, jpk))

      ALLOCATE(sf_ice(1)%fnow(jpi, jpj, 1))
      ALLOCATE(sf_rnf(1)%fnow(jpi, jpj, 1))
      !
      IF( trc_sms_my_trc_alloc /= 0 ) CALL ctl_warn('trc_sms_my_trc_alloc : failed to allocate arrays')
      !
   END FUNCTION trc_sms_my_trc_alloc


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_my_trc( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_my_trc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_my_trc
#endif

   !!======================================================================
END MODULE trcsms_my_trc
