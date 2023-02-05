module mod_julfort

    use iso_fortran_env

contains

    SUBROUTINE mydfree(NV,R0,RV0,E,K,RO,PO,QO,PHASEO,DELTAO)
        USE CONSTANTS
        IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)

    !C  ****  Maximum radial grid dimension.
    !    INTEGER*4, PARAMETER :: NDIM=25000
    !C  ****  Maximum number of terms in asymptotic series.
    !    INTEGER*4, PARAMETER :: MNT=50
    !C  ----  Speed of light (1/alpha).
    !    DOUBLE PRECISION, PARAMETER :: SL=137.035999139D0
    !C  ----  Bohr radius (cm).
    !    DOUBLE PRECISION, PARAMETER :: A0B=5.2917721067D-9
    !C  ----  Hartree energy (eV).
    !    DOUBLE PRECISION, PARAMETER :: HREV=27.21138602D0
    !C  ----  Electron rest energy (eV).
    !    DOUBLE PRECISION, PARAMETER :: REV=510.9989461D3
        
        CHARACTER FILEN*16,VFNAME*20
        PARAMETER (PI=3.1415926535897932D0)
        DIMENSION DR0(NDIM)  ! Output from SGRID.
    !C  ****  Potential.
        DIMENSION R0(NDIM),RV0(NDIM)
    !C  ****  Output radial functions.
        COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
    !C  ****  Coulomb wave function parameters.
        COMMON/OCOUL/WAVNUM,ETA,DELTA

        DOUBLE PRECISION, intent(out)::PHASEO, DELTAO
        DOUBLE PRECISION, intent(out)::RO(NDIM), PO(NDIM), QO(NDIM)
        DOUBLE PRECISION, intent(in)::E
        INTEGER*4, intent(in)::K

    !Alok
    !  ****  Input potential.
    !      VFNAME='dhfs079.tab'

    !      OPEN(3,FILE=VFNAME)
    !      DO I=1,NDIM
    ! 1      CONTINUE
    !        READ(3,*,ERR=1,END=10) R0(I),RV0(I)
    !        NV=I
    !      ENDDO
    ! 10   CONTINUE
    !      CLOSE(3)
    !ALOK


    !CALOK      WRITE(6,'(/A,I5)') ' # Potential grid. Number of radii =',NV
        CALL SPLERR(R0,RV0,0.0D0,0.0D0,ERR,NV,1)
    !CALOK      WRITE(6,'(A,1P,E9.1)') ' # Spline interpolation error =',ERR
    !
    !  ****  Spline interpolation of the potential.
        CALL VINT(R0,RV0,NV)
        RANGE=VRANGE()

    !
    !  ****  High-energy limit of the Dirac inner phase shift.
    !
        CALL DELINF(HEDEL)
        WRITE(6,'(/,'' # Delta_infty = '',1P,E22.15)') HEDEL
        Z=RV0(NV)

    ! If mannual grid is to be provided
!ALOK
    !    FILEN='n'       !enter the filename of user grid

    !    IF(FILEN.NE.'N'.AND.FILEN.NE.'n') THEN
    !        IGRID=1
    !        OPEN(7,FILE=FILEN)
    !        DO I=1,NDIM
    !        READ(7,*) RAD(I)
    !        NGP=I
    !        ENDDO
    !        CLOSE(7)
    !    ELSE
    !        IGRID=0
    !    ENDIF

    ! If mannual grid is to be provided

    !      FILEN='y'   !'n'    !enter the filename of user grid

    !      IF(FILEN.NE.'N'.AND.FILEN.NE.'n') THEN
    !        IGRID=1
    !        OPEN(7,FILE=FILEN)
    !        DO I=1,NDIM   !3318
    !          READ(7,*,END=10) RAD(I)
    !          NGP=I
    !        ENDDO
    !10        CONTINUE
    !        CLOSE(7)
    !      ELSE
    !        IGRID=0
    !      ENDIF

! Calculating with manual grid
        FILEN='y'   !'n'    !enter the filename of user grid

        IF(FILEN.NE.'N'.AND.FILEN.NE.'n') THEN
          IGRID=1
          DO I=1,NV
            RAD(I)=R0(I)
            NGP=I
          ENDDO
        ELSE
          IGRID=0
        ENDIF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    !      K=2
    !      E=100
        EPS=1e-15
        IF(K.LT.0) THEN
            L=-K-1
        ELSE
            L=K
        ENDIF
        IF(IGRID.EQ.0) THEN
            NGP=2000
            WAVEL=2.0D0*PI/SQRT(E*(2.0D0+E/SL**2))
            DRN=WAVEL/40.0D0
            RN=DRN*DBLE(NGP-300)
    !        write(*,*) DRN, RN, SL
            CALL SGRID(RAD,DR0,RN,1.0D-6,DRN,NGP,NDIM,IERS)
            IF(IERS.NE.0) STOP 'Error in the grid definition (DF).'
        ENDIF
        CALL DFREE(E,EPS,PHASE,K,1)
    !      write(*,*) PHASE
!        PHASEO=PHASE
!        DELTAO=DELTA
        IF(IER.NE.0) THEN
            WRITE(6,'(A,I3)') 'Error in DFREE. IER =',IER
    !        GO TO 20
        ENDIF
    !C
    !      WRITE(6,1401) VFNAME,E,K,EPS,PHASE,DELTA,ETA
    !      WRITE(8,1401) VFNAME,E,K,EPS,PHASE,DELTA,ETA
    !1401   FORMAT(1X,1P,'# **** Dirac equation. Numerical potential.',
    !   1    /' #',6X,'Input potential file: ',A20,
    !   2    /' #',6X,'Free state: E=',E13.6,', K=',I4,'  (EPS=',E8.1,')'
    !   3    /' #',6X,'  Inner phase shift=',E22.15,
    !   4    /' #',6X,'Coulomb phase shift=',E22.15,'  (ETA=',E13.6,')')
    !      WRITE(6,1202) ILAST,RAD(ILAST)
        FILEN='dirac.dat'

        
    !C  ****  Radial wave functions printed in output file.
    !C
    !*     READ(5,'(A)') FILEN  ! Uncomment to change the output filename.
        OPEN(10,FILE=FILEN)
    ! 1501 FORMAT(1X,'# Radial wave functions calculated by RADIAL.')
    !      WRITE(10,1501)
    !      IF(IOPT.EQ.1) THEN
    !        WRITE(10,1101) VFNAME,N,L,EPS,E
    !      ELSE IF(IOPT.EQ.2) THEN
    !        WRITE(10,1201) VFNAME,E,L,EPS,PHASE,DELTA,ETA
    !        WRITE(10,1202) ILAST,RAD(ILAST)
    !      ELSE IF(IOPT.EQ.3) THEN
    !        WRITE(10,1301) VFNAME,N,K,EPS,E
    !      ELSE IF(IOPT.EQ.4) THEN
    !        WRITE(10,1401) VFNAME,E,K,EPS,PHASE,DELTA,ETA
    !        WRITE(10,1202) ILAST,RAD(ILAST)
    !      ENDIF
    !C
        NTAB=NGP
        DO I=NGP,1,-1
            IF(ABS(P(I)).GT.1.0D-35) THEN
            NTAB=I
            GO TO 30
            ENDIF
        ENDDO
    30   CONTINUE
    !C

    !Added by Alok
        WRITE(10,*) 'Inner phase shift'
        WRITE(10,*) PHASE
        WRITE(10,*) 'Coulomb phase shift'
        WRITE(10,*) DELTA

        WRITE(10,1502)
        1502  FORMAT(1X,'#',7X,'R',14X,'P(R)',12X,'Q(R)')
        DO I=1,NTAB
    !C  ----  Do not print values less than 1.0D-99  ------------------------
            IF(ABS(P(I)).LT.1.0D-98) P(I)=0.0D0
            IF(ABS(Q(I)).LT.1.0D-98) Q(I)=0.0D0
    !C  ---------------------------------------------------------------------
            WRITE(10,'(1X,1P,3E16.8)') RAD(I),P(I),Q(I)
        ENDDO
        CLOSE(10)

        RO=RAD
        PO=P
        QO=Q
        PHASEO=PHASE
        DELTAO=DELTA
        
    end SUBROUTINE mydfree


end module

!********************************************************************************
!********************************************************************************
module mod_julfort_wrapper

    use iso_c_binding

    use :: mod_julfort

contains

    subroutine mydfree_wrapper(NV,R0,RV0,energy,kappa,r,p,q,phase,delta) bind(C, name="mydfree")
        USE CONSTANTS
        implicit none
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
        !DIMENSION R0(NDIM),RV0(NDIM)
        integer :: NV, kappa
        real(c_double) :: R0(NDIM),RV0(NDIM), r(NDIM), p(NDIM), q(NDIM)
        real(c_double) :: energy
        real(c_double), intent(out):: phase, delta

        CALL mydfree(NV,R0,RV0,energy,kappa,r,p,q,phase,delta)

        write(*,*) phase, kappa

    end subroutine mydfree_wrapper

end module