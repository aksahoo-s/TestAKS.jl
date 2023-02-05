program mynew
    USE CONSTANTS
    IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
    DIMENSION R0(NDIM),RV0(NDIM)


    CHARACTER VFNAME*20

!    VFNAME1='dhfs079.tab'
    energy=200.
    kappa=80

    VFNAME='dhfs079.tab'

    OPEN(3,FILE=VFNAME)
    DO I=1,NDIM
1      CONTINUE
      READ(3,*,ERR=1,END=10) R0(I),RV0(I)
      NV=I
    ENDDO
10   CONTINUE
    CLOSE(3)

    CALL mydfree(NV,R0,RV0,energy,kappa,0)
    
end program mynew