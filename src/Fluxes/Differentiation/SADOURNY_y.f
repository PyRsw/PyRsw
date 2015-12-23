C     File: SADOURNY_x.f

C
C     DDY
C
      SUBROUTINE DDY(DF,F,DY,NX,NY)

      use omp_lib
      INTEGER :: NX, NY
      REAL*8, DIMENSION(NX, NY)   :: F
      REAL*8, DIMENSION(NX, NY-1) :: DF
      REAL*8 :: DY

Cf2py intent(in) F
Cf2py intent(out) DF
Cf2py intent(in) DY
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      call OMP_SET_NUM_THREADS(2)
      !$OMP PARALLEL DO
      DO I = 1,NX
        DO J = 1,NY-1
            DF(I,J) = (F(I,J+1) - F(I,J))/DY
        END DO
      END DO
      !$OPM END PARALLEL DO

      END

C
C     AVY
C
      SUBROUTINE AVY(AF,F,NX,NY)

      use omp_lib
      INTEGER :: NX, NY
      REAL*8, DIMENSION(NX, NY)   :: F
      REAL*8, DIMENSION(NX, NY-1) :: AF

Cf2py intent(in) F
Cf2py intent(out) AF
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      call OMP_SET_NUM_THREADS(2)
      !$OMP PARALLEL DO
      DO I = 1,NX
        DO J = 1,NY-1
            AF(I,J) = 0.5*(F(I,J+1) + F(I,J))
        END DO
      END DO
      !$OPM END PARALLEL DO

      END

C
C     DDY_PERIODIC
C
      SUBROUTINE DDY_PERIODIC(DF,F,DY,NX,NY)

      use omp_lib
      INTEGER :: NX, NY
      REAL*8, DIMENSION(NX, NY)   :: F
      REAL*8, DIMENSION(NX, NY+1) :: DF
      REAL*8 :: DY

Cf2py intent(in) F
Cf2py intent(out) DF
Cf2py intent(in) DY
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      call OMP_SET_NUM_THREADS(2)
      !$OMP PARALLEL DO
      DO I = 1,NX
        DF(I,1) = (F(I,1) - F(I,NY))/DY
        DO J = 2,NY
          DF(I,J) = (F(I,J) - F(I,J-1))/DY
        END DO
        DF(I,NY+1) = (F(I,1) - F(I,NY))/DY
      END DO
      !$OPM END PARALLEL DO

      END

C
C     AVY_PERIODIC
C
      SUBROUTINE AVY_PERIODIC(AF,F,NX,NY)

      use omp_lib
      INTEGER :: NX, NY
      REAL*8, DIMENSION(NX, NY)   :: F
      REAL*8, DIMENSION(NX, NY+1) :: AF

Cf2py intent(in) F
Cf2py intent(out) AF
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      call OMP_SET_NUM_THREADS(2)
      !$OMP PARALLEL DO
      DO I = 1,NX
        AF(I,1) = 0.5*(F(I,1) + F(I,NY))
        DO J = 2,NY
          AF(I,J) = 0.5*(F(I,J) + F(I,J-1))
        END DO
        AF(I,NY+1) = 0.5*(F(I,1) + F(I,NY))
      END DO
      !$OPM END PARALLEL DO

      END

C
C     DDY_WALLS
C
      SUBROUTINE DDY_WALLS(DF,F,DY,NX,NY)

      use omp_lib
      INTEGER :: NX, NY
      REAL*8, DIMENSION(NX, NY)   :: F
      REAL*8, DIMENSION(NX, NY+1) :: DF
      REAL*8 :: DY

Cf2py intent(in) F
Cf2py intent(out) DF
Cf2py intent(in) DY
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      call OMP_SET_NUM_THREADS(2)
      !$OMP PARALLEL DO
      DO I = 1,NX
        DF(I,1) = (F(I,1) + F(I,1))/DY
        DO J = 2,NY
          DF(I,J) = (F(I,J) - F(I,J-1))/DY
        END DO
        DF(I,NY+1) = -(F(I,NY) + F(I,NY))/DY
      END DO
      !$OPM END PARALLEL DO

      END


C
C     AVY_WALLS
C
      SUBROUTINE AVY_WALLS(AF,F,NX,NY)

      use omp_lib
      INTEGER :: NX, NY
      REAL*8, DIMENSION(NX, NY)   :: F
      REAL*8, DIMENSION(NX, NY+1) :: AF

Cf2py intent(in) F
Cf2py intent(out) AF
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      call OMP_SET_NUM_THREADS(2)
      !$OMP PARALLEL DO
      DO I = 1,NX
        AF(I,1) = 0.0
        DO J = 2,NY
          AF(I,J) = 0.5*(F(I,J) + F(I,J-1))
        END DO
        AF(I,NY+1) = 0.0
      END DO
      !$OPM END PARALLEL DO

      END
