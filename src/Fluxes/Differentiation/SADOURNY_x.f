C     File: SADOURNY_x.f

C
C     DDX
C
      SUBROUTINE DDX(DF,F,DX,NX,NY)

      implicit none

      INTEGER :: NX, NY, nthreads
      INTEGER :: OMP_GET_NUM_THREADS
      REAL*8, DIMENSION(NX,  NY) :: F
      REAL*8, DIMENSION(NX-1,NY) :: DF
      REAL*8 :: DX

Cf2py intent(in) F
Cf2py intent(out) DF
Cf2py intent(in) DX
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      nthreads = OMP_GET_NUM_THREADS()
      call OMP_SET_NUM_THREADS(nthreads)
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(SHARED) PROVATE(I,J)

      DO I = 1,NX-1
        DO J = 1,NY
            DF(I,J) = (F(I+1,J) - F(I,J))/DX
        END DO
      END DO
      !$OPM END PARALLEL DO

      END

C
C     AVX
C
      SUBROUTINE AVX(AF,F,NX,NY)

      implicit none

      INTEGER :: NX, NY, nthreads
      INTEGER :: OMP_GET_NUM_THREADS
      REAL*8, DIMENSION(NX,  NY) :: F
      REAL*8, DIMENSION(NX-1,NY) :: AF

Cf2py intent(in) F
Cf2py intent(out) AF
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      nthreads = OMP_GET_NUM_THREADS()
      call OMP_SET_NUM_THREADS(nthreads)
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(SHARED) PROVATE(I,J)

      DO I = 1,NX-1
        DO J = 1,NY
            AF(I,J) = 0.5*(F(I+1,J) + F(I,J))
        END DO
      END DO
      !$OPM END PARALLEL DO

      END

C
C     DDX_PERIODIC
C
      SUBROUTINE DDX_PERIODIC(DF,F,DX,NX,NY)

      implicit none

      INTEGER :: NX, NY, nthreads
      INTEGER :: OMP_GET_NUM_THREADS
      REAL*8, DIMENSION(NX,  NY) :: F
      REAL*8, DIMENSION(NX+1,NY) :: DF
      REAL*8 :: DX

Cf2py intent(in) F
Cf2py intent(out) DF
Cf2py intent(in) DX
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      nthreads = OMP_GET_NUM_THREADS()
      call OMP_SET_NUM_THREADS(nthreads)
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(SHARED) PROVATE(I,J)

      DO J = 1,NY
        DF(1,J) = (F(1,J) - F(NX,J))/DX
        DO I = 2,NX
          DF(I,J) = (F(I,J) - F(I-1,J))/DX
        END DO
        DF(NX+1,J) = (F(1,J) - F(NX,J))/DX
      END DO
      !$OPM END PARALLEL DO

      END

C
C     AVX_PERIODIC
C
      SUBROUTINE AVX_PERIODIC(AF,F,NX,NY)

      implicit none

      INTEGER :: NX, NY, nthreads
      INTEGER :: OMP_GET_NUM_THREADS
      REAL*8, DIMENSION(NX,  NY) :: F
      REAL*8, DIMENSION(NX+1,NY) :: AF

Cf2py intent(in) F
Cf2py intent(out) AF
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      nthreads = OMP_GET_NUM_THREADS()
      call OMP_SET_NUM_THREADS(nthreads)
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(SHARED) PROVATE(I,J)

      DO J = 1,NY
        AF(1,J) = 0.5*(F(1,J) + F(NX,J))
        DO I = 2,NX
          AF(I,J) = 0.5*(F(I,J) + F(I-1,J))
        END DO
        AF(NX+1,J) = 0.5*(F(1,J) + F(NX,J))
      END DO
      !$OPM END PARALLEL DO

      END


C
C     DDX_WALLS
C
      SUBROUTINE DDX_WALLS(DF,F,DX,NX,NY)

      implicit none

      INTEGER :: NX, NY, nthreads
      INTEGER :: OMP_GET_NUM_THREADS
      REAL*8, DIMENSION(NX,   NY) :: F
      REAL*8, DIMENSION(NX+1, NY) :: DF
      REAL*8 :: DX

Cf2py intent(in) F
Cf2py intent(out) DF
Cf2py intent(in) DX
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      nthreads = OMP_GET_NUM_THREADS()
      call OMP_SET_NUM_THREADS(nthreads)
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(SHARED) PROVATE(I,J)

      DO J = 1,NY
        DF(1,J) = (F(1,J) + F(1,J))/DX
        DO I = 2,NX
          DF(I,J) = (F(I,J) - F(I-1,J))/DX
        END DO
        DF(I,NY+1) = -(F(NX,J) + F(NX,J))/DX
      END DO
      !$OPM END PARALLEL DO

      END


C
C     AVX_WALLS
C
      SUBROUTINE AVX_WALLS(AF,F,NX,NY)

      implicit none

      INTEGER :: NX, NY, nthreads
      INTEGER :: OMP_GET_NUM_THREADS
      REAL*8, DIMENSION(NX,   NY) :: F
      REAL*8, DIMENSION(NX+1, NY) :: AF

Cf2py intent(in) F
Cf2py intent(out) AF
Cf2py intent(in) NX
Cf2py intent(in) NY

C     LOOP INDEX
      INTEGER :: I,J

      nthreads = OMP_GET_NUM_THREADS()
      call OMP_SET_NUM_THREADS(nthreads)
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(SHARED) PROVATE(I,J)

      DO J = 1,NY
        AF(1,J) = 0.0
        DO I = 2,NX
          AF(I,J) = 0.5*(F(I,J) + F(I-1,J))
        END DO
        AF(NX+1,J) = 0.0
      END DO
      !$OPM END PARALLEL DO

      END

