C     File: Fortran_Code.f

      SUBROUTINE FLUX_SPLIT_X(FM,FP,H,UH,VH,G,NX,NY)

C     Split fluxes along x

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8 :: G
      REAL*8, Dimension(NX,NY) :: H,UH,VH
      REAL*8, Dimension(3,NX,NY) :: FM,FP
Cf2py intent(out) FM
Cf2py intent(out) FP
Cf2py intent(in) H
Cf2py intent(in) UH
Cf2py intent(in) VH
Cf2py intent(in) G
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Create variables
      REAL*8, Dimension(3,NX,NY) :: F
      REAL*8, Dimension(3,3) :: ALPHA

C     Some doubles
      REAL*8 :: AQ1, AQ2, AQ3, LAM1, LAM2, LAM3
      REAL*8 :: U, V, D, C, LB1, LB2, LB3, TMP 
      REAL*8 :: EPS

C     Eigenvalues
      REAL*8, Dimension(NX,NY) :: LAM1_ALL, LAM2_ALL, LAM3_ALL
      REAL*8, Dimension(NX,NY) :: LAM1_BAR, LAM2_BAR, LAM3_BAR

      integer :: Iu, Iv, Ih

C     LOOP INDEX
      integer :: I,J

C     REAL*8 :: print_tmp
C     print_tmp = 1.0
C     write(*,100) 'pt=', print_tmp
C 100 format(A,F)

C     Variable indices
      Iu = 1; Iv = 2; Ih = 3;
      EPS = 1E-6

C     First compute the eigenvalue arrays
      LB1 = 0.0; LB2 = 0.0; LB3 = 0.0
      do I = 1,NX
        do J = 1,NY

          U = UH(I,J)
          D = H(I,J)

C         The eigenvalues
          LAM1_ALL(I,J) = U/(D+EPS)
          LAM2_ALL(I,J) = U/(D+EPS) + SQRT(G*D)
          LAM3_ALL(I,J) = U/(D+EPS) - SQRT(G*D)

          LB1 = MAX(ABS(LAM1_ALL(I,J)),LB1)
          LB2 = MAX(ABS(LAM2_ALL(I,J)),LB2)
          LB3 = MAX(ABS(LAM3_ALL(I,J)),LB3)

        end do
      end do

C     print_tmp = 2.0
C     write(*,200) 'pt=', print_tmp
C 200 format(A,F)

C     Now compute the lambda bars
C     Using global LF
      do I = 1,NX
        do J = 1,NY

          LAM1_BAR(I,J) = LB1
          LAM2_BAR(I,J) = LB2
          LAM3_BAR(I,J) = LB3

        end do
      end do

C     print_tmp = 3.0
C     write(*,300) 'pt=', print_tmp
C 300 format(A,F)
      
C     First compute the eigenvalue arrays
      do I = 1,NX
        do J = 1,NY

          U = UH(I,J)
          V = VH(I,J)
          D = H(I,J)
          C = SQRT(G*D)

C         The eigenvalues
          LAM1 = LAM1_ALL(I,J)
          LAM2 = LAM2_ALL(I,J)
          LAM3 = LAM3_ALL(I,J)
          LB1 = LAM1_BAR(I,J)
          LB2 = LAM2_BAR(I,J)
          LB3 = LAM3_BAR(I,J)

C         print_tmp = 4.0
C         write(*,400) 'pt=', print_tmp
C 400     format(A,F)

C         Compute the alpha matrix
          TMP = 1.0/(LAM2-LAM3)
          ALPHA(Ih,Ih) = -(LB2*LAM3-LB3*LAM2)*TMP
          ALPHA(Ih,Iu) =  (LB2-LB3)*TMP
          ALPHA(Ih,Iv) = 0.0
          ALPHA(Iu,Ih) = -LAM2*LAM3*(LB2-LB3)*TMP
          ALPHA(Iu,Iu) =  (LB2*LAM2-LB3*LAM3)*TMP
          ALPHA(Iu,Iv) = 0.0
          TMP = V/((LAM2-LAM3)*(U-LAM2*D)*(U-LAM3*D))
          ALPHA(Iv,Ih) = -TMP*(
     +          (LAM2*(LB1-LB3)+LAM3*(LB2-LB1))*(U*U+D*LAM2*LAM3) 
     +         + U*(D*LAM2*LAM2*(LB3-LB1) + D*LAM3*LAM3*(LB1-LB2) 
     +              + LAM2*LAM3*(LB3-LB1) ) )
          ALPHA(Iv,Iu) =  TMP*(
     +          U*LB1*(LAM2-LAM3)*(1.0-D) + U*U*(LB2-LB3) 
     +        + D*(LB2*LAM3*(LAM3-U) + LB3*LAM2*(U-LAM3))
     +        + U*(LB3*LAM3 - LB2*LAM2))
          ALPHA(Iv,Iv) = LB1

C         print_tmp = 5.0
C         write(*,500) 'pt=', print_tmp
C 500     format(A,F)

C         Compute alpha*q
          AQ1 = ALPHA(Ih,Ih)*D + ALPHA(Ih,Iu)*U + ALPHA(Ih,Iv)*V 
          AQ2 = ALPHA(Iu,Ih)*D + ALPHA(Iu,Iu)*U + ALPHA(Iu,Iv)*V 
          AQ3 = ALPHA(Iv,Ih)*D + ALPHA(Iv,Iu)*U + ALPHA(Iv,Iv)*V 

C         print_tmp = 6.0
C         write(*,600) 'pt=', print_tmp
C 600     format(A,F)

C         Compute F, Fplus, and Fminus
          F(Ih,I,J) = U
          F(Iu,I,J) = U*U/(D+EPS) + 0.5*G*D*D
          F(Iv,I,J) = U*V/(D+EPS)

          FP(Ih,I,J) = 0.5*(F(Ih,I,J) + AQ1)
          FP(Iu,I,J) = 0.5*(F(Iu,I,J) + AQ2)
          FP(Iv,I,J) = 0.5*(F(Iv,I,J) + AQ3)

          FM(Ih,I,J) = 0.5*(F(Ih,I,J) - AQ1)
          FM(Iu,I,J) = 0.5*(F(Iu,I,J) - AQ2)
          FM(Iv,I,J) = 0.5*(F(Iv,I,J) - AQ3)

C         print_tmp = 7.0
C         write(*,700) 'pt=', print_tmp
C 700     format(A,F)

        end do
      end do

      end

      SUBROUTINE FLUX_SPLIT_Y(FM,FP,H,UH,VH,G,NX,NY)

C     Split fluxes along x

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8 :: G
      REAL*8, Dimension(NX,NY) :: H,UH,VH
      REAL*8, Dimension(3,NX,NY) :: FM,FP
Cf2py intent(out) FM
Cf2py intent(out) FP
Cf2py intent(in) H
Cf2py intent(in) UH
Cf2py intent(in) VH
Cf2py intent(in) G
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Create variables
      REAL*8, Dimension(3,NX,NY) :: F
      REAL*8, Dimension(3,3) :: ALPHA

C     Some doubles
      REAL*8 :: AQ1, AQ2, AQ3, LAM1, LAM2, LAM3
      REAL*8 :: U, V, D, C, LB1, LB2, LB3, TMP 
      REAL*8 :: EPS

C     Eigenvalues
      REAL*8, Dimension(NX,NY) :: LAM1_ALL, LAM2_ALL, LAM3_ALL
      REAL*8, Dimension(NX,NY) :: LAM1_BAR, LAM2_BAR, LAM3_BAR

      integer :: Iu, Iv, Ih

C     LOOP INDEX
      integer :: I,J

C     REAL*8 :: print_tmp
C     print_tmp = 1.0
C     write(*,100) 'pt=', print_tmp
C 100 format(A,F)

C     Variable indices
      Iu = 1; Iv = 2; Ih = 3;
      EPS = 1E-6

C     First compute the eigenvalue arrays
      LB1 = 0.0; LB2 = 0.0; LB3 = 0.0
      do I = 1,NX
        do J = 1,NY

          V = VH(I,J)
          D = H(I,J)

C         The eigenvalues
          LAM1_ALL(I,J) = V/(D+EPS)
          LAM2_ALL(I,J) = V/(D+EPS) + SQRT(G*D)
          LAM3_ALL(I,J) = V/(D+EPS) - SQRT(G*D)

          LB1 = MAX(ABS(LAM1_ALL(I,J)),LB1)
          LB2 = MAX(ABS(LAM2_ALL(I,J)),LB2)
          LB3 = MAX(ABS(LAM3_ALL(I,J)),LB3)

        end do
      end do

C     print_tmp = 2.0
C     write(*,200) 'pt=', print_tmp
C 200 format(A,F)

C     Now compute the lambda bars
C     Using global LF
      do I = 1,NX
        do J = 1,NY

          LAM1_BAR(I,J) = LB1
          LAM2_BAR(I,J) = LB2
          LAM3_BAR(I,J) = LB3

        end do
      end do

C     print_tmp = 3.0
C     write(*,300) 'pt=', print_tmp
C 300 format(A,F)
      
C     First compute the eigenvalue arrays
      do I = 1,NX
        do J = 1,NY

          U = UH(I,J)
          V = VH(I,J)
          D = H(I,J)
          C = SQRT(G*D)

C         The eigenvalues
          LAM1 = LAM1_ALL(I,J)
          LAM2 = LAM2_ALL(I,J)
          LAM3 = LAM3_ALL(I,J)
          LB1 = LAM1_BAR(I,J)
          LB2 = LAM2_BAR(I,J)
          LB3 = LAM3_BAR(I,J)

C         print_tmp = 4.0
C         write(*,400) 'pt=', print_tmp
C 400     format(A,F)

C         Compute the alpha matrix
          TMP = 1.0/(LAM2-LAM3)
          ALPHA(Ih,Ih) = -(LB2*LAM3-LB3*LAM2)*TMP
          ALPHA(Ih,Iv) =  (LB2-LB3)*TMP
          ALPHA(Ih,Iu) = 0.0
          ALPHA(Iv,Ih) = -LAM2*LAM3*(LB2-LB3)*TMP
          ALPHA(Iv,Iv) =  (LB2*LAM2-LB3*LAM3)*TMP
          ALPHA(Iv,Iu) = 0.0
          TMP = U/((LAM2-LAM3)*(V-LAM2*D)*(V-LAM3*D))
          ALPHA(Iu,Ih) = -TMP*(
     +          (LAM2*(LB1-LB3)+LAM3*(LB2-LB1))*(V*V+D*LAM2*LAM3) 
     +         + V*(D*LAM2*LAM2*(LB3-LB1) + D*LAM3*LAM3*(LB1-LB2) 
     +              + LAM2*LAM3*(LB3-LB1) ) )
          ALPHA(Iu,Iv) =  TMP*(
     +          V*LB1*(LAM2-LAM3)*(1.0-D) + V*V*(LB2-LB3) 
     +        + D*(LB2*LAM3*(LAM3-V) + LB3*LAM2*(V-LAM3))
     +        + V*(LB3*LAM3 - LB2*LAM2))
          ALPHA(Iu,Iu) = LB1

C         print_tmp = 5.0
C         write(*,500) 'pt=', print_tmp
C 500     format(A,F)

C         Compute alpha*q
          AQ1 = ALPHA(Ih,Ih)*D + ALPHA(Ih,Iu)*U + ALPHA(Ih,Iv)*V 
          AQ2 = ALPHA(Iu,Ih)*D + ALPHA(Iu,Iu)*U + ALPHA(Iu,Iv)*V 
          AQ3 = ALPHA(Iv,Ih)*D + ALPHA(Iv,Iu)*U + ALPHA(Iv,Iv)*V 

C         print_tmp = 6.0
C         write(*,600) 'pt=', print_tmp
C 600     format(A,F)

C         Compute F, Fplus, and Fminus
          F(Ih,I,J) = V
          F(Iu,I,J) = U*V/(D+EPS)
          F(Iv,I,J) = V*V/(D+EPS) + 0.5*G*D*D

          FP(Ih,I,J) = 0.5*(F(Ih,I,J) + AQ1)
          FP(Iu,I,J) = 0.5*(F(Iu,I,J) + AQ2)
          FP(Iv,I,J) = 0.5*(F(Iv,I,J) + AQ3)

          FM(Ih,I,J) = 0.5*(F(Ih,I,J) - AQ1)
          FM(Iu,I,J) = 0.5*(F(Iu,I,J) - AQ2)
          FM(Iv,I,J) = 0.5*(F(Iv,I,J) - AQ3)

C         print_tmp = 7.0
C         write(*,700) 'pt=', print_tmp
C 700     format(A,F)

        end do
      end do

      end

