C     File: FV_x_WENO.f

      SUBROUTINE AX(AXM,AXP,F,NX,NY)

C     Given F, compute WENO 3-5 interpolations
C        along the first dimension. 

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8, Dimension(NX,NY) :: AXM, AXP, F
Cf2py intent(in) F
Cf2py intent(out) AXM
Cf2py intent(out) AXP
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Temporary arrays for storing slices
      REAL*8, Dimension(NY) :: TMM, TM, T, TP, TPP

C     Temporary arrays for storing smoothness values
      REAL*8, Dimension(NY) :: BL, BC, BR

C     Temporary arrays for storing WENO interpolants
      REAL*8, Dimension(NY) :: PL, PC, PR

C     Optimal coefficients
      REAL*8 :: G1, G2, G3

C     Computed coefficients
      REAL*8, Dimension(NY) :: WPL, WPC, WPR, SP
      REAL*8 :: EPS

C     LOOP INDEX
      integer :: i

C     SOME CONSTANTS
      parameter (G1 = 0.1, G2 = 0.6, G3 = 0.3, EPS = 1.0E-6)

      do I = 1,NX

C         Determine which slices to use
C            i.e. boundary conditions
          if (I.eq.1) then
              TMM = 5*TM + 2*T
              TM  = 0.5*TP - 2.5*T
              T   = F(1,:)
              TP  = F(2,:)
              TPP = F(3,:) 
          else if (I.eq.2) then
              TMM = 5*TM + 2*T
              TM  = F(1,:)
              T   = F(2,:)
              TP  = F(3,:)
              TPP = F(4,:) 
          else if (I.eq.NX) then
              TMM = F(NX-2,:)
              TM  = F(NX-1,:)
              T   = F(NX,:)
              TP  = F(1,:)
              TPP = F(2,:) 
          else if (I.eq.(NX-1)) then
              TMM = F(NX-3,:)
              TM  = F(NX-2,:)
              T   = F(NX-1,:)
              TP  = F(NX,:)
              TPP = 5*TP + 2*T
          else
              TMM = F(I-2,:)
              TM  = F(I-1,:)
              T   = F(I,:)
              TP  = 0.5*TM - 2.5*T
              TPP = 5*TP + 2*T
          end if

C         Smoothness coefficients
          BL = (13.0/12.0)*((TMM - 2.0*TM +     T)**2)
     +         + (1.0/4.0)*((TMM - 4.0*TM + 3.0*T)**2)

          BC = (13.0/12.0)*((TM - 2.0*T + TP)**2) 
     +         + (1.0/4.0)*((TM         - TP)**2)

          BR = (13.0/12.0)*((    T - 2.0*TP + TPP)**2) 
     +         + (1.0/4.0)*((3.0*T - 4.0*TP + TPP)**2)

C         Compute interpolations onto i+1/2
          PL =  (2.0/6.0)*TMM - (7.0/6.0)*TM + (11.0/6.0)*T
          PC = -(1.0/6.0)*TM  + (5.0/6.0)*T  +  (2.0/6.0)*TP
          PR =  (2.0/6.0)*T   + (5.0/6.0)*TP -  (1.0/6.0)*TPP

          WPL = G1/((EPS+BL)**2)
          WPC = G2/((EPS+BC)**2)
          WPR = G3/((EPS+BR)**2)

          SP = WPL + WPC + WPR

          WPL = WPL/SP
          WPC = WPC/SP
          WPR = WPR/SP

          AXP(I,:) = WPL*PL + WPC*PC + WPR*PR

C         Compute interpolation onto i-1/2
          PL = -(1.0/6.0)*TMM + (5.0/6.0)*TM + (2.0/6.0)*T
          PC =  (2.0/6.0)*TM  + (5.0/6.0)*T  - (1.0/6.0)*TP
          PR = (11.0/6.0)*T   - (7.0/6.0)*TP + (2.0/6.0)*TPP

          WPL = G3/((EPS+BL)**2)
          WPC = G2/((EPS+BC)**2)
          WPR = G1/((EPS+BR)**2)

          SP = WPL + WPC + WPR

          WPL = WPL/SP
          WPC = WPC/SP
          WPR = WPR/SP

          AXM(I,:) = WPL*PL + WPC*PC + WPR*PR
      end do

      end



      SUBROUTINE DDX(DXOUT,F,DX,NX,NY)

C     Given F, compute WENO 3-5 derivatives
C        along the first dimension. 

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8, Dimension(NX,NY) :: DXOUT, F
      REAL*8, Dimension(2,1)   :: DX
Cf2py intent(in) F
Cf2py intent(out) DXOUT
Cf2py intent(in) DX
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Temporary arrays for storing slices
      REAL*8, Dimension(NY) :: TMM, TM, T, TP, TPP

C     Temporary arrays for storing smoothness values
      REAL*8, Dimension(NY) :: BL, BC, BR

C     Temporary arrays for storing WENO interpolants
      REAL*8, Dimension(NY) :: PL, PC, PR

C     Optimal coefficients
      REAL*8 :: G1, G2, G3

C     Computed coefficients
      REAL*8, Dimension(NY) :: WPL, WPC, WPR, SP
      REAL*8 :: EPS

C     LOOP INDEX
      integer :: I

C     SOME CONSTANTS
      parameter (G1 = 1.0/6.0, G2 = 4.0/6.0, G3 = 1.0/6.0, EPS = 1.0E-6)

      do I = 1,NX

C         Determine which slices to use
C            i.e. boundary conditions
          if (I.eq.1) then
              TMM = 5*TM + 2*T
              TM  = 0.5*TP - 2.5*T
              T   = F(1,:)
              TP  = F(2,:)
              TPP = F(3,:) 
          else if (I.eq.2) then
              TMM = 5*TM + 2*T
              TM  = F(1,:)
              T   = F(2,:)
              TP  = F(3,:)
              TPP = F(4,:) 
          else if (I.eq.NX) then
              TMM = F(NX-2,:)
              TM  = F(NX-1,:)
              T   = F(NX,:)
              TP  = F(1,:)
              TPP = F(2,:) 
          else if (I.eq.(NX-1)) then
              TMM = F(NX-3,:)
              TM  = F(NX-2,:)
              T   = F(NX-1,:)
              TP  = F(NX,:)
              TPP = 5*TP + 2*T
          else
              TMM = F(I-2,:)
              TM  = F(I-1,:)
              T   = F(I,:)
              TP  = 0.5*TM - 2.5*T
              TPP = 5*TP + 2*T
          end if

C         Smoothness coefficients
          BL = (13.0/12.0)*((TMM - 2.0*TM +     T)**2)
     +         + (1.0/4.0)*((TMM - 4.0*TM + 3.0*T)**2)

          BC = (13.0/12.0)*((TM - 2.0*T + TP)**2) 
     +         + (1.0/4.0)*((TM         - TP)**2)

          BR = (13.0/12.0)*((    T - 2.0*TP + TPP)**2) 
     +         + (1.0/4.0)*((3.0*T - 4.0*TP + TPP)**2)

C         Compute average over cell I
          PL =  0.5*TMM - 2.0*TM + 1.5*T
          PC = -0.5*TM           + 0.5*TP
          PR = -1.5*T   + 2.0*TP - 0.5*TPP

          WPL = G1/((EPS+BL)**2)
          WPC = G2/((EPS+BC)**2)
          WPR = G3/((EPS+BR)**2)

          SP = WPL + WPC + WPR

          WPL = WPL/SP
          WPC = WPC/SP
          WPR = WPR/SP

          DXOUT(I,:) = (WPL*PL + WPC*PC + WPR*PR)/DX(1,1)
C          DXOUT(I,:) = WPL*PL + WPC*PC + WPR*PR

      end do

      end
