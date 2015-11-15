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
              TMM = F(NX-1,:)
              TM  = F(NX,:)
              T   = F(1,:)
              TP  = F(2,:)
              TPP = F(3,:) 
          else if (I.eq.2) then
              TMM = F(NX,:)
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
              TPP = F(1,:) 
          else
              TMM = F(I-2,:)
              TM  = F(I-1,:)
              T   = F(I,:)
              TP  = F(I+1,:)
              TPP = F(I+2,:) 
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

      SUBROUTINE DDX(DXM,DXP,F,DX,NX,NY)

C     Given F, compute WENO 3-5 derivatives
C        along the first dimension. 

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8, Dimension(NX,NY) :: DXP, DXM, F
      REAL*8, Dimension(2,1)   :: DX
Cf2py intent(in) F
Cf2py intent(out) DXP
Cf2py intent(out) DXM
Cf2py intent(in) DX
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Temporary arrays for storing slices
      REAL*8, Dimension(NY) :: TMMM, TMM, TM, T, TP, TPP, TPPP

C     Temporary arrays for storing smoothness values
      REAL*8, Dimension(NY) :: BLL, BL, BR, BRR

C     Temporary arrays for storing WENO interpolants
      REAL*8, Dimension(NY) :: PLL, PL, PR, PRR

C     Optimal coefficients
      REAL*8 :: G1, G2, G3, G4

C     Computed coefficients
      REAL*8, Dimension(NY) :: WPLL, WPL, WPR, WPRR, SP
      REAL*8 :: EPS

C     LOOP INDEX
      integer :: I

C     SOME CONSTANTS
      parameter (G1 = 0.0, G2 = -2.0/15.0)
      parameter (G3 = 19.0/15.0, G4 = -2.0/15.0)
      parameter (EPS = 1.0E-6)

      do I = 1,NX

C         Determine which slices to use
C            i.e. boundary conditions
          if (I.eq.1) then
              TMMM = F(NX-2,:)
              TMM  = F(NX-1,:)
              TM   = F(NX,:)
              T    = F(I,:)
              TP   = F(I+1,:)
              TPP  = F(I+2,:) 
              TPPP = F(I+3,:) 
          else if (I.eq.2) then
              TMMM = F(NX-1,:)
              TMM  = F(NX,:)
              TM   = F(I-1,:)
              T    = F(I,:)
              TP   = F(I+1,:)
              TPP  = F(I+2,:) 
              TPPP = F(I+3,:) 
          else if (I.eq.3) then
              TMMM = F(NX,:)
              TMM  = F(I-2,:)
              TM   = F(I-1,:)
              T    = F(I,:)
              TP   = F(I+1,:)
              TPP  = F(I+2,:) 
              TPPP = F(I+3,:) 
          else if (I.eq.NX) then
              TMMM = F(I-3,:)
              TMM  = F(I-2,:)
              TM   = F(I-1,:)
              T    = F(I,:)
              TP   = F(1,:)
              TPP  = F(2,:) 
              TPPP = F(3,:) 
          else if (I.eq.(NX-1)) then
              TMMM = F(I-3,:)
              TMM  = F(I-2,:)
              TM   = F(I-1,:)
              T    = F(I,:)
              TP   = F(I+1,:)
              TPP  = F(1,:) 
              TPPP = F(2,:) 
          else if (I.eq.(NX-2)) then
              TMMM = F(I-3,:)
              TMM  = F(I-2,:)
              TM   = F(I-1,:)
              T    = F(I,:)
              TP   = F(I+1,:)
              TPP  = F(I+2,:) 
              TPPP = F(1,:) 
          else
              TMMM = F(I-3,:)
              TMM  = F(I-2,:)
              TM   = F(I-1,:)
              T    = F(I,:)
              TP   = F(I+1,:)
              TPP  = F(I+2,:) 
              TPPP = F(I+3,:) 
          end if

C         Smoothness coefficients
          BLL = TMMM*(547.0*TMMM - 3882.0*TMM +  4642.0*TM - 1854.0*T)
     +        +  TMM*(             7043.0*TMM - 17246.0*TM + 7042.0*T)
     +        +   TM*(                          11003.0*TM - 9402.0*T)
     +        +    T*(                                       2107.0*T)

          BL  = TMM*(267.0*TMM - 1642.0*TM + 1602.0*T -  494.0*TP)
     +        +  TM*(            2843.0*TM - 5966.0*T + 1922.0*TP)
     +        +   T*(                        3443.0*T - 2522.0*TP)
     +        +  TP*(                                    547.0*TP)

          BR  =  TM*(547.0*TM - 2522.0*T + 1922.0*TP -  494.0*TPP)
     +        +   T*(           3443.0*T - 5966.0*TP + 1602.0*TPP)
     +        +  TP*(                      2843.0*TP - 1642.0*TPP)
     +        + TPP*(                                   267.0*TPP)

          BRR =    T*(2107.0*T -  9402.0*TP +  7042.0*TPP - 1854.0*TPPP)
     +        +   TP*(           11003.0*TP - 17246.0*TPP + 4642.0*TPPP)
     +        +  TPP*(                         7043.0*TPP - 3882.0*TPPP)
     +        + TPPP*(                                       547.0*TPPP)

C         Compute interpolations onto i+1/2
          PLL =  (1.0/12.0)*(-11.0*TMMM + 45.0*TMM - 69.0*TM + 35.0*T)
          PL  = -(1.0/12.0)*(-1.0*TMM + 3.0*TM + 9.0*T - 11.0*TP)
          PR  = -(1.0/12.0)*(-1.0*TM + 15.0*T - 15.0*TP + TPP)
          PRR = -(1.0/12.0)*(11.0*T - 9.0*TP - 3.0*TPP + TPPP)

          WPLL = G1/((EPS+BLL)**2)
          WPL  = G2/((EPS+BL)**2)
          WPR  = G3/((EPS+BR)**2)
          WPRR = G4/((EPS+BRR)**2)

          SP = WPLL + WPL + WPR + WPRR

          WPLL = WPLL/SP
          WPL  = WPL/SP
          WPR  = WPR/SP
          WPRR = WPRR/SP

          DXP(I,:) = (WPLL*PLL + WPL*PL + WPR*PR + WPRR*PRR)/DX(1,1)

C         Compute interpolations onto i+1/2
          PLL =  (1.0/12.0)*(TMMM - 3.0*TMM - 9.0*TM + 11.0*T)
          PL  =  (1.0/12.0)*(TMM - 15.0*TM + 15.0*T - TP)
          PR  =  (1.0/12.0)*(-11.0*TM + 9.0*T + 3.0*TP - TPP)
          PRR = -(1.0/12.0)*(35.0*T - 69.0*TP + 45.0*TPP - 11.0*TPPP)

          WPLL = G4/((EPS+BLL)**2)
          WPL  = G3/((EPS+BL)**2)
          WPR  = G2/((EPS+BR)**2)
          WPRR = G1/((EPS+BRR)**2)

          SP = WPLL + WPL + WPR + WPRR

          WPLL = WPLL/SP
          WPL  = WPL/SP
          WPR  = WPR/SP
          WPRR = WPRR/SP

          DXM(I,:) = (WPLL*PLL + WPL*PL + WPR*PR + WPRR*PRR)/DX(1,1)

      end do

      end
