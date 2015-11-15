include '../Differentiation/FV_x_WENO.f'
include '../Differentiation/FV_y_WENO.f'

C     File: FV_SW.f

      SUBROUTINE FV_SW_FLUX(FLUX,H,HU,HV,NX,NY,NZ)

C     Given F, compute WENO 3-5 interpolations
C        along the first dimension. 

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8, Dimension(3,NX,NY,NZ) :: FLUX
      REAL*8, Dimension(NX,NY) :: H,HU,HV
Cf2py intent(inout) FLUX
Cf2py intent(in) H
Cf2py intent(in) HU
Cf2py intent(in) HV
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Interpolants
      REAL*8, Dimension(NX,NY) :: HMX,HPX,HMY,HPY
      REAL*8, Dimension(NX,NY) :: UMX,UPX,UMY,UPY
      REAL*8, Dimension(NX,NY) :: VMX,VPX,VMY,VPY

C     Flux values
      REAL*8 :: FP,FM
      REAL*8 :: FXP,FXM,FYP,FYM

C     Winds
      integer :: EW,WW,NW,SW

C     Loop Indices
      integer :: I,J

C     Compute the interpolants
      HMX,HPX = AX(H)
      UMX,UPX = AX(UH)
      VMX,VPX = AX(VH)

      HMY,HPY = AY(H)
      UMY,UPY = AY(UH)
      VMY,VPY = AY(VH)

C     Loop through the points, computing the flux
      do I = 1,NX
          do J = 1,NY

C            Compute winds
             if (UH(I,J).GE.0.0) then
                 EW = 1
                 WW = 0
             else
                 EW = 0
                 WW = 1
             end if
                 
             if (VH(I,J).GE.0.0) then
                 NW = 1
                 SW = 0
             else
                 NW = 0
                 SW = 1
             end if
                 

C            U Flux
             if (EW.eq.1) then
                 FP = UPX(I,J)*UPX(I,J)/HPX(I,J)
     +                + 0.5*G*HMX(I,J)*HMX(I,J)
                 if (I.eq.1) then
                     FM = UPX(NX,J)*UPX(NX,J)/HPX(NX,J)
                 else
                     FM = UPX(I-1,J)*UPX(I-1,J)/HPX(I-1,J)
                 end if
             else
                 FM = UMX(I,J)*UMX(I,J)/HMX(I,J)
     +                + 0.5*G*HPX(I,J)*HPX(I,J)
                 if (I.eq.NX) then
                     FP = UPX(1,J)*UPX(1,J)/HPX(1,J)
                 else
                     FP = UPX(I+1,J)*UPX(I+1,J)/HPX(I+1,J)
                 end if
             end if

             


          end do
      end do

      end

