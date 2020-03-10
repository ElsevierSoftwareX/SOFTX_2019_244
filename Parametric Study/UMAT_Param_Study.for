!   This file is part of MicroFEA 1.0.

!   MicroFEA 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   MicroFEA 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.

!   You should have received a copy of the GNU General Public License
!   along with MicroFEA 1.0 if not, see <https://www.gnu.org/licenses/>.
!
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,       &
     & DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,        &
     & DTEMP,PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS,    &
     & NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL,      &
     & NPT, LAYER, KSPT, KSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*8 CMNAME

      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),    &
     & DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),       &
     & PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),       &
     & DFGRD0(3, 3), DFGRD1(3, 3), TEMP(1), ETHERM(6)

!     ------------------------------------------------------------------
!                     Functionally Graded Material UMAT
!                Volume Fraction Variation based on B-Spline
!                 Written by Marcelo S. Medeiros 01/25/2018
!     ------------------------------------------------------------------

      INTEGER K1, K2, N, P, S, CNT
      DOUBLE PRECISION UK, RAD, B, PC, Z, CONST_MAT
      DOUBLE PRECISION THICKNESS, Vf, E_MAT, E_INC, NU_MAT, NU_INC
      DOUBLE PRECISION K_MAT, K_INC, ETHERM
      DIMENSION UK(13), PC(9), B(4), CONST_MAT(6, 6)
	  
!     ------------------------------------------------------------------
!                             INPUT PARAMETERS
!     ------------------------------------------------------------------
      THICKNESS = 1.0D-2               ! Plate Thickness
      
!      DO CNT=1, 13
!        UK(CNT) = PROPS(CNT)          ! Knot Vector
!      END DO

!      DO CNT=1,9
!        PC(CNT) = PROPS(CNT+13)       ! Control Point Vector
!      END DO
      
!      N = PROPS(23)
!      P = PROPS(24)

      E_MAT = 70.0D9				  ! Matrix Elastic Modulus (Pa)
      E_INC = 380.0D9				  ! Inclusion's Elastic Modulus (Pa)
      NU_MAT = 0.3D0			      ! Matrix Poisson's Ratio
      NU_INC = 0.3D0     			  ! Inclusion's Poisson's Ratio	  
      CTE_MAT = 15.D-6                ! CTE of the Matrix
	  CTE_INC = 10.D-6                ! CTE of the Inclusions
      Z = COORDS(3)                   ! 'Z' coordinate of the Gauss Point

      RAD = Z/THICKNESS               ! Normalized Coord.

!     ------------------------------------------------------------------
!		  BULK MODULI FOR COEFF. OF THERMAL EXPANSION CALCULATION

      K_MAT = E_MAT/(3.D0*(1.D0 - 2.D0*NU_MAT))
      K_INC = E_INC/(3.D0*(1.D0 - 2.D0*NU_INC))

!     ------------------------------------------------------------------
!           CALCULATE THE VOLUME FRACTION AT THE CURRENT POSITION    
      
!      CALL NURBS (UK, N, S, P, PC, B, RAD, Vf)

      Vf = (Z/THICKNESS)**2.0D0
!     ------------------------------------------------------------------
! 				 CALL SUBROUTINE FOR CTE HOMOGENIZATION

      CALL TURNER (K_MAT, K_INC,CTE_MAT, CTE_INC, Vf, CTE)
	  
!     ------------------------------------------------------------------
!                     CALL HOMOGENIZATION SUBROUTINE
!              Uncomment the Homogenization model to be used

      CALL MORI_TANAKA (E_MAT, E_INC, NU_MAT, NU_INC, Vf, CONST_MAT)
!      CALL VOIGT (E_MAT, E_INC, NU_MAT, NU_INC, Vf, CONST_MAT)


!     ------------------------------------------------------------------
!                            CALCULATE STRESS
      DO K1=1, NDI
        ETHERM(K1) = CTE*(TEMP(1) + DTEMP)
        ETHERM(K1+3) = 0.0D0
      END DO
	  
	  DO K1=1,NTENS
	    DO K2=1,NTENS
	      DDSDDE(K1,K2) = CONST_MAT(K1,K2)
		END DO
	  END DO
	  
      DO K1=1,NTENS
        DO K2=1,NTENS
          STRESS(K2)=STRESS(K2)+CONST_MAT(K2,K1)*((STRAN(K1) +          &
     &      DSTRAN(K1)) - ETHERM(K1))
        END DO
      END DO

      RETURN
      END SUBROUTINE
!     ------------------------------------------------------------------
!                CALCULATE Vf BASED ON B-SPLINE INPUT
!     ------------------------------------------------------------------

      SUBROUTINE NURBS (UK, N, S, P, PC, B, RAD, Vf)
      
      DOUBLE PRECISION UK, PC, B, RAD, Vf
      INTEGER N, S, P, CNT

      DIMENSION UK(N+P+1), B(P+1), PC(10)

      CALL FINDSPAN(N, P, RAD, UK, S)
      CALL BASISFUN(S, RAD, P, N, UK, B)
      Vf = 0.D0
      DO CNT = 1, P + 1
        Vf = Vf + (B(CNT)*PC(S-P+CNT))
      END DO

      RETURN
      END SUBROUTINE

!     ------------------------------------------------------------------
!                                FIND SPAN
!     ------------------------------------------------------------------

      SUBROUTINE FINDSPAN(N, P, U, UK, S)

      INTEGER S, N, P, LOW, MID, HIGH
      DOUBLE PRECISION UK, U
      DIMENSION UK(N+P+1)

      IF (U == UK(N+2)) THEN
        S = N
      END IF
      LOW = P
      HIGH = N+1
      MID = (LOW+HIGH)/2
      DO WHILE (U < UK(MID+1) .OR. U >= UK(MID+2))
        IF (U < UK(MID+1)) THEN
            HIGH = MID
        ELSE
            LOW = MID
        END IF
        MID = (LOW+HIGH)/2
      END DO
      S = MID

      RETURN
      END SUBROUTINE


!     ------------------------------------------------------------------
!                              BASIS FUNCTION
!     ------------------------------------------------------------------

      SUBROUTINE BASISFUN(IV, UV, P, N, UK, B)

      INTEGER P, R, J, K1, CNT, I, IV, N
      DOUBLE PRECISION UV, U, UK, B
      DOUBLE PRECISION LEFT, RIGHT, SAVED, TEMP

      DIMENSION B(P+1)
      DIMENSION LEFT(P+1), RIGHT(P+1), UK(N+P+1)

      DO K1 = 1, P+1
        B(K1) = 0.D0
        LEFT(K1) = 0.D0
        RIGHT(K1) = 0.D0
      END DO

      I = IV + 1
      U = UV
      B(1) = 1
      DO J=1, P
        LEFT(J+1) = U - UK(I+1-J)
        RIGHT(J+1) = UK(I+J) - U
        SAVED = 0
        DO R=0, J-1
          TEMP = B(R+1)/(RIGHT(R+2) + LEFT(J-R+1))
          B(R+1) = SAVED + RIGHT(R+2)*TEMP
          SAVED = LEFT(J-R+1)*TEMP
        END DO
        B(J+1) = SAVED
      END DO


      RETURN
      END SUBROUTINE
	  
      SUBROUTINE MORI_TANAKA (E_MAT,E_INC,NU_MAT,NU_INC,Vf,CONST_MAT)
!     ------------------------------------------------------------------
!            Implementation of Mori-Tanaka Homogenization Scheme
!         As per "A New Approach to the Application of Mori-Tanaka's
!             Theory in Composite Materials", Y. Benveniste (1987)
!                     DOI: 10.1016/0167-6636(87)90005-6
!            Written by Marcelo Medeiros on February 3rd, 2018
!
!     ------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION ALPHA, ESHELBY, COMPL, IDENT, Vm, Vf
      DOUBLE PRECISION E_MAT, NU_MAT, E_INC, NU_INC, C_MAT, C_INC, B
      DOUBLE PRECISION MAT_DIFF, AUX1, AUX2, CONST_MAT
      INTEGER K1, K2
      DIMENSION ESHELBY (6,6), COMPL(6,6), C_MAT(6,6), C_INC(6,6)
      DIMENSION MAT_DIFF(6,6), AUX1(6,6), AUX2(6,6), B(6,6), IDENT(6,6)
      DIMENSION CONST_MAT(6,6)

      DATA IDENT /36*0.0D0/, AUX1 /36*0.D0/, AUX2 /36*0.D0/
!      DATA CONST_MAT /36*0.0D0/

      ALPHA = 1.0D0               ! Inclusion's aspect ratio

      Vm = 1.0D0 - Vf             ! Volume Fraction of the Matrix
      
      DO K1=1, 6
        IDENT(K1, K1) = 1.D0      ! Creates Identity Matrix
      END DO

      CALL STIFF_MATRIX (E_MAT, NU_MAT, C_MAT)     ! Matrix Stiffnes
      CALL STIFF_MATRIX (E_INC, NU_INC, C_INC)     ! Inclusion Stiffness
      CALL ESHELBY_SPHERE (NU_MAT, ALPHA, ESHELBY) ! Eshelby tensor
      

!     ------- Calculates the Dilute Strain Concentration Tensor --------
!                            (See Equation 7a)

      CALL INV_MATRIX(C_MAT, COMPL)         ! Compliance Tensor
      
      DO K1=1, 6
        DO K2=1, 6
          MAT_DIFF(K1, K2) = C_INC(K1, K2) - C_MAT(K1, K2)
        END DO
      END DO

      CALL MUL_MATRIX (ESHELBY, COMPL, AUX1)
      CALL MUL_MATRIX (AUX1, MAT_DIFF, AUX2)

      DO K1=1, 6
        AUX2(K1,K1) = AUX2(K1,K1) + IDENT(K1,K1)
      END DO
      
      CALL INV_MATRIX(AUX2, B)      ! B = Strain Concentration Tensor

!     ------------ Calculates the Effective Elastic Tensor -------------
!                            (See Equation 14a)

      DO K1=1, 6
        DO K2=1, 6
          AUX1(K1,K2) = (Vm*IDENT(K1,K2))+(Vf*B(K1,K2))
        END DO
      END DO
      
      CALL INV_MATRIX(AUX1, AUX2)
      CALL MUL_MATRIX(MAT_DIFF, B, AUX1)
      
      DO K1=1, 6
        DO K2=1, 6
          AUX1(K1,K2) = AUX1(K1,K2)*Vf
        END DO
      END DO
      
      CALL MUL_MATRIX(AUX1, AUX2, CONST_MAT)
      
      DO K1=1, 6
        DO K2=1, 6
          CONST_MAT(K1,K2) = CONST_MAT(K1,K2) + C_MAT(K1,K2)
        END DO
      END DO
      
      RETURN
      END SUBROUTINE

      SUBROUTINE VOIGT (E_MAT,E_INC,NU_MAT,NU_INC,Vf,CONST_MAT)
!     ------------------------------------------------------------------
!            Implementation of Voigt Homogenization Scheme
!            Written by Marcelo Medeiros on February 3rd, 2018
!
!     ------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION E_MAT, NU_MAT, E_INC, NU_INC, C_MAT, C_INC
      DOUBLE PRECISION CONST_MAT, Vm, Vf
      INTEGER K1, K2
      DIMENSION C_MAT(6,6), C_INC(6,6), CONST_MAT(6,6) 


      Vm = 1.0D0 - Vf             ! Volume Fraction of the Matrix
      
      CALL STIFF_MATRIX (E_MAT, NU_MAT, C_MAT)     ! Matrix Stiffnes
      CALL STIFF_MATRIX (E_INC, NU_INC, C_INC)     ! Inclusion Stiffness

      DO K1=1, 6
        DO K2=1, 6
          C_MAT(K1, K2) = C_MAT(K1, K2)*Vm
		  C_INC(K1, K2) = C_INC(K1, K2)*Vf
		  CONST_MAT(K1,K2) = C_MAT(K1,K2) + C_INC(K1,K2)
        END DO
      END DO

      RETURN
      END SUBROUTINE      


      SUBROUTINE ESHELBY_SPHERE (NU_MAT, ALPHA, ESHELBY)
!     ------------------------------------------------------------------
!         Implementation of the Eshelby Tensor for Isotropic Spheres
!       Based on "Micromechanics of Defects in Solids", T. Mura (1987)
!          Second Edition,  ISBN 90-247-3343-X,  Page 79, Eq 11.21
!                   Case 1: a1 = a2 = a3 (alhpa = 1)
!     ------------------------------------------------------------------


      DOUBLE PRECISION NU_MAT, ALPHA, ESHELBY
      INTEGER K1, K2
      DIMENSION ESHELBY (6,6)
      
      DO K1=1, 6
        DO K2=1, 6
          ESHELBY(K1,K2) = 0.0D0
        END DO
      END DO

      IF (ALPHA .NE. 1.D0) THEN
        PRINT *, '   ONLY WORKS FOR SPHERICAL INCLUSIONS (ALPHA = 1)'
        RETURN
      END IF
      
      ESHELBY(1,1) = (7.D0-(5.D0*NU_MAT))/(15.D0*(1.D0-NU_MAT))
      ESHELBY(2,2) = ESHELBY(1,1)
      ESHELBY(3,3) = ESHELBY(1,1)
      ESHELBY(1,2) = ((5.D0*NU_MAT)-1)/(15.D0*(1.D0-NU_MAT))
      ESHELBY(2,3) = ESHELBY(1,2);
      ESHELBY(3,1) = ESHELBY(1,2);
      ESHELBY(1,3) = ESHELBY(1,2);
      ESHELBY(2,1) = ESHELBY(1,2);
      ESHELBY(3,2) = ESHELBY(1,2);
      ESHELBY(4,4) = (4.D0-(5.D0*NU_MAT))/(15.D0*(1.D0-NU_MAT))
      ESHELBY(5,5) = ESHELBY(4,4);
      ESHELBY(6,6) = ESHELBY(4,4);
      
      RETURN
      END SUBROUTINE


      SUBROUTINE INV_MATRIX(A, C)
!     ------------------------------------------------------------------
!               Invert Matrix using LU decomposition algorithm.
!                     A is the input and C is the output
!     ------------------------------------------------------------------


      DOUBLE PRECISION A(6,6), C(6,6), L(6,6), U(6,6), B(6), D(6), X(6)
      DOUBLE PRECISION COEFF
      INTEGER K1, K2, K3

      DO K1=1, 6
        B(K1)=0.0D0
        DO K2=1, 6
          L(K1,K2)=0.0D0
          U(K1,K2)=0.0D0
        END DO
      END DO

      DO K1=1, 5
        DO K2=K1+1, 6
          COEFF=A(K2,K1)/A(K1,K1)
          L(K2,K1) = COEFF
            DO K3=K1+1, 6
              A(K2,K3) = A(K2,K3)-COEFF*A(K1,K3)
            END DO
        END DO
      END DO

      DO K1=1, 6
        L(K1,K1) = 1.0D0
        DO K2=1,K1
          U(K2,K1) = A(K2,K1)
        END DO
      END DO

      DO K1=1, 6
        B(K1)=1.0D0
        D(1) = B(1)
        DO K2=2, 6
          D(K2)=B(K2)
          DO K3=1, K2-1
            D(K2) = D(K2) - L(K2,K3)*D(K3)
          END DO
        END DO
        X(6)=D(6)/U(6,6)
        DO K2 = 5,1,-1
          X(K2) = D(K2)
          DO K3=6, K2+1,-1
            X(K2) = X(K2)-U(K2,K3)*X(K3)
          END DO
          X(K2) = X(K2)/U(K2,K2)
        END DO
        DO K2=1, 6
          C(K2,K1) = X(K2)
        END DO
        B(K1)=0.0D0
      END DO

      RETURN
      END SUBROUTINE
      
      
      SUBROUTINE STIFF_MATRIX (E, NU, MATR)
!     ------------------------------------------------------------------
!                  Assembles Isotropic Stiffness Matrix.
!              based on Young's Modulus and Poisson's Ratio
!                 Shear Components in Engineering Notation
!     ------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION E, NU, G, LAMBDA, MATR
      INTEGER K1, K2
      DIMENSION MATR(6,6)
      
      G = E/(2.D0*(1.D0+NU))
      LAMBDA = (E*NU)/((1.D0+NU)*(1.D0-2.D0*NU))
      
      DO K1=1, 3
        DO K2=1, 3
          IF (K1 == K2) THEN
            MATR(K1,K2) = LAMBDA + 2.D0*G
			MATR(K1+3,K2+3) = G
          ELSE
            MATR(K1,K2) = LAMBDA
			MATR(K1+3,K2+3) = 0.0D0
          END IF
		  MATR(K1+3, K2) = 0.0D0
		  MATR(K1, K2+3) = 0.0D0
        END DO
      END DO

      RETURN
      END SUBROUTINE

      
      SUBROUTINE MUL_MATRIX (G, H, R)
!     ------------------------------------------------------------------
!                      Multiplies two 6 by 6 Matrices.
!                G and H are the inputs and R is the output
!     ------------------------------------------------------------------
      INTEGER K1, K2, K3
      DOUBLE PRECISION SUMM, G, H, R
      DIMENSION G(6,6), H(6,6), R(6,6)
      
      DO K1=1, 6
        DO K2=1, 6
          SUMM = 0.0D0
          DO K3=1, 6
            SUMM = SUMM + (G(K1,K3)*H(K3,K2))
          END DO
          R(K1, K2) = SUMM
        END DO
      END DO
    
      RETURN	
      END SUBROUTINE

      SUBROUTINE TURNER (K_MAT, K_INC,CTE_MAT, CTE_INC, Vf, CTE)
!     ------------------------------------------------------------------
!             Calculate the Coeff of Thermal Expansion CTE
!          "Thermal-expansion stresses in reinforced plastics"
!                    DOI: 10.6028/jres.037.015
!     ------------------------------------------------------------------
      DOUBLE PRECISION K_MAT, K_INC,CTE_MAT, CTE_INC, Vf, CTE

      CTE = (K_MAT*CTE_MAT*(1.D0 - Vf) + K_INC*CTE_INC*VF)/             &
	 &        (K_MAT*(1.D0 - Vf) + K_INC*Vf)
    
      RETURN	
      END SUBROUTINE

