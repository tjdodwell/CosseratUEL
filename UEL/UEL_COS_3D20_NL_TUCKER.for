C*************************************************
C**
C** written by Amir Hosein Sakhaei
C**
C** a.sakhaei@exeter.ac.uk  -OR- sakhaei.a@gmail.com 
C**
C** START DATE   : 01/AUGUST/2018
C** last modified: 20/AUGUST/2018
C** University of Exeter
C**
C**     #####
C**    #### _\_  ________
C**    ##=-[.].]| \      \
C**    #(    _\ |  |------|
C**     #   __| |  ||||||||
C**      \  _/  |  ||||||||
C**   .--'--'-. |  | ____ |
C**  / __      `|__|[o__o]|
C**_(____nm_______ /____\____
C*************************************************
C**
C**  Version-1 ::  UEL_COS_3D20_NL_TUCKER
C**  
C**  NOTE 01: THIS VERSION IS THE SECOND VERSION (01 AUGUST 2018)
C**
C*************************************************
C 
       SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     + PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     + KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     + LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
C********************************************************************
C
C      USER SUBROUTINE UEL
C
C********************************************************************
C
      	INCLUDE 'ABA_PARAM.INC'
C
      	DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     + SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     + DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     + JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     + PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
	 
C	 	user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT 
C
C--------------------------------------------------------------------
C		INPUT VARIABLE DEFINITIONS:
C	1)	PROPS: MATERIAL PROPERTIES (READING FROM INPUT FILE)
C	2)	UPDATE: Logical variable: If true, save stress values
C	3)	LTAN: Logical variable: If true, calculate the global stiffness matrix
C	4)	NE: Total number of elements
C	5)	NDOF: Dimension of problem (3)
C	6)	XYZ: (3,NNODE) : Coordinates of all nodes
C	7)	LE: (8,NE) : Element connectivity
C
C**********************************************************************
C             DEFINITION OF PARAMETERS
C**********************************************************************     
C 
C
		INTEGER I, J, K, L, INTN, LX, LY, LZ, COL,INTN_QUAD,LX2, LY2, LZ2
C     
          REAL*8   AIDENT(3,3), PROPERTIES(14), K_CONSTRAIN,  
     +        U_NODES_QUAD(60,1), THETA_NODES_QUAD(24,1),
     +        AMATRX_OLD_QUAD(84,84),RHS_OLD_QUAD(84,1),
     +        AMATRX_QUAD(60,60), RHS_QUAD(60,1),
     +        AMATRX_LIN(24,24), RHS_LIN(24,1),  
     +        AMATRX_LIN_PEN(84,84), RHS_LIN_PEN(84,1),
     +        AMATRX_NEW_QUAD(84,84), RHS_NEW_QUAD(84,1)

C         
		PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0, 
     +          ONE_FOURTH=1.D0/4.D0, NINE=9.D0 )
C   
		CALL ONEM(AIDENT)
          
!          PRINT *, "$$$$$$  I AM HERE AT TIME = $$$$$$" , TIME(1), 
!     +                     " AND TIME STEP = ", DTIME
          
!          PRINT *, "MLVARX" , MLVARX
!          PRINT *, "NDOFEL" , NDOFEL
!          PRINT *, "AIDENT = "
!          CALL PRTMAT(AIDENT,3,3)
        
          PROPERTIES(1) = PROPS(1)
          PROPERTIES(2) = PROPS(2)
          PROPERTIES(3) = PROPS(3)
          PROPERTIES(4) = PROPS(4)
          PROPERTIES(5) = PROPS(5)
          PROPERTIES(6) = PROPS(6)
          PROPERTIES(7) = PROPS(7)
          PROPERTIES(8) = PROPS(8)
          PROPERTIES(9) = PROPS(9)
          PROPERTIES(10) = PROPS(10)
          PROPERTIES(11) = PROPS(11)
          PROPERTIES(12) = PROPS(12)
          PROPERTIES(13) = PROPS(13)
          PROPERTIES(14) = PROPS(14)
        
          K_CONSTRAIN = PROPS(15)

!C
C%***********************************************************************
C%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX AND RESIDUAL FORCE FOR
C%  HYPERELASTIC MATERIAL MODELS
C%***********************************************************************
C        
C		MAKE 20*3 ELEMENT DOF (DISPLACEMENT) MATRIX ELDPS FROM U(NDOFEL) 
C		HERE 	ELDPS(20,3) AND U(84)
C
          CALL U_SEPARATOR_C3D20_COS(U,U_NODES_QUAD,THETA_NODES_QUAD)
          
		DO I = 1, 84
    	    DO J = 1, 84
	        AMATRX(I,J) = ZERO
	        RHS(I,1) = ZERO
              AMATRX_OLD_QUAD(I,J) = ZERO
	        RHS_OLD_QUAD(I,1) = ZERO
		END DO
		END DO
          
          CALL QUAD_PART_1(COORDS,U_NODES_QUAD,PROPERTIES, 
     +       AMATRX_QUAD,RHS_QUAD)
          
          CALL LIN_PART_2(COORDS,U_NODES_QUAD,THETA_NODES_QUAD,
     +        PROPERTIES,K_CONSTRAIN,AMATRX_LIN,RHS_LIN,AMATRX_LIN_PEN,
     +        RHS_LIN_PEN)          
C
          !print *, 'RHS_QUAD:'
          !CALL PRTMAT(RHS_QUAD,60,1)
          
          !print *, 'AMATRX_QUAD:'
          !CALL PRTMAT(AMATRX_QUAD,60,60)
          
          !print *, 'RHS_LIN:'
          !CALL PRTMAT(RHS_LIN,24,1)
          !
          !print *, 'RHS_LIN_PEN:'
          !CALL PRTMAT(RHS_LIN_PEN,84,1)
          
          DO I = 1, 60
    	    DO J = 1, 60
              AMATRX_OLD_QUAD(I,J) = AMATRX_QUAD(I,J)
		END DO
		END DO
          
          DO I = 1, 60
	        RHS_OLD_QUAD(I,1) = RHS_QUAD(I,1)
		END DO
C   
          DO I = 1, 24
    	    DO J = 1, 24
              AMATRX_OLD_QUAD(I+60,J+60) = AMATRX_LIN(I,J)
		END DO
		END DO
          
          DO I = 1, 24
	        RHS_OLD_QUAD(I+60,1) = RHS_LIN(I,1)
		END DO
C          
          DO I = 1, 84
    	    DO J = 1, 84
              AMATRX_OLD_QUAD(I,J) = AMATRX_OLD_QUAD(I,J) + 
     +            AMATRX_LIN_PEN(I,J)
		END DO
		END DO
          
          DO I = 1, 84
	        RHS_OLD_QUAD(I,1) = RHS_OLD_QUAD(I,1) + RHS_LIN_PEN(I,1)
		END DO
          
C          
C         NEED TO REARRANGE THE AMATRX AND RHS MATRIXS              
C          
C         !!!! WE ASSUNE THE FOLLOWING SUBROUTINE IS RIGHT !!!
C
          CALL AMATRX_RHS_CONVERSION_QUAD(AMATRX_OLD_QUAD,RHS_OLD_QUAD,
     +         AMATRX_NEW_QUAD,RHS_NEW_QUAD)
C
          
              DO I = 1, 84
              DO J = 1, 84
	    	    AMATRX(I,J) = AMATRX_NEW_QUAD(I,J)
			END DO
              END DO
          
              DO I = 1, 84
	    	    RHS(I,1) = RHS_NEW_QUAD(I,1)
			END DO
              
!          print *, 'RHS:'
!          CALL PRTMAT(RHS,84,1)
          
		RETURN
		END
				
C
C******************************************************************************************************
C******************************************************************************************************
C******************************************************************************************************
C******************************************************************************************************
C******************************************************************************************************
C******************************************************************************************************
C
C**************************************************************
C	SUBROUTINE APPLIED IN THIS ROUTINE
C
C**************************************************************
C
C          
C%*************************************************************************
C% SEPARATING TO CALCULATE AMATRX AND RHS IN QUAD DISPLACEMENT 20 NODE
C%*************************************************************************
C%%            
          SUBROUTINE QUAD_PART_1(COORDS,U_NODES_QUAD,PROPERTIES, 
     +       AMATRX_QUAD,RHS_QUAD)
C
      
		INTEGER I, J, K, INTN_QUAD, LX, LY, LZ
C     
		REAL*8   AIDENT(3,3),
     +        COORDS(3,20), U_NODES_QUAD(60,1), PROPERTIES(14),    
     +        XG_QUAD(1,3), WGT_QUAD(1,3), ELXYZ_QUAD(20,3), 
     +        ELDPS_QUAD(20,3), AMATRX_QUAD(60,60), RHS_QUAD(60,1),
     +        E1, E2, E3, XI(1,3), DETGJ_QUAD, GJ_QUAD(3,3), FAC,
     +        SHPD_QUAD(3,20), SHPD_TRANS_QUAD(20,3), DUDX_QUAD(3,3),
     +        ELDPS_TRANS_QUAD(3,20), F_QUAD(3,3), STRESS_VECT(6),
     +        KAPPA_CURVE_QUAD(3,3), MATERIAL_MATRIX(6,6), 
     +        COSSERAT_MATRIX(9,9), DWDKAPPA(3,3), BN_QUAD(6,60), 
     +        BG_QUAD(9,60), STRESS_MATRIX(3,3), SHEAD(9,9), 
     +        BN_TRANS_QUAD(60,6), BG_TRANS_QUAD(60,9), Z1_QUAD(60,6),
     +        EKF1_QUAD(60,60), EKF2_QUAD(60,60), Z2_QUAD(60,9),
     +        EKF_QUAD(60,60), RF1_QUAD(60,1)

C         
		PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, SIX=6.D0, ONE_FOURTH=1.D0/4.D0 , 
     +          FIVE=5.D0, EIGHT=8.D0, NINE=9.D0)          
C
C
C%***********************************************************************
C%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX AND RESIDUAL FORCE FOR
C%  HYPERELASTIC MATERIAL MODELS , 3D AND 20 NODES
C%***********************************************************************
C
C		global DISPTD FORCE GKF SIGMA ----> CONVERT TO FORTRAN
C
C		 Integration points and weights
  		XG_QUAD(1,1) = -0.7745966692414D0
  		XG_QUAD(1,2) = ZERO
  		XG_QUAD(1,3) = 0.7745966692414D0
  		WGT_QUAD(1,1) = FIVE/NINE
  		WGT_QUAD(1,2) = EIGHT/NINE
  		WGT_QUAD(1,3) = FIVE/NINE
          
          CALL ONEM(AIDENT)
C
C		 MAKE 20*3 ELEMENT COORDINATE MATRIX 		
  		DO I=1, 20
		DO J=1, 3
			ELXYZ_QUAD(I,J) = COORDS(J,I)	
		END DO
		END DO
C
		DO I=1, 20
		DO J=1, 3
			ELDPS_QUAD(I,J) = U_NODES_QUAD(3*(I-1)+J,1)	
		END DO
		END DO          
C
C	 Index for history variables (each integration pt)
  		INTN_QUAD=0          
C
          DO I = 1, 60
    	    DO J = 1, 60
	        AMATRX_QUAD(I,J) = ZERO
	        RHS_QUAD(I,1) = ZERO
		END DO
		END DO
C          
C
C		%LOOP OVER INTEGRATION POINTS
C
		DO LX=1, 3
		DO LY=1, 3
		DO LZ=1, 3
C
			E1 = XG_QUAD(1,LX)
			E2 = XG_QUAD(1,LY)
			E3 = XG_QUAD(1,LZ)
			XI(1,1) = E1 
			XI(1,2) = E2
			XI(1,3) = E3
C			
			INTN_QUAD = INTN_QUAD + 1
!              PRINT *, 'INTN_QUAD = ', INTN_QUAD
C
              CALL SHAPEL_HEX20(XI, ELXYZ_QUAD, DETGJ_QUAD, GJ_QUAD,
     +             SHPD_QUAD)
C
			FAC = WGT_QUAD(1,LX)*WGT_QUAD(1,LY)*WGT_QUAD(1,LZ)*DETGJ_QUAD
              
              !print *, 'DETGJ_QUAD:' , DETGJ_QUAD
              !
              !print *, 'GJ_QUAD:'
              !CALL PRTMAT(GJ_QUAD,3,3)
              !
              !print *, 'SHPD_QUAD:'
              !CALL PRTMAT(SHPD_QUAD,3,20)
C              
C		MAKE SHPD_TRANS(20,3) && ELDPS_TRANS(3,20)
			DO I = 1, 20
    	        DO J = 1, 3
	    	    SHPD_TRANS_QUAD(I,J) = SHPD_QUAD(J,I)
			END DO
			END DO
C
			DO I = 1, 3
    	        DO J = 1, 20
	    	    ELDPS_TRANS_QUAD(I,J) = ELDPS_QUAD(J,I)
			END DO
			END DO
C
C		MAKE Displacement gradient DU/DX = ELDPS_TRANS(3,20)*SHPD_TRANS(20,3)
C			
			DO I = 1, 3
    	        DO J = 1, 3
	    	    DUDX_QUAD(I,J) = ZERO
	    	    DO K = 1, 20
	    	    	DUDX_QUAD(I,J) = DUDX_QUAD(I,J) + 
     +                    ELDPS_TRANS_QUAD(I,K)*SHPD_TRANS_QUAD(K,J)
	    	    END DO
			END DO
			END DO
      
C
C		MAKE DEFORMATION gradient F(3,3) = DUDX(I,J) + AIDENT(I,J)
C
			DO I = 1, 3
    	        DO J = 1, 3
	    	    F_QUAD(I,J) = DUDX_QUAD(I,J) + AIDENT(I,J)
			END DO
			END DO
C
C
C		% Compute stress and tangent stiffness BY CALLING THE UMAT
C
              CALL ZEROM(KAPPA_CURVE_QUAD)
              
              !print *, 'SHPD_QUAD:'
              !CALL PRTMAT(SHPD_QUAD,3,20)
              !
              !print *, 'ELDPS_QUAD:'
              !CALL PRTMAT(ELDPS_QUAD,20,3)
              !
              !print *, 'SHPD_TRANS_QUAD:'
              !CALL PRTMAT(SHPD_TRANS_QUAD,20,3)
              !
              !print *, 'ELDPS_TRANS_QUAD:'
              !CALL PRTMAT(ELDPS_TRANS_QUAD,3,20)
              
              !print *, 'AIDENT:'
              !CALL PRTMAT(AIDENT,3,3)
              !
              !print *, 'DUDX_QUAD:'
              !CALL PRTMAT(DUDX_QUAD,3,3)
              !
              !print *, 'F_QUAD:'
              !CALL PRTMAT(F_QUAD,3,3)
              
              !print *, 'KAPPA_CURVE_QUAD:'
              !CALL PRTMAT(KAPPA_CURVE_QUAD,3,3)
              
              CALL UMAT_COSSERAT_TUCKER(PROPERTIES, F_QUAD, 
     +         KAPPA_CURVE_QUAD, STRESS_VECT, MATERIAL_MATRIX,
     +         COSSERAT_MATRIX, DWDKAPPA)
C
C ---> WE CAN ALSO STORE STRAIN AND STRESS AND USED THEM LATER FOR POST PROCESSING
C
C		CALCULATE "TANGENT STIFFNESS MATRIX" AND "RESIDUAL FORCE"
C			
C
C		CALCULATE BN NA D BG   ---> NEEDS MORE UNDERSTANDING
C			
			CALL MAKER_BNBG_HEX20(SHPD_QUAD,F_QUAD,BN_QUAD,BG_QUAD)
C
C		CALCULATE TANGENT STIFFNESS
C
			STRESS_MATRIX(1,1) = STRESS_VECT(1)
			STRESS_MATRIX(2,2) = STRESS_VECT(2)
			STRESS_MATRIX(3,3) = STRESS_VECT(3)
			
			STRESS_MATRIX(1,2) = STRESS_VECT(4)
			STRESS_MATRIX(2,1) = STRESS_VECT(4)
			
			STRESS_MATRIX(2,3) = STRESS_VECT(5)
			STRESS_MATRIX(3,2) = STRESS_VECT(5)
			
			STRESS_MATRIX(1,3) = STRESS_VECT(6)
			STRESS_MATRIX(3,1) = STRESS_VECT(6)
C			
			DO I = 1, 9
    	        DO J = 1, 9
	    	    SHEAD(I,J) = ZERO
			END DO
			END DO            

              DO I = 1, 3
    	        DO J = 1, 3
	    	    SHEAD(I,J) = STRESS_MATRIX(I,J)
	    	    SHEAD(I+3,J+3) = STRESS_MATRIX(I,J)
	    	    SHEAD(I+6,J+6) = STRESS_MATRIX(I,J)
			END DO
			END DO
C
C		MAKE BN_TRANS(60,6) && BG_TRANS(60,9)
C
			DO I = 1, 60
    	        DO J = 1, 6
	    	    BN_TRANS_QUAD(I,J) = BN_QUAD(J,I)
			END DO
			END DO
C			
			DO I = 1, 60
    	        DO J = 1, 9
	    	    BG_TRANS_QUAD(I,J) = BG_QUAD(J,I)
			END DO
			END DO
C					
			CALL MPROD_GEN(BN_TRANS_QUAD,60,6,MATERIAL_MATRIX,6,6,Z1_QUAD)
			CALL MPROD_GEN(Z1_QUAD,60,6,BN_QUAD,6,60,EKF1_QUAD)
C			
			CALL MPROD_GEN(BG_TRANS_QUAD,60,9,SHEAD,9,9,Z2_QUAD)
			CALL MPROD_GEN(Z2_QUAD,60,9,BG_QUAD,9,60,EKF2_QUAD)
C              
C
			DO I = 1, 60
    	        DO J = 1, 60
	    	    EKF_QUAD(I,J) = EKF1_QUAD(I,J) + EKF2_QUAD(I,J)
			END DO
			END DO             
C
C		CALCULATE Tangent stiffness BY INTEGRATION
C
			DO I = 1, 60
    	        DO J = 1, 60
	    	    AMATRX_QUAD(I,J) = AMATRX_QUAD(I,J) + FAC*EKF_QUAD(I,J)
			END DO
			END DO              
C
C		CALCULATE RESIDUAL FORCE BY INTEGRATION
C
			CALL MPROD_GEN(BN_TRANS_QUAD,60,6,STRESS_VECT,6,1,RF1_QUAD)
C              
			DO I = 1, 60
	    	    RHS_QUAD(I,1) = RHS_QUAD(I,1) - FAC*RF1_QUAD(I,1)
			END DO
C
          END DO
		END DO
		END DO
         
C
          RETURN
          END         
C          
C%*************************************************************************
C% SEPARATING TO CALCULATE AMATRX AND RHS IN LINEAR ROTATION 8 NODE
C%*************************************************************************
C%%       
          
          SUBROUTINE LIN_PART_2(COORDS,U_NODES_QUAD,THETA_NODES_QUAD,
     +        PROPERTIES,K_CONSTRAIN,AMATRX_LIN,RHS_LIN,AMATRX_LIN_PEN,
     +        RHS_LIN_PEN)
C
      
		INTEGER I, J, K, INTN, LX, LY, LZ
C     
		REAL*8   AIDENT(3,3), K_CONSTRAIN, THETA_NODES_QUAD(24,1),
     +        COORDS(3,20), U_NODES_QUAD(60,1), PROPERTIES(14),     
     +        XG(1,2), WGT(1,2), ELXYZ_QUAD(20,3), ELXYZ_LIN(8,3), 
     +        ELDPS_QUAD(20,3), ELROT_LIN(8,3),     
     +        AMATRX_LIN(24,24), RHS_LIN(24,1),
     +        E1, E2, E3, XI(1,3), DETGJ, GJ(3,3), SHPD(3,8), SF(8,1),
     +        DETGJ_QUAD, GJ_QUAD(3,3), FAC,
     +        SHPD_QUAD(3,20), SHPD_TRANS(8,3), ELROT_TRANS_LIN(3,8),
     +        SHPD_TRANS_QUAD(20,3), DUDX_QUAD(3,3), DTHETADX(3,3),
     +        ELDPS_TRANS_QUAD(3,20), F_QUAD(3,3), STRESS_VECT(6),
     +        KAPPA_CURVE(3,3), MATERIAL_MATRIX(6,6), 
     +        COSSERAT_MATRIX(9,9), DWDKAPPA(3,3), BN_QUAD(6,60), 
     +        BG_QUAD(9,60), STRESS_MATRIX(3,3), SHEAD(9,9), 
     +        BN_TRANS_QUAD(60,6), BG_TRANS_QUAD(60,9), Z1_QUAD(60,6),
     +        EKF1_QUAD(60,60), EKF2_QUAD(60,60), Z2_QUAD(60,9),
     +        EKF_QUAD(60,60), RF1_QUAD(60,1), Z2(24,9),
     +        B1_COS(9,24), B2_COS(9,24), B3_COS(3,24), B4_COS(3,24),
     +		B5_COS(9,24), B6_COS(9,24), B1_COS_TRANS(24,9),
     +		B2_COS_TRANS(24,9), B3_COS_TRANS(24,3), B4_COS_TRANS(24,3),
     +		B5_COS_TRANS(24,9), B6_COS_TRANS(24,9), Z7_CHECKER(3,3),
     +        B4_COS_QUAD(3,60), B4_COS_TRANS_QUAD(60,3), Z6,
     +        B_PEN_QUAD(3,84), B_PEN_TRANS_QUAD(84,3), F_INV_NEW(3,3),
     +        THETA_VEC2(24,1), THETA_CHECK(3,1), ERROR_CONST(3,1),
     +        EKF_COS(24,24), Z5_QUAD(84,84), K_PEN_QUAD(84,84), 
     +        AMATRX_LIN_PEN(84,84), COS_STRESS_VEC(9,1), RF1_COS(24,1),
     +        U_PEN_QUAD(84,1), RF_PEN_QUAD(84,1), RHS_LIN_PEN(84,1)

C         
		PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, SIX=6.D0, ONE_FOURTH=1.D0/4.D0 , 
     +          FIVE=5.D0, EIGHT=8.D0, NINE=9.D0)          
C
C
C%***********************************************************************
C%  MAIN PROGRAM COMPUTING GLOBAL STIFFNESS MATRIX AND RESIDUAL FORCE FOR
C%  HYPERELASTIC MATERIAL MODELS , 3D AND 20 NODES
C%***********************************************************************
C		global DISPTD FORCE GKF SIGMA ----> CONVERT TO FORTRAN
C
C		 Integration points and weights
  		XG(1,1) = -0.57735026918963D0
  		XG(1,2) = 0.57735026918963D0
  		WGT(1,1) = ONE
  		WGT(1,2) = ONE
          
          CALL ONEM(AIDENT)
C
C		 MAKE 20*3 ELEMENT COORDINATE MATRIX 
C
  		DO I=1, 20
		DO J=1, 3
			ELXYZ_QUAD(I,J) = COORDS(J,I)	
		END DO
		END DO
C
          DO J=1, 3
			ELXYZ_LIN(1,J) = ELXYZ_QUAD(1,J)
              ELXYZ_LIN(2,J) = ELXYZ_QUAD(3,J)
              ELXYZ_LIN(3,J) = ELXYZ_QUAD(5,J)
              ELXYZ_LIN(4,J) = ELXYZ_QUAD(7,J)
              ELXYZ_LIN(5,J) = ELXYZ_QUAD(13,J)
              ELXYZ_LIN(6,J) = ELXYZ_QUAD(15,J)
              ELXYZ_LIN(7,J) = ELXYZ_QUAD(17,J)
              ELXYZ_LIN(8,J) = ELXYZ_QUAD(19,J)
		END DO
C
          DO I=1, 20
		DO J=1, 3
			ELDPS_QUAD(I,J) = U_NODES_QUAD(3*(I-1)+J,1)	
		END DO
		END DO
C
          DO I=1, 8
		DO J=1, 3
              ELROT_LIN(I,J) = THETA_NODES_QUAD(3*(I-1)+J,1)
		END DO
		END DO
C	 Index for history variables (each integration pt)
  		INTN=0
C
          DO I = 1, 24
    	    DO J = 1, 24
	        AMATRX_LIN(I,J) = ZERO
	        RHS_LIN(I,1) = ZERO
		END DO
		END DO 
          
          DO I = 1, 84
    	    DO J = 1, 84
	        AMATRX_LIN_PEN(I,J) = ZERO
	        RHS_LIN_PEN(I,1) = ZERO
		END DO
		END DO
C
C		%LOOP OVER INTEGRATION POINTS
C
		DO LX=1, 2
		DO LY=1, 2
		DO LZ=1, 2
C
			E1 = XG(1,LX)
			E2 = XG(1,LY)
			E3 = XG(1,LZ)
			XI(1,1) = E1 
			XI(1,2) = E2
			XI(1,3) = E3
C			
			INTN = INTN + 1
!              PRINT *, 'INTN = ', INTN
C       %
C       % Determinant and shape function derivatives
C		SHPD(3,8) = GLOBAL DERIVATION OF SHAPE FUNCTIONS
C		DETGJ	= DETERMINANT OF JACOUBIAN MATRIX
C		GJ(3,3) = JACOUBIAN MATRIX
C		SF(8,1) = SHAPE FUNCTIONS N1 TO N8 AT EACH GAUSS POINT
C		%             
			CALL SHAPEL_HEX8(XI, ELXYZ_LIN, DETGJ, GJ, SHPD, SF)  
              
              FAC = WGT(1,LX)*WGT(1,LY)*WGT(1,LZ)*DETGJ
C
              CALL SHAPEL_HEX20(XI, ELXYZ_QUAD, DETGJ_QUAD, GJ_QUAD,
     +             SHPD_QUAD)
C
C		MAKE SHPD_TRANS(8,3) && ELDPS_COS_TRANS(3,8) && ELROT_COS_TRANS(3,8)
C            
			DO I = 1, 8
    	        DO J = 1, 3
	    	    SHPD_TRANS(I,J) = SHPD(J,I)
			END DO
			END DO
              
C		MAKE SHPD_TRANS(20,3) && ELDPS_TRANS(3,20)
			DO I = 1, 20
    	        DO J = 1, 3
	    	    SHPD_TRANS_QUAD(I,J) = SHPD_QUAD(J,I)
			END DO
			END DO
C
C
			DO I = 1, 3
    	        DO J = 1, 20
	    	    ELDPS_TRANS_QUAD(I,J) = ELDPS_QUAD(J,I)
			END DO
			END DO
C
			DO I = 1, 3
    	        DO J = 1, 8
	    	    ELROT_TRANS_LIN(I,J) = ELROT_LIN(J,I)
			END DO
			END DO
C
C		MAKE Displacement gradient DU/DX = ELDPS_TRANS(3,20)*SHPD_TRANS(20,3)
C			
			DO I = 1, 3
    	        DO J = 1, 3
	    	    DUDX_QUAD(I,J) = ZERO
	    	    DO K = 1, 20
	    	    	DUDX_QUAD(I,J) = DUDX_QUAD(I,J) + 
     +                    ELDPS_TRANS_QUAD(I,K)*SHPD_TRANS_QUAD(K,J)
	    	    END DO
			END DO
			END DO
C
C		MAKE DEFORMATION gradient F(3,3) = DUDX(I,J) + AIDENT(I,J)
C
			DO I = 1, 3
    	        DO J = 1, 3
	    	    F_QUAD(I,J) = DUDX_QUAD(I,J) + AIDENT(I,J)
			END DO
			END DO              
C
C		MAKE Displacement gradient DU/DX = ELDPS_TRANS(3,8)*SHPD_TRANS(8,3)
C		&& 	ROTATION GRADIENT	DTHETADX = ELROT_COS_TRANS(3,8)*SHPD_TRANS(8,3)
C			
			DO I = 1, 3
    	        DO J = 1, 3
              DTHETADX(I,J) = ZERO
	    	    DO K = 1, 8
	    	    	DTHETADX(I,J) = DTHETADX(I,J) + 
     +                        ELROT_TRANS_LIN(I,K)*SHPD_TRANS(K,J)
	    	    END DO
			END DO
			END DO              
C
C		MAKE THE CURVATURE TENSOR KAPPA_CURVE(3,3) = DTHETADX(I,J)
C			
			DO I = 1, 3
    	        DO J = 1, 3
	    	    KAPPA_CURVE(I,J) = DTHETADX(I,J) 
			END DO
			END DO              
C
C		% Compute stress and tangent stiffness BY CALLING THE UMAT
C             
              CALL UMAT_COSSERAT_TUCKER(PROPERTIES, F_QUAD, KAPPA_CURVE,
     +       STRESS_VECT, MATERIAL_MATRIX, COSSERAT_MATRIX, DWDKAPPA)
C
C ---> WE CAN ALSO STORE STRAIN AND STRESS AND USED THEM LATER FOR POST PROCESSING
C
C		CALCULATE "TANGENT STIFFNESS MATRIX" AND "RESIDUAL FORCE"
C
C
C		CALCULATE B MATRIX FOR COSSERAT MODEL
C			
			CALL MAKER_B_MATRIXS_COS_HEX8(SF,SHPD,F_QUAD,B1_COS,B2_COS,B3_COS,
     +		  B4_COS,B5_COS,B6_COS,B1_COS_TRANS,B2_COS_TRANS,B3_COS_TRANS,
     +		  B4_COS_TRANS,B5_COS_TRANS,B6_COS_TRANS) 
              
              CALL MAKER_B4_COS_QUAD_HEX20(SHPD_QUAD,F_QUAD,
     +		  B4_COS_QUAD, B4_COS_TRANS_QUAD)               
        
C
              DO I=1,3  
              DO J=1,60               
                  B_PEN_QUAD(I,J) = ONE_HALF*B4_COS_QUAD(I,J) 
              END DO
              END DO
              
              DO I=1,3  
              DO J=1,24  
                  B_PEN_QUAD(I,J+60) = B3_COS(I,J)         
              END DO
              END DO
              
              DO I = 1, 84
    	        DO J = 1, 3
	    	    B_PEN_TRANS_QUAD(I,J) = B_PEN_QUAD(J,I)
			END DO
			END DO
C
CCCCCCCCC CHECK CONSTRAIN CCCCCCCCCCC
              CALL MATINV(F_QUAD,F_INV_NEW,Z6)
              
              CALL MPROD(DUDX_QUAD,F_INV_NEW,Z7_CHECKER)
              
              DO I = 1, 3
    	        DO J = 1, 3
	    	    Z7_CHECKER(I,J) = ONE_HALF*Z7_CHECKER(I,J)
			END DO
			END DO
              
              DO I=1, 24
                  THETA_VEC2(I,1) = THETA_NODES_QUAD(I,1)
		    END DO
               
              CALL MPROD_GEN(B3_COS,3,24,THETA_VEC2,24,1,THETA_CHECK)
              
              ERROR_CONST(1,1) = THETA_CHECK(1,1) + Z7_CHECKER(2,3)
     +                     - Z7_CHECKER(3,2)
              
              ERROR_CONST(2,1) = THETA_CHECK(2,1) + Z7_CHECKER(3,1)
     +                     - Z7_CHECKER(1,3)
              
              ERROR_CONST(3,1) = THETA_CHECK(3,1) + Z7_CHECKER(1,2)
     +                     - Z7_CHECKER(2,1)
              
              !PRINT *, "THETA_CHECK = "
              !CALL PRTMAT(THETA_CHECK,3,1)
              !PRINT *, "Z7_CHECKER = "
              !CALL PRTMAT(Z7_CHECKER,3,3)
              !PRINT *, "DUDX_QUAD = "
              !CALL PRTMAT(DUDX_QUAD,3,3)
              !PRINT *, "ERROR_CONST = "
              !CALL PRTMAT(ERROR_CONST,3,1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC              
C
C		CALCULATING COSSERAT PART OF THE STIFFNESS MATRIX
C			
			CALL MPROD_GEN(B2_COS_TRANS,24,9,COSSERAT_MATRIX,9,9,Z2)
			CALL MPROD_GEN(Z2,24,9,B1_COS,9,24,EKF_COS)
C
C		CALCULATING CONSTRAINT PART OF THE STIFFNESS MATRIX
C
C         NEW PENALTY CONSTRAIN
C
              CALL MPROD_GEN(B_PEN_TRANS_QUAD,84,3,B_PEN_QUAD,
     +                3,84,Z5_QUAD)
              
              DO I = 1, 84
    	        DO J = 1, 84
	    	    K_PEN_QUAD(I,J) =  K_CONSTRAIN*Z5_QUAD(I,J)
			END DO
			END DO              
C
              DO I = 1, 84    
    	        DO J = 1, 84         
                  AMATRX_LIN_PEN(I,J) = AMATRX_LIN_PEN(I,J) +
     +             FAC*K_PEN_QUAD(I,J)
              END DO
              END DO
C
              DO I = 1, 24
    	        DO J = 1, 24                 
                  AMATRX_LIN(I,J) = AMATRX_LIN(I,J) + FAC*EKF_COS(I,J)
              END DO
              END DO
C
C		CALCULATE RESIDUAL FORCE BY INTEGRATION (COSSERAT PART)
C
			COS_STRESS_VEC(1,1) = DWDKAPPA(1,1)
			COS_STRESS_VEC(2,1) = DWDKAPPA(2,2)
			COS_STRESS_VEC(3,1) = DWDKAPPA(3,3)
			COS_STRESS_VEC(4,1) = DWDKAPPA(1,2)
			COS_STRESS_VEC(5,1) = DWDKAPPA(2,1)
			COS_STRESS_VEC(6,1) = DWDKAPPA(2,3)
			COS_STRESS_VEC(7,1) = DWDKAPPA(3,2)
			COS_STRESS_VEC(8,1) = DWDKAPPA(1,3)
			COS_STRESS_VEC(9,1) = DWDKAPPA(3,1)
              
              CALL MPROD_GEN(B2_COS_TRANS,24,9,COS_STRESS_VEC,9,1,
     +                RF1_COS)
C
			DO I = 1, 24
	    	    RHS_LIN(I,1) = RHS_LIN(I,1) - FAC*RF1_COS(I,1)
			END DO              
C
C		CALCULATE RESIDUAL FORCE BY INTEGRATION (CONSTRAINT PART)
C 
              DO I=1,60
                  U_PEN_QUAD(I,1) = U_NODES_QUAD(I,1)
              END DO
              
              DO I=1,24
                  U_PEN_QUAD(I+60,1) = THETA_NODES_QUAD(I,1)
              END DO
              
              CALL MPROD_GEN(K_PEN_QUAD,84,84,U_PEN_QUAD,84,1,
     +                             RF_PEN_QUAD)
              
              DO I = 1, 84
	    	    RHS_LIN_PEN(I,1) = RHS_LIN_PEN(I,1) - FAC*RF_PEN_QUAD(I,1)
			END DO
C
          END DO
		END DO
		END DO

C
          RETURN
          END
C          
C%*************************************************************************
C% SEPARATING U VECTOR IN QUADRATIC COSERAT ELEMENT
C%*************************************************************************
C%%             
      SUBROUTINE U_SEPARATOR_C3D20_COS(U,U_NODES_QUAD,THETA_NODES_QUAD)
      
      IMPLICIT REAL*8 (A-H,O-Z)
        
      INTEGER I, J, K
          
      REAL*8 U(84,1), U_NODES_QUAD(60,1), THETA_NODES_QUAD(24,1)
      
      PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )
        
      U_NODES_QUAD(1,1) = U(1,1)  ! NODE 1
      U_NODES_QUAD(2,1) = U(2,1)
      U_NODES_QUAD(3,1) = U(3,1)
      THETA_NODES_QUAD(1,1) = U(4,1)
      THETA_NODES_QUAD(2,1) = U(5,1)
      THETA_NODES_QUAD(3,1) = U(6,1)
      
      U_NODES_QUAD(4,1) = U(7,1)  ! NODE 2
      U_NODES_QUAD(5,1) = U(8,1)
      U_NODES_QUAD(6,1) = U(9,1)
      
      U_NODES_QUAD(7,1) = U(10,1)  ! NODE 3
      U_NODES_QUAD(8,1) = U(11,1)
      U_NODES_QUAD(9,1) = U(12,1)
      THETA_NODES_QUAD(4,1) = U(13,1)
      THETA_NODES_QUAD(5,1) = U(14,1)
      THETA_NODES_QUAD(6,1) = U(15,1)
      
      U_NODES_QUAD(10,1) = U(16,1)   ! NODE 4
      U_NODES_QUAD(11,1) = U(17,1)
      U_NODES_QUAD(12,1) = U(18,1)
      
      U_NODES_QUAD(13,1) = U(19,1)  ! NODE 5
      U_NODES_QUAD(14,1) = U(20,1)
      U_NODES_QUAD(15,1) = U(21,1)
      THETA_NODES_QUAD(7,1) = U(22,1)
      THETA_NODES_QUAD(8,1) = U(23,1)
      THETA_NODES_QUAD(9,1) = U(24,1)
      
      U_NODES_QUAD(16,1) = U(25,1)   ! NODE 6
      U_NODES_QUAD(17,1) = U(26,1)
      U_NODES_QUAD(18,1) = U(27,1)
      
      U_NODES_QUAD(19,1) = U(28,1)  ! NODE 7
      U_NODES_QUAD(20,1) = U(29,1)
      U_NODES_QUAD(21,1) = U(30,1)
      THETA_NODES_QUAD(10,1) = U(31,1)
      THETA_NODES_QUAD(11,1) = U(32,1)
      THETA_NODES_QUAD(12,1) = U(33,1)
      
      U_NODES_QUAD(22,1) = U(34,1)   ! NODE 8
      U_NODES_QUAD(23,1) = U(35,1)
      U_NODES_QUAD(24,1) = U(36,1)
      
      U_NODES_QUAD(25,1) = U(37,1)   ! NODE 9
      U_NODES_QUAD(26,1) = U(38,1)
      U_NODES_QUAD(27,1) = U(39,1)
      
      U_NODES_QUAD(28,1) = U(40,1)   ! NODE 10
      U_NODES_QUAD(29,1) = U(41,1)
      U_NODES_QUAD(30,1) = U(42,1)
      
      U_NODES_QUAD(31,1) = U(43,1)   ! NODE 11
      U_NODES_QUAD(32,1) = U(44,1)
      U_NODES_QUAD(33,1) = U(45,1)
      
      U_NODES_QUAD(34,1) = U(46,1)   ! NODE 12
      U_NODES_QUAD(35,1) = U(47,1)
      U_NODES_QUAD(36,1) = U(48,1)
      
      U_NODES_QUAD(37,1) = U(49,1)   ! NODE 13
      U_NODES_QUAD(38,1) = U(50,1)
      U_NODES_QUAD(39,1) = U(51,1)
      THETA_NODES_QUAD(13,1) = U(52,1)
      THETA_NODES_QUAD(14,1) = U(53,1)
      THETA_NODES_QUAD(15,1) = U(54,1)
      
      U_NODES_QUAD(40,1) = U(55,1)   ! NODE 14
      U_NODES_QUAD(41,1) = U(56,1)
      U_NODES_QUAD(42,1) = U(57,1)
      
      U_NODES_QUAD(43,1) = U(58,1)   ! NODE 15
      U_NODES_QUAD(44,1) = U(59,1)
      U_NODES_QUAD(45,1) = U(60,1)
      THETA_NODES_QUAD(16,1) = U(61,1)
      THETA_NODES_QUAD(17,1) = U(62,1)
      THETA_NODES_QUAD(18,1) = U(63,1)
      
      U_NODES_QUAD(46,1) = U(64,1)   ! NODE 16
      U_NODES_QUAD(47,1) = U(65,1)
      U_NODES_QUAD(48,1) = U(66,1)
      
      U_NODES_QUAD(49,1) = U(67,1)   ! NODE 17
      U_NODES_QUAD(50,1) = U(68,1)
      U_NODES_QUAD(51,1) = U(69,1)
      THETA_NODES_QUAD(19,1) = U(70,1)
      THETA_NODES_QUAD(20,1) = U(71,1)
      THETA_NODES_QUAD(21,1) = U(72,1)
      
      U_NODES_QUAD(52,1) = U(73,1)   ! NODE 18
      U_NODES_QUAD(53,1) = U(74,1)
      U_NODES_QUAD(54,1) = U(75,1)
      
      U_NODES_QUAD(55,1) = U(76,1)   ! NODE 19
      U_NODES_QUAD(56,1) = U(77,1)
      U_NODES_QUAD(57,1) = U(78,1)
      THETA_NODES_QUAD(22,1) = U(79,1)
      THETA_NODES_QUAD(23,1) = U(80,1)
      THETA_NODES_QUAD(24,1) = U(81,1)
      
      U_NODES_QUAD(58,1) = U(82,1)   ! NODE 20
      U_NODES_QUAD(59,1) = U(83,1)
      U_NODES_QUAD(60,1) = U(84,1)
      
      RETURN 
     	END             
C          
C%*************************************************************************
C% Compute BN AND BG MATRIX FROM DERIVATIVE OF SHAPE FUNCTIONS AND DEFORMATION GRADIENTS
C%*************************************************************************
C%%
        SUBROUTINE MAKER_BNBG_HEX20(SHPD,F,BN,BG)
		
        IMPLICIT REAL*8 (A-H,O-Z)
        
        INTEGER I, J, K, COL1, COL2, COL3
          
        REAL*8 SHPD(3,20), F(3,3), BN(6,60), BG(9,60)
        
        PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )
     
     		    DO I = 1, 6
    	        DO J = 1, 60
	    	    BN(I,J) = ZERO
              END DO
              END DO			
C
              DO I = 1, 9
    	        DO J = 1, 60
                  BG(I,J) = ZERO
              END DO
              END DO
C
			DO I = 1, 20
				COL1 = (I-1)*3+1
				COL2 = (I-1)*3+2
				COL3 = (I-1)*3+3
				
				BN(1,COL1) = SHPD(1,I)*F(1,1)
				BN(2,COL1) = SHPD(2,I)*F(1,2)
				BN(3,COL1) = SHPD(3,I)*F(1,3)
				BN(4,COL1) = SHPD(1,I)*F(1,2)+SHPD(2,I)*F(1,1)
				BN(5,COL1) = SHPD(2,I)*F(1,3)+SHPD(3,I)*F(1,2)
				BN(6,COL1) = SHPD(1,I)*F(1,3)+SHPD(3,I)*F(1,1)
				
				BN(1,COL2) = SHPD(1,I)*F(2,1)
				BN(2,COL2) = SHPD(2,I)*F(2,2)
				BN(3,COL2) = SHPD(3,I)*F(2,3)
				BN(4,COL2) = SHPD(1,I)*F(2,2)+SHPD(2,I)*F(2,1)
				BN(5,COL2) = SHPD(2,I)*F(2,3)+SHPD(3,I)*F(2,2)
				BN(6,COL2) = SHPD(1,I)*F(2,3)+SHPD(3,I)*F(2,1)
				
				BN(1,COL3) = SHPD(1,I)*F(3,1)
				BN(2,COL3) = SHPD(2,I)*F(3,2)
				BN(3,COL3) = SHPD(3,I)*F(3,3)
				BN(4,COL3) = SHPD(1,I)*F(3,2)+SHPD(2,I)*F(3,1)
				BN(5,COL3) = SHPD(2,I)*F(3,3)+SHPD(3,I)*F(3,2)
				BN(6,COL3) = SHPD(1,I)*F(3,3)+SHPD(3,I)*F(3,1)
C----------------------------------
				BG(1,COL1) = SHPD(1,I)	
				BG(2,COL1) = SHPD(2,I)	
				BG(3,COL1) = SHPD(3,I)		
				
				BG(4,COL2) = SHPD(1,I)	
				BG(5,COL2) = SHPD(2,I)	
				BG(6,COL2) = SHPD(3,I)	
				
				BG(7,COL3) = SHPD(1,I)	
				BG(8,COL3) = SHPD(2,I)	
				BG(9,COL3) = SHPD(3,I)	
				
			END DO     	
     	
     	RETURN 
     	END
C          
C%*************************************************************************
C% Compute shape function, derivatives, and determinant of 20 NODE hexahedron element
C%*************************************************************************
C%%
      SUBROUTINE SHAPEL_HEX20(XI,ELXYZ_QUAD,DETGJ_QUAD,GJ_QUAD,
     +             SHPD_QUAD)

        !IMPLICIT REAL*8 (A-H,O-Z)
        
        INTEGER I, J, K
          
        REAL*8 XNODE(3,20), QUAR,SF(20,1),DSF(3,20),SHPD_QUAD(3,20),
     +    XI(1,3),ELXYZ_QUAD(20,3), XP, YP, ZP, XI0(1,3),GJ_QUAD(3,3),
     +    GJINV(3,3), DETGJ_QUAD
        
        PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0,
     +          ONE_FOURTH=1.D0/4.D0  )
        
        !PRINT *, "ELXYZ_QUAD = "
        !  CALL PRTMAT(ELXYZ_QUAD,20,3)
        
        XNODE(1,1) = -ONE
        XNODE(1,2) = ZEO
        XNODE(1,3) = ONE
        XNODE(1,4) = ONE
        XNODE(1,5) = ONE
        XNODE(1,6) = ZERO
        XNODE(1,7) = -ONE
        XNODE(1,8) = -ONE        
        XNODE(1,9) = -ONE
        XNODE(1,10) = ONE
        XNODE(1,11) = ONE
        XNODE(1,12) = -ONE        
        XNODE(1,13) = -ONE
        XNODE(1,14) = ZERO
        XNODE(1,15) = ONE
        XNODE(1,16) = ONE
        XNODE(1,17) = ONE
        XNODE(1,18) = ZERO
        XNODE(1,19) = -ONE
        XNODE(1,20) = -ONE
        
        XNODE(2,1) = -ONE
        XNODE(2,2) = -ONE
        XNODE(2,3) = -ONE
        XNODE(2,4) = ZERO
        XNODE(2,5) = ONE
        XNODE(2,6) = ONE
        XNODE(2,7) = ONE
        XNODE(2,8) = ZERO
        XNODE(2,9) = -ONE
        XNODE(2,10) = -ONE
        XNODE(2,11) = ONE
        XNODE(2,12) = ONE       
        XNODE(2,13) = -ONE
        XNODE(2,14) = -ONE
        XNODE(2,15) = -ONE
        XNODE(2,16) = ZERO
        XNODE(2,17) = ONE
        XNODE(2,18) = ONE
        XNODE(2,19) = ONE
        XNODE(2,20) = ZERO
        
        XNODE(3,1) = -ONE
        XNODE(3,2) = -ONE
        XNODE(3,3) = -ONE
        XNODE(3,4) = -ONE
        XNODE(3,5) = -ONE
        XNODE(3,6) = -ONE
        XNODE(3,7) = -ONE
        XNODE(3,8) = -ONE
        XNODE(3,9) = ZERO
        XNODE(3,10) = ZERO
        XNODE(3,11) = ZERO
        XNODE(3,12) = ZERO
        XNODE(3,13) = ONE
        XNODE(3,14) = ONE
        XNODE(3,15) = ONE
        XNODE(3,16) = ONE
        XNODE(3,17) = ONE
        XNODE(3,18) = ONE
        XNODE(3,19) = ONE
        XNODE(3,20) = ONE
        
        QUAR = 0.125D0

        DO I = 1, 20
	        SF(I,1) = ZERO
	        DSF(1,I) = ZERO
	        DSF(2,I) = ZERO
	        DSF(3,I) = ZERO
        END DO
		
          DO I = 1, 20
			XP = XNODE(1,I)
    		    YP = XNODE(2,I)
    		    ZP = XNODE(3,I)
C
    		    XI0(1,1) = XI(1,1)*XP   	
    		    XI0(1,2) = XI(1,2)*YP		
    		    XI0(1,3) = XI(1,3)*ZP
C		
C 		! ---> THIS IS SHAPE FUNCTION AT EACH GAUSS POINTS
C
			SF(I,1) = QUAR*(ONE+XI0(1,1))*(ONE+XI0(1,2))*(ONE+XI0(1,3))*(XI0(1,1)
     +         + XI0(1,2)+XI0(1,3)-TWO)
			
			IF ( (I.EQ.2) .OR. (I.EQ.6) .OR. (I.EQ.14) .OR. (I.EQ.18) ) THEN			
				SF(I,1) = ONE_FOURTH*(ONE-XI(1,1)*XI(1,1))*
     +                 (ONE+XI0(1,2))*(ONE+XI0(1,3))				
			ELSE IF ( (I.EQ.4) .OR. (I.EQ.8) .OR. (I.EQ.16) .OR. (I.EQ.20) ) THEN
				SF(I,1) = ONE_FOURTH*(ONE-XI(1,2)*XI(1,2))*
     +                 (ONE+XI0(1,1))*(ONE+XI0(1,3))			
			ELSE IF ((I.EQ.9) .OR. (I.EQ.10) .OR. (I.EQ.11) .OR. (I.EQ.12)) THEN
				SF(I,1) = ONE_FOURTH*(ONE-XI(1,3)*XI(1,3))*
     +                 (ONE+XI0(1,1))*(ONE+XI0(1,2))				
			END IF

C		
C 		! ---> THESE ARE THE DERIVATIVE OF SHAPE FUNCTION W.R.T NETURAL COORDS AT EACH GAUSS POINTS
C			
    		    DSF(1,I) = QUAR*XP*(ONE+XI0(1,2))*(ONE+XI0(1,3))*(XI0(1,1)+
     +             XI0(1,2)+XI0(1,3)-TWO) + 
     +             QUAR*(ONE+XI0(1,1))*(ONE+XI0(1,2))*(ONE+XI0(1,3))*XP
    		    DSF(2,I) = QUAR*YP*(ONE+XI0(1,1))*(ONE+XI0(1,3))*(XI0(1,1)+
     +             XI0(1,2)+XI0(1,3)-TWO) + 
     +             QUAR*(ONE+XI0(1,1))*(ONE+XI0(1,2))*(ONE+XI0(1,3))*YP
    		    DSF(3,I) = QUAR*ZP*(ONE+XI0(1,1))*(ONE+XI0(1,2))*(XI0(1,1)+
     +             XI0(1,2)+XI0(1,3)-TWO) + 
     +             QUAR*(ONE+XI0(1,1))*(ONE+XI0(1,2))*(ONE+XI0(1,3))*ZP	
    		
    		    IF ((I.EQ.2) .OR. (I.EQ.6) .OR. (I.EQ.14) .OR. (I.EQ.18)) THEN
				DSF(1,I) = -ONE_HALF*XI(1,1)*(ONE+XI0(1,2))*(ONE+XI0(1,3))
				DSF(2,I) = ONE_FOURTH*YP*(ONE-XI(1,1)*XI(1,1))*(ONE+XI0(1,3))
				DSF(3,I) = ONE_FOURTH*ZP*(ONE-XI(1,1)*XI(1,1))*(ONE+XI0(1,2))
			ELSE IF ((I.EQ.4) .OR. (I.EQ.8) .OR. (I.EQ.16) .OR. (I.EQ.20)) THEN
				DSF(1,I) = ONE_FOURTH*XP*(ONE-XI(1,2)*XI(1,2))*(ONE+XI0(1,3))
				DSF(2,I) = -ONE_HALF*XI(1,2)*(ONE+XI0(1,1))*(ONE+XI0(1,3))
				DSF(3,I) = ONE_FOURTH*ZP*(ONE-XI(1,2)*XI(1,2))*(ONE+XI0(1,1))			
			ELSE IF ((I.EQ.9) .OR. (I.EQ.10) .OR. (I.EQ.11) .OR. (I.EQ.12)) THEN
				DSF(1,I) = ONE_FOURTH*XP*(ONE-XI(1,3)*XI(1,3))*(ONE+XI0(1,2))
				DSF(2,I) = ONE_FOURTH*YP*(ONE-XI(1,3)*XI(1,3))*(ONE+XI0(1,1))
				DSF(3,I) = -ONE_HALF*XI(1,3)*(ONE+XI0(1,1))*(ONE+XI0(1,2))				
			END IF

          END DO
        
C      CALCULATE GLOBAL JACOUBIAN :=  GJ_QUAD(3,3) = DSF(3,20)*ELXY(20,3) 
		DO I = 1,3
		DO J = 1,3
			GJ_QUAD(I,J) = ZERO
			DO K = 1,20
				GJ_QUAD(I,J) = GJ_QUAD(I,J) + DSF(I,K)*ELXYZ_QUAD(K,J)
			END DO
		END DO
		END DO
		
C		CALCULATE DETERMINATE GJ_QUAD AND JACOUBIAN INVERSE :=  DETGJ_QUAD = det(GJ) ;  GJINV=inv(GJ)
		CALL MATINV(GJ_QUAD,GJINV,DETGJ_QUAD)

C      CALCULATE GLOBAL DERIVATION OF SHAPE FUNCTIONS :=   SHPD_QUAD(3,20) = GJINV(3,3)*DSF(3,20)
		DO I = 1,3
		DO J = 1,20
			SHPD_QUAD(I,J) = ZERO
			DO K = 1,3
				SHPD_QUAD(I,J) = SHPD_QUAD(I,J) + GJINV(I,K)*DSF(K,J)
			END DO
		END DO
		END DO
C	
		RETURN
      	END  
C          
C%*************************************************************************
C% CONVERT AMATRX_OLD_QUAD,RHS_OLD_QUAD TO AMATRX,RHS
C%*************************************************************************
C%%              
      SUBROUTINE AMATRX_RHS_CONVERSION_QUAD(AMATRX_OLD,RHS_OLD,
     +         AMATRX_NEW,RHS_NEW)
            
      IMPLICIT REAL*8 (A-H,O-Z)
        
      INTEGER I, J, K
          
      REAL*8 AMATRX_OLD(84,84), AMATRX_NEW(84,84), RHS_OLD(84,1), 
     +        AMATRX_PART1(84,60), AMATRX_PART2(84,24), 
     +        AMATRX_SEMI_OLD(84,84), RHS_NEW(84,1),
     +        AMATRX_PART3(60,84), AMATRX_PART4(24,84),
     +        RHS_PART1(60,1),RHS_PART2(24,1) 
      
      PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )

C          
C         /////////////////////////////////////////// 
C
          
          DO I = 1, 84
          DO J = 1, 60
              AMATRX_PART1(I,J) = ZERO
              AMATRX_PART3(J,I) = ZERO
              RHS_PART1(J,1)=ZERO
          END DO
          END DO
          
          DO I = 1, 84
          DO J = 1, 24
              AMATRX_PART2(I,J) = ZERO
              AMATRX_PART4(J,I) = ZERO
              RHS_PART2(J,1)=ZERO
          END DO
          END DO
          
          
          DO I = 1, 84
          DO J = 1, 84
              AMATRX_SEMI_OLD(I,J) = ZERO
              AMATRX_NEW(I,J) = ZERO
          END DO
          END DO
      
          DO I = 1, 84
              RHS_NEW(I,1)=ZERO
          END DO
C          
C         /////////////////////////////////////////// 
C
          DO I = 1, 84
          DO J = 1, 60
              AMATRX_PART1(I,J) = AMATRX_OLD(I,J)
          END DO
          END DO
          
          DO I = 1, 84
          DO J = 1, 24
              AMATRX_PART2(I,J) = AMATRX_OLD(I,J+60)
          END DO
          END DO
          
          DO I = 1, 84
              ! NODE 1
              AMATRX_SEMI_OLD(I,1) = AMATRX_PART1(I,1)
              AMATRX_SEMI_OLD(I,2) = AMATRX_PART1(I,2)
              AMATRX_SEMI_OLD(I,3) = AMATRX_PART1(I,3)
              AMATRX_SEMI_OLD(I,4) = AMATRX_PART2(I,1)
              AMATRX_SEMI_OLD(I,5) = AMATRX_PART2(I,2)
              AMATRX_SEMI_OLD(I,6) = AMATRX_PART2(I,3)
              ! NODE 2
              AMATRX_SEMI_OLD(I,7) = AMATRX_PART1(I,4)
              AMATRX_SEMI_OLD(I,8) = AMATRX_PART1(I,5)
              AMATRX_SEMI_OLD(I,9) = AMATRX_PART1(I,6)
              ! NODE 3
              AMATRX_SEMI_OLD(I,10) = AMATRX_PART1(I,7)
              AMATRX_SEMI_OLD(I,11) = AMATRX_PART1(I,8)
              AMATRX_SEMI_OLD(I,12) = AMATRX_PART1(I,9)
              AMATRX_SEMI_OLD(I,13) = AMATRX_PART2(I,4)
              AMATRX_SEMI_OLD(I,14) = AMATRX_PART2(I,5)
              AMATRX_SEMI_OLD(I,15) = AMATRX_PART2(I,6)
              ! NODE 4
              AMATRX_SEMI_OLD(I,16) = AMATRX_PART1(I,10)
              AMATRX_SEMI_OLD(I,17) = AMATRX_PART1(I,11)
              AMATRX_SEMI_OLD(I,18) = AMATRX_PART1(I,12)
              ! NODE 5
              AMATRX_SEMI_OLD(I,19) = AMATRX_PART1(I,13)
              AMATRX_SEMI_OLD(I,20) = AMATRX_PART1(I,14)
              AMATRX_SEMI_OLD(I,21) = AMATRX_PART1(I,15)
              AMATRX_SEMI_OLD(I,22) = AMATRX_PART2(I,7)
              AMATRX_SEMI_OLD(I,23) = AMATRX_PART2(I,8)
              AMATRX_SEMI_OLD(I,24) = AMATRX_PART2(I,9)
              ! NODE 6
              AMATRX_SEMI_OLD(I,25) = AMATRX_PART1(I,16)
              AMATRX_SEMI_OLD(I,26) = AMATRX_PART1(I,17)
              AMATRX_SEMI_OLD(I,27) = AMATRX_PART1(I,18)
              ! NODE 7
              AMATRX_SEMI_OLD(I,28) = AMATRX_PART1(I,19)
              AMATRX_SEMI_OLD(I,29) = AMATRX_PART1(I,20)
              AMATRX_SEMI_OLD(I,30) = AMATRX_PART1(I,21)
              AMATRX_SEMI_OLD(I,31) = AMATRX_PART2(I,10)
              AMATRX_SEMI_OLD(I,32) = AMATRX_PART2(I,11)
              AMATRX_SEMI_OLD(I,33) = AMATRX_PART2(I,12)
              ! NODE 8
              AMATRX_SEMI_OLD(I,34) = AMATRX_PART1(I,22)
              AMATRX_SEMI_OLD(I,35) = AMATRX_PART1(I,23)
              AMATRX_SEMI_OLD(I,36) = AMATRX_PART1(I,24)
              !!!! MIDDLE NODES !!! 9 , 10 , 11 , 12
              ! NODE 9
              AMATRX_SEMI_OLD(I,37) = AMATRX_PART1(I,25)
              AMATRX_SEMI_OLD(I,38) = AMATRX_PART1(I,26)
              AMATRX_SEMI_OLD(I,39) = AMATRX_PART1(I,27)
              ! NODE 10
              AMATRX_SEMI_OLD(I,40) = AMATRX_PART1(I,28)
              AMATRX_SEMI_OLD(I,41) = AMATRX_PART1(I,29)
              AMATRX_SEMI_OLD(I,42) = AMATRX_PART1(I,30)
              ! NODE 11
              AMATRX_SEMI_OLD(I,43) = AMATRX_PART1(I,31)
              AMATRX_SEMI_OLD(I,44) = AMATRX_PART1(I,32)
              AMATRX_SEMI_OLD(I,45) = AMATRX_PART1(I,33)
              ! NODE 12
              AMATRX_SEMI_OLD(I,46) = AMATRX_PART1(I,34)
              AMATRX_SEMI_OLD(I,47) = AMATRX_PART1(I,35)
              AMATRX_SEMI_OLD(I,48) = AMATRX_PART1(I,36)
              !!!!!!!!!!!!!!!
              ! NODE 13
              AMATRX_SEMI_OLD(I,49) = AMATRX_PART1(I,37)
              AMATRX_SEMI_OLD(I,50) = AMATRX_PART1(I,38)
              AMATRX_SEMI_OLD(I,51) = AMATRX_PART1(I,39)
              AMATRX_SEMI_OLD(I,52) = AMATRX_PART2(I,13)
              AMATRX_SEMI_OLD(I,53) = AMATRX_PART2(I,14)
              AMATRX_SEMI_OLD(I,54) = AMATRX_PART2(I,15)
              ! NODE 14
              AMATRX_SEMI_OLD(I,55) = AMATRX_PART1(I,40)
              AMATRX_SEMI_OLD(I,56) = AMATRX_PART1(I,41)
              AMATRX_SEMI_OLD(I,57) = AMATRX_PART1(I,42)
              ! NODE 15
              AMATRX_SEMI_OLD(I,58) = AMATRX_PART1(I,43)
              AMATRX_SEMI_OLD(I,59) = AMATRX_PART1(I,44)
              AMATRX_SEMI_OLD(I,60) = AMATRX_PART1(I,45)
              AMATRX_SEMI_OLD(I,61) = AMATRX_PART2(I,16)
              AMATRX_SEMI_OLD(I,62) = AMATRX_PART2(I,17)
              AMATRX_SEMI_OLD(I,63) = AMATRX_PART2(I,18)
              ! NODE 16
              AMATRX_SEMI_OLD(I,64) = AMATRX_PART1(I,46)
              AMATRX_SEMI_OLD(I,65) = AMATRX_PART1(I,47)
              AMATRX_SEMI_OLD(I,66) = AMATRX_PART1(I,48)
              ! NODE 17
              AMATRX_SEMI_OLD(I,67) = AMATRX_PART1(I,49)
              AMATRX_SEMI_OLD(I,68) = AMATRX_PART1(I,50)
              AMATRX_SEMI_OLD(I,69) = AMATRX_PART1(I,51)
              AMATRX_SEMI_OLD(I,70) = AMATRX_PART2(I,19)
              AMATRX_SEMI_OLD(I,71) = AMATRX_PART2(I,20)
              AMATRX_SEMI_OLD(I,72) = AMATRX_PART2(I,21)
              ! NODE 18
              AMATRX_SEMI_OLD(I,73) = AMATRX_PART1(I,52)
              AMATRX_SEMI_OLD(I,74) = AMATRX_PART1(I,53)
              AMATRX_SEMI_OLD(I,75) = AMATRX_PART1(I,54)
              ! NODE 19
              AMATRX_SEMI_OLD(I,76) = AMATRX_PART1(I,55)
              AMATRX_SEMI_OLD(I,77) = AMATRX_PART1(I,56)
              AMATRX_SEMI_OLD(I,78) = AMATRX_PART1(I,57)
              AMATRX_SEMI_OLD(I,79) = AMATRX_PART2(I,22)
              AMATRX_SEMI_OLD(I,80) = AMATRX_PART2(I,23)
              AMATRX_SEMI_OLD(I,81) = AMATRX_PART2(I,24)
              ! NODE 20
              AMATRX_SEMI_OLD(I,82) = AMATRX_PART1(I,58)
              AMATRX_SEMI_OLD(I,83) = AMATRX_PART1(I,59)
              AMATRX_SEMI_OLD(I,84) = AMATRX_PART1(I,60)

          END DO
C SECOND PHASE OF CONVERSION          
          DO J = 1, 84
          DO I = 1, 60
              AMATRX_PART3(I,J) = AMATRX_SEMI_OLD(I,J)
          END DO
          END DO
          
          DO J = 1, 84
          DO I = 1, 24
              AMATRX_PART4(I,J) = AMATRX_SEMI_OLD(I+60,J)
          END DO
          END DO
          
          DO J = 1, 84

              ! NODE 1
              AMATRX_NEW(1,J) = AMATRX_PART3(1,J)
              AMATRX_NEW(2,J) = AMATRX_PART3(2,J)
              AMATRX_NEW(3,J) = AMATRX_PART3(3,J)
              AMATRX_NEW(4,J) = AMATRX_PART4(1,J)
              AMATRX_NEW(5,J) = AMATRX_PART4(2,J)
              AMATRX_NEW(6,J) = AMATRX_PART4(3,J)
              ! NODE 2
              AMATRX_NEW(7,J) = AMATRX_PART3(4,J)
              AMATRX_NEW(8,J) = AMATRX_PART3(5,J)
              AMATRX_NEW(9,J) = AMATRX_PART3(6,J)
              ! NODE 3
              AMATRX_NEW(10,J) = AMATRX_PART3(7,J)
              AMATRX_NEW(11,J) = AMATRX_PART3(8,J)
              AMATRX_NEW(12,J) = AMATRX_PART3(9,J)
              AMATRX_NEW(13,J) = AMATRX_PART4(4,J)
              AMATRX_NEW(14,J) = AMATRX_PART4(5,J)
              AMATRX_NEW(15,J) = AMATRX_PART4(6,J)
              ! NODE 4
              AMATRX_NEW(16,J) = AMATRX_PART3(10,J)
              AMATRX_NEW(17,J) = AMATRX_PART3(11,J)
              AMATRX_NEW(18,J) = AMATRX_PART3(12,J)
              ! NODE 5
              AMATRX_NEW(19,J) = AMATRX_PART3(13,J)
              AMATRX_NEW(20,J) = AMATRX_PART3(14,J)
              AMATRX_NEW(21,J) = AMATRX_PART3(15,J)
              AMATRX_NEW(22,J) = AMATRX_PART4(7,J)
              AMATRX_NEW(23,J) = AMATRX_PART4(8,J)
              AMATRX_NEW(24,J) = AMATRX_PART4(9,J)
              ! NODE 6
              AMATRX_NEW(25,J) = AMATRX_PART3(16,J)
              AMATRX_NEW(26,J) = AMATRX_PART3(17,J)
              AMATRX_NEW(27,J) = AMATRX_PART3(18,J)
              ! NODE 7
              AMATRX_NEW(28,J) = AMATRX_PART3(19,J)
              AMATRX_NEW(29,J) = AMATRX_PART3(20,J)
              AMATRX_NEW(30,J) = AMATRX_PART3(21,J)
              AMATRX_NEW(31,J) = AMATRX_PART4(10,J)
              AMATRX_NEW(32,J) = AMATRX_PART4(11,J)
              AMATRX_NEW(33,J) = AMATRX_PART4(12,J)
              ! NODE 8
              AMATRX_NEW(34,J) = AMATRX_PART3(22,J)
              AMATRX_NEW(35,J) = AMATRX_PART3(23,J)
              AMATRX_NEW(36,J) = AMATRX_PART3(24,J)
              !!!! MIDDLE NODES !!! 9 , 10 , 11 , 12
              ! NODE 9
              AMATRX_NEW(37,J) = AMATRX_PART3(25,J)
              AMATRX_NEW(38,J) = AMATRX_PART3(26,J)
              AMATRX_NEW(39,J) = AMATRX_PART3(27,J)
              ! NODE 10
              AMATRX_NEW(40,J) = AMATRX_PART3(28,J)
              AMATRX_NEW(41,J) = AMATRX_PART3(29,J)
              AMATRX_NEW(42,J) = AMATRX_PART3(30,J)
              ! NODE 11
              AMATRX_NEW(43,J) = AMATRX_PART3(31,J)
              AMATRX_NEW(44,J) = AMATRX_PART3(32,J)
              AMATRX_NEW(45,J) = AMATRX_PART3(33,J)
              ! NODE 12
              AMATRX_NEW(46,J) = AMATRX_PART3(34,J)
              AMATRX_NEW(47,J) = AMATRX_PART3(35,J)
              AMATRX_NEW(48,J) = AMATRX_PART3(36,J)
              !!!!!!!!!!!!!!!
              ! NODE 13
              AMATRX_NEW(49,J) = AMATRX_PART3(37,J)
              AMATRX_NEW(50,J) = AMATRX_PART3(38,J)
              AMATRX_NEW(51,J) = AMATRX_PART3(39,J)
              AMATRX_NEW(52,J) = AMATRX_PART4(13,J)
              AMATRX_NEW(53,J) = AMATRX_PART4(14,J)
              AMATRX_NEW(54,J) = AMATRX_PART4(15,J)
              ! NODE 14
              AMATRX_NEW(55,J) = AMATRX_PART3(40,J)
              AMATRX_NEW(56,J) = AMATRX_PART3(41,J)
              AMATRX_NEW(57,J) = AMATRX_PART3(42,J)
              ! NODE 15
              AMATRX_NEW(58,J) = AMATRX_PART3(43,J)
              AMATRX_NEW(59,J) = AMATRX_PART3(44,J)
              AMATRX_NEW(60,J) = AMATRX_PART3(45,J)
              AMATRX_NEW(61,J) = AMATRX_PART4(16,J)
              AMATRX_NEW(62,J) = AMATRX_PART4(17,J)
              AMATRX_NEW(63,J) = AMATRX_PART4(18,J)
              ! NODE 16
              AMATRX_NEW(64,J) = AMATRX_PART3(46,J)
              AMATRX_NEW(65,J) = AMATRX_PART3(47,J)
              AMATRX_NEW(66,J) = AMATRX_PART3(48,J)
              ! NODE 17
              AMATRX_NEW(67,J) = AMATRX_PART3(49,J)
              AMATRX_NEW(68,J) = AMATRX_PART3(50,J)
              AMATRX_NEW(69,J) = AMATRX_PART3(51,J)
              AMATRX_NEW(70,J) = AMATRX_PART4(19,J)
              AMATRX_NEW(71,J) = AMATRX_PART4(20,J)
              AMATRX_NEW(72,J) = AMATRX_PART4(21,J)
              ! NODE 18
              AMATRX_NEW(73,J) = AMATRX_PART3(52,J)
              AMATRX_NEW(74,J) = AMATRX_PART3(53,J)
              AMATRX_NEW(75,J) = AMATRX_PART3(54,J)
              ! NODE 19
              AMATRX_NEW(76,J) = AMATRX_PART3(55,J)
              AMATRX_NEW(77,J) = AMATRX_PART3(56,J)
              AMATRX_NEW(78,J) = AMATRX_PART3(57,J)
              AMATRX_NEW(79,J) = AMATRX_PART4(22,J)
              AMATRX_NEW(80,J) = AMATRX_PART4(23,J)
              AMATRX_NEW(81,J) = AMATRX_PART4(24,J)
              ! NODE 20
              AMATRX_NEW(82,J) = AMATRX_PART3(58,J)
              AMATRX_NEW(83,J) = AMATRX_PART3(59,J)
              AMATRX_NEW(84,J) = AMATRX_PART3(60,J)

          END DO
C
C NOW IT IS TIME FOR RHS          

          DO I = 1, 60
              RHS_PART1(I,1) = RHS_OLD(I,1)
          END DO
          DO I = 1, 24
              RHS_PART2(I,1) = RHS_OLD(I+60,1)
          END DO
          
          DO I = 1, 84
              ! NODE 1
              RHS_NEW(1,1) = RHS_PART1(1,1)
              RHS_NEW(2,1) = RHS_PART1(2,1)
              RHS_NEW(3,1) = RHS_PART1(3,1)
              RHS_NEW(4,1) = RHS_PART2(1,1)
              RHS_NEW(5,1) = RHS_PART2(2,1)
              RHS_NEW(6,1) = RHS_PART2(3,1)
              ! NODE 2
              RHS_NEW(7,1) = RHS_PART1(4,1)
              RHS_NEW(8,1) = RHS_PART1(5,1)
              RHS_NEW(9,1) = RHS_PART1(6,1)
              ! NODE 3
              RHS_NEW(10,1) = RHS_PART1(7,1)
              RHS_NEW(11,1) = RHS_PART1(8,1)
              RHS_NEW(12,1) = RHS_PART1(9,1)
              RHS_NEW(13,1) = RHS_PART2(4,1)
              RHS_NEW(14,1) = RHS_PART2(5,1)
              RHS_NEW(15,1) = RHS_PART2(6,1)
              ! NODE 4
              RHS_NEW(16,1) = RHS_PART1(10,1)
              RHS_NEW(17,1) = RHS_PART1(11,1)
              RHS_NEW(18,1) = RHS_PART1(12,1)
              ! NODE 5
              RHS_NEW(19,1) = RHS_PART1(13,1)
              RHS_NEW(20,1) = RHS_PART1(14,1)
              RHS_NEW(21,1) = RHS_PART1(15,1)
              RHS_NEW(22,1) = RHS_PART2(7,1)
              RHS_NEW(23,1) = RHS_PART2(8,1)
              RHS_NEW(24,1) = RHS_PART2(9,1)
              ! NODE 6
              RHS_NEW(25,1) = RHS_PART1(16,1)
              RHS_NEW(26,1) = RHS_PART1(17,1)
              RHS_NEW(27,1) = RHS_PART1(18,1)
              ! NODE 7
              RHS_NEW(28,1) = RHS_PART1(19,1)
              RHS_NEW(29,1) = RHS_PART1(20,1)
              RHS_NEW(30,1) = RHS_PART1(21,1)
              RHS_NEW(31,1) = RHS_PART2(10,1)
              RHS_NEW(32,1) = RHS_PART2(11,1)
              RHS_NEW(33,1) = RHS_PART2(12,1)
              ! NODE 8
              RHS_NEW(34,1) = RHS_PART1(22,1)
              RHS_NEW(35,1) = RHS_PART1(23,1)
              RHS_NEW(36,1) = RHS_PART1(24,1)
              !!!! MIDDLE NODES !!! 9 , 10 , 11 , 12
              ! NODE 9
              RHS_NEW(37,1) = RHS_PART1(25,1)
              RHS_NEW(38,1) = RHS_PART1(26,1)
              RHS_NEW(39,1) = RHS_PART1(27,1)
              ! NODE 10
              RHS_NEW(40,1) = RHS_PART1(28,1)
              RHS_NEW(41,1) = RHS_PART1(29,1)
              RHS_NEW(42,1) = RHS_PART1(30,1)
              ! NODE 11
              RHS_NEW(43,1) = RHS_PART1(31,1)
              RHS_NEW(44,1) = RHS_PART1(32,1)
              RHS_NEW(45,1) = RHS_PART1(33,1)
              ! NODE 12
              RHS_NEW(46,1) = RHS_PART1(34,1)
              RHS_NEW(47,1) = RHS_PART1(35,1)
              RHS_NEW(48,1) = RHS_PART1(36,1)
              !!!!!!!!!!!!!!!
              ! NODE 13
              RHS_NEW(49,1) = RHS_PART1(37,1)
              RHS_NEW(50,1) = RHS_PART1(38,1)
              RHS_NEW(51,1) = RHS_PART1(39,1)
              RHS_NEW(52,1) = RHS_PART2(13,1)
              RHS_NEW(53,1) = RHS_PART2(14,1)
              RHS_NEW(54,1) = RHS_PART2(15,1)
              ! NODE 14
              RHS_NEW(55,1) = RHS_PART1(40,1)
              RHS_NEW(56,1) = RHS_PART1(41,1)
              RHS_NEW(57,1) = RHS_PART1(42,1)
              ! NODE 15
              RHS_NEW(58,1) = RHS_PART1(43,1)
              RHS_NEW(59,1) = RHS_PART1(44,1)
              RHS_NEW(60,1) = RHS_PART1(45,1)
              RHS_NEW(61,1) = RHS_PART2(16,1)
              RHS_NEW(62,1) = RHS_PART2(17,1)
              RHS_NEW(63,1) = RHS_PART2(18,1)
              ! NODE 16
              RHS_NEW(64,1) = RHS_PART1(46,1)
              RHS_NEW(65,1) = RHS_PART1(47,1)
              RHS_NEW(66,1) = RHS_PART1(48,1)
              ! NODE 17
              RHS_NEW(67,1) = RHS_PART1(49,1)
              RHS_NEW(68,1) = RHS_PART1(50,1)
              RHS_NEW(69,1) = RHS_PART1(51,1)
              RHS_NEW(70,1) = RHS_PART2(19,1)
              RHS_NEW(71,1) = RHS_PART2(20,1)
              RHS_NEW(72,1) = RHS_PART2(21,1)
              ! NODE 18
              RHS_NEW(73,1) = RHS_PART1(52,1)
              RHS_NEW(74,1) = RHS_PART1(53,1)
              RHS_NEW(75,1) = RHS_PART1(54,1)
              ! NODE 19
              RHS_NEW(76,1) = RHS_PART1(55,1)
              RHS_NEW(77,1) = RHS_PART1(56,1)
              RHS_NEW(78,1) = RHS_PART1(57,1)
              RHS_NEW(79,1) = RHS_PART2(22,1)
              RHS_NEW(80,1) = RHS_PART2(23,1)
              RHS_NEW(81,1) = RHS_PART2(24,1)
              ! NODE 20
              RHS_NEW(82,1) = RHS_PART1(58,1)
              RHS_NEW(83,1) = RHS_PART1(59,1)
              RHS_NEW(84,1) = RHS_PART1(60,1)
          
          END DO
C
C
          
      RETURN 
     	END
C          
C%*************************************************************************
C% CONVERT AMATRX_OLD,RHS_OLD TO AMATRX,RHS
C%*************************************************************************
C%%              
      SUBROUTINE AMATRX_RHS_CONVERSION(AMATRX_OLD,RHS_OLD,AMATRX_NEW,
     +         RHS_NEW)
            
      IMPLICIT REAL*8 (A-H,O-Z)
        
      INTEGER I, J, K
          
      REAL*8 AMATRX_OLD(48,48), AMATRX_NEW(48,48), RHS_OLD(48,1), 
     +        AMATRX_PART1(48,24), AMATRX_PART2(48,24), 
     +        AMATRX_SEMI_OLD(48,48), RHS_NEW(48,1),
     +        AMATRX_PART3(24,48), AMATRX_PART4(24,48),
     +        RHS_PART1(24,1),RHS_PART2(24,1) 
      
      PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )
      
          !print *, 'AMATRX_OLD:'
          !CALL PRTMAT(AMATRX_OLD,48,48)
          
          DO I = 1, 48
          DO J = 1, 24
              AMATRX_PART1(I,J) = ZERO
              AMATRX_PART2(I,J) = ZERO
              AMATRX_PART3(J,I) = ZERO
              AMATRX_PART4(J,I) = ZERO
              RHS_PART1(J,1)=ZERO
              RHS_PART2(J,1)=ZERO
          END DO
          END DO
          
          
          DO I = 1, 48
          DO J = 1, 48
              AMATRX_SEMI_OLD(I,J) = ZERO
              AMATRX_NEW(I,J) = ZERO
          END DO
          END DO
      
          
          DO I = 1, 24
              RHS_PART1(I,1)=ZERO
              RHS_PART2(I,1)=ZERO
          END DO
          
          DO I = 1, 48
              RHS_NEW(I,1)=ZERO
          END DO
C          
C         ////////////////////// 
C
          DO I = 1, 48
          DO J = 1, 24
              AMATRX_PART1(I,J) = AMATRX_OLD(I,J)
              AMATRX_PART2(I,J) = AMATRX_OLD(I,J+24)
          END DO
          END DO
          
          DO I = 1, 48
          DO J = 1, 8
              AMATRX_SEMI_OLD(I,6*(J-1)+1) = AMATRX_PART1(I,3*(J-1)+1)
              AMATRX_SEMI_OLD(I,6*(J-1)+2) = AMATRX_PART1(I,3*(J-1)+2)
              AMATRX_SEMI_OLD(I,6*(J-1)+3) = AMATRX_PART1(I,3*(J-1)+3)
              AMATRX_SEMI_OLD(I,6*(J-1)+4) = AMATRX_PART2(I,3*(J-1)+1)
              AMATRX_SEMI_OLD(I,6*(J-1)+5) = AMATRX_PART2(I,3*(J-1)+2)
              AMATRX_SEMI_OLD(I,6*(J-1)+6) = AMATRX_PART2(I,3*(J-1)+3)
          END DO
          END DO
          
          DO J = 1, 48
          DO I = 1, 24
              AMATRX_PART3(I,J) = AMATRX_SEMI_OLD(I,J)
              AMATRX_PART4(I,J) = AMATRX_SEMI_OLD(I+24,J)
          END DO
          END DO
          
          DO J = 1, 48
          DO I = 1, 8
              AMATRX_NEW(6*(I-1)+1,J) = AMATRX_PART3(3*(I-1)+1,J)
              AMATRX_NEW(6*(I-1)+2,J) = AMATRX_PART3(3*(I-1)+2,J)
              AMATRX_NEW(6*(I-1)+3,J) = AMATRX_PART3(3*(I-1)+3,J)
              AMATRX_NEW(6*(I-1)+4,J) = AMATRX_PART4(3*(I-1)+1,J)
              AMATRX_NEW(6*(I-1)+5,J) = AMATRX_PART4(3*(I-1)+2,J)
              AMATRX_NEW(6*(I-1)+6,J) = AMATRX_PART4(3*(I-1)+3,J)
          END DO
          END DO
          

          DO I = 1, 24
              RHS_PART1(I,1) = RHS_OLD(I,1)
              RHS_PART2(I,1) = RHS_OLD(I+24,1)
          END DO
          
          DO I = 1, 8
              RHS_NEW(6*(I-1)+1,1) = RHS_PART1(3*(I-1)+1,1)
              RHS_NEW(6*(I-1)+2,1) = RHS_PART1(3*(I-1)+2,1)
              RHS_NEW(6*(I-1)+3,1) = RHS_PART1(3*(I-1)+3,1)
              RHS_NEW(6*(I-1)+4,1) = RHS_PART2(3*(I-1)+1,1)
              RHS_NEW(6*(I-1)+5,1) = RHS_PART2(3*(I-1)+2,1)
              RHS_NEW(6*(I-1)+6,1) = RHS_PART2(3*(I-1)+3,1)
          END DO
C
C
C
          !print *, 'RHS_NEW:'
          !CALL PRTMAT(RHS_NEW,48,1)
          !
          !print *, 'AMATRX_NEW:'
          !CALL PRTMAT(AMATRX_NEW,48,48)
          !print *, 'CHOS GOOOZ AKHAR:'
          
      RETURN 
     	END
          
C          
C%*************************************************************************
C% Compute B1_COS,B2_COS,B3_COS,B4_COS,B5_COS
C%*************************************************************************
C%%
		SUBROUTINE MAKER_B_MATRIXS_COS_HEX8(SF,SHPD,F,B1_COS,B2_COS,B3_COS,
     +		  B4_COS,B5_COS,B6_COS,B1_COS_TRANS,B2_COS_TRANS,B3_COS_TRANS,
     +		  B4_COS_TRANS,B5_COS_TRANS,B6_COS_TRANS)  
     		
        IMPLICIT REAL*8 (A-H,O-Z)
        
        INTEGER I, J, K, COL1, COL2, COL3
          
        REAL*8 SHPD(3,8), F(3,3), F_INV(3,3), Z3, 
     +		B1_COS(9,24), B2_COS(9,24), B3_COS(3,24), B4_COS(3,24),
     +		B5_COS(9,24), B6_COS(9,24), B1_COS_TRANS(24,9),
     +		B2_COS_TRANS(24,9), B3_COS_TRANS(24,3), B4_COS_TRANS(24,3),
     +		B5_COS_TRANS(24,9), B6_COS_TRANS(24,9), SF(8,1)
        
        PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )
     
     	CALL MATINV(F,F_INV,Z3)
     	
     	DO I = 1, 9
    	DO J = 1, 24
	    	B1_COS(I,J) = ZERO
	    	B2_COS(I,J) = ZERO
	    	B5_COS(I,J) = ZERO
	    	B6_COS(I,J) = ZERO
		END DO
		END DO
		
		DO I = 1, 8
			COL1 = (I-1)*3+1
			COL2 = (I-1)*3+2
			COL3 = (I-1)*3+3
				
			B1_COS(1,COL1) = SHPD(1,I)
			B1_COS(2,COL2) = SHPD(2,I)			
			B1_COS(3,COL3) = SHPD(3,I)
			B1_COS(4,COL1) = SHPD(2,I)
			B1_COS(5,COL2) = SHPD(1,I)
			B1_COS(6,COL2) = SHPD(3,I)
			B1_COS(7,COL3) = SHPD(2,I)
			B1_COS(8,COL1) = SHPD(3,I)
			B1_COS(9,COL3) = SHPD(1,I)
			
			B2_COS(1,COL1) = F(1,1)*SHPD(1,I)
			B2_COS(1,COL2) = F(2,1)*SHPD(1,I)
			B2_COS(1,COL3) = F(3,1)*SHPD(1,I)
			
			B2_COS(2,COL1) = F(1,2)*SHPD(2,I)
			B2_COS(2,COL2) = F(2,2)*SHPD(2,I)
			B2_COS(2,COL3) = F(3,2)*SHPD(2,I)
			
			B2_COS(3,COL1) = F(1,3)*SHPD(3,I)
			B2_COS(3,COL2) = F(2,3)*SHPD(3,I)
			B2_COS(3,COL3) = F(3,3)*SHPD(3,I)
			
			B2_COS(4,COL1) = F(1,1)*SHPD(2,I)
			B2_COS(4,COL2) = F(2,1)*SHPD(2,I)
			B2_COS(4,COL3) = F(3,1)*SHPD(2,I)
			
			B2_COS(5,COL1) = F(1,2)*SHPD(1,I)
			B2_COS(5,COL2) = F(2,2)*SHPD(1,I)
			B2_COS(5,COL3) = F(3,2)*SHPD(1,I)
			
			B2_COS(6,COL1) = F(1,2)*SHPD(3,I)
			B2_COS(6,COL2) = F(2,2)*SHPD(3,I)
			B2_COS(6,COL3) = F(3,2)*SHPD(3,I)
			
			B2_COS(7,COL1) = F(1,3)*SHPD(2,I)
			B2_COS(7,COL2) = F(2,3)*SHPD(2,I)
			B2_COS(7,COL3) = F(3,3)*SHPD(2,I)
			
			B2_COS(8,COL1) = F(1,1)*SHPD(3,I)
			B2_COS(8,COL2) = F(2,1)*SHPD(3,I)
			B2_COS(8,COL3) = F(3,1)*SHPD(3,I)
			
			B2_COS(9,COL1) = F(1,3)*SHPD(1,I)
			B2_COS(9,COL2) = F(2,3)*SHPD(1,I)
			B2_COS(9,COL3) = F(3,3)*SHPD(1,I)
			
			B5_COS(1,COL1) = SHPD(1,I)*F_INV(1,1) + 
     +			SHPD(2,I)*F_INV(2,1) + SHPD(3,I)*F_INV(3,1)
     
			B5_COS(2,COL2) = SHPD(1,I)*F_INV(1,2) + 
     +			SHPD(2,I)*F_INV(2,2) + SHPD(3,I)*F_INV(3,2)	
	 	
			B5_COS(3,COL3) = SHPD(1,I)*F_INV(1,3) + 
     +			SHPD(2,I)*F_INV(2,3) + SHPD(3,I)*F_INV(3,3)
	 
			B5_COS(4,COL1) = SHPD(1,I)*F_INV(1,2) + 
     +			SHPD(2,I)*F_INV(2,2) + SHPD(3,I)*F_INV(3,2)
	 
			B5_COS(5,COL2) = SHPD(1,I)*F_INV(1,1) + 
     +			SHPD(2,I)*F_INV(2,1) + SHPD(3,I)*F_INV(3,1)
	 
			B5_COS(6,COL2) = SHPD(1,I)*F_INV(1,3) + 
     +			SHPD(2,I)*F_INV(2,3) + SHPD(3,I)*F_INV(3,3)
	 
			B5_COS(7,COL3) = SHPD(1,I)*F_INV(1,2) + 
     +			SHPD(2,I)*F_INV(2,2) + SHPD(3,I)*F_INV(3,2)
	 
			B5_COS(8,COL1) = SHPD(1,I)*F_INV(1,3) + 
     +			SHPD(2,I)*F_INV(2,3) + SHPD(3,I)*F_INV(3,3)
	 
			B5_COS(9,COL3) = SHPD(1,I)*F_INV(1,1) + 
     +			SHPD(2,I)*F_INV(2,1) + SHPD(3,I)*F_INV(3,1)
			
			B6_COS(1,COL1) = SHPD(1,I)*F_INV(1,1) + 
     +			SHPD(2,I)*F_INV(2,1) + SHPD(3,I)*F_INV(3,1)
	 
			B6_COS(2,COL2) = SHPD(1,I)*F_INV(1,2) + 
     +			SHPD(2,I)*F_INV(2,2) + SHPD(3,I)*F_INV(3,2)	
	 	
			B6_COS(3,COL3) = SHPD(1,I)*F_INV(1,3) + 
     +			SHPD(2,I)*F_INV(2,3) + SHPD(3,I)*F_INV(3,3)
	 
			B6_COS(4,COL2) = SHPD(1,I)*F_INV(1,1) + 
     +			SHPD(2,I)*F_INV(2,1) + SHPD(3,I)*F_INV(3,1)
	 
			B6_COS(5,COL1) = SHPD(1,I)*F_INV(1,2) + 
     +			SHPD(2,I)*F_INV(2,2) + SHPD(3,I)*F_INV(3,2)
	 
			B6_COS(6,COL3) = SHPD(1,I)*F_INV(1,2) +
     +			SHPD(2,I)*F_INV(2,2) + SHPD(3,I)*F_INV(3,2)
	 
			B6_COS(7,COL2) = SHPD(1,I)*F_INV(1,3) + 
     +			SHPD(2,I)*F_INV(2,3) + SHPD(3,I)*F_INV(3,3)
	 
			B6_COS(8,COL3) = SHPD(1,I)*F_INV(1,1) + 
     +			SHPD(2,I)*F_INV(2,1) + SHPD(3,I)*F_INV(3,1)
	 
			B6_COS(9,COL1) = SHPD(1,I)*F_INV(1,3) + 
     +			SHPD(2,I)*F_INV(2,3) + SHPD(3,I)*F_INV(3,3)
		
			
		END DO
		
		DO I = 1, 3
    	DO J = 1, 24
	    	B3_COS(I,J) = ZERO
	    	B4_COS(I,J) = ZERO
		END DO
		END DO
		
		DO I = 1, 8
			COL1 = (I-1)*3+1
			COL2 = (I-1)*3+2
			COL3 = (I-1)*3+3
				
			B3_COS(1,COL1) = SF(I,1)
			B3_COS(2,COL2) = SF(I,1)			
			B3_COS(3,COL3) = SF(I,1)
			
			B4_COS(1,COL2) = SHPD(1,I)*F_INV(1,3) + 
     +			SHPD(2,I)*F_INV(2,3) + SHPD(3,I)*F_INV(3,3)
	 
			B4_COS(1,COL3) = -SHPD(1,I)*F_INV(1,2) - 
     +			SHPD(2,I)*F_INV(2,2) - SHPD(3,I)*F_INV(3,2)
			
			B4_COS(2,COL1) = -SHPD(1,I)*F_INV(1,3) - 
     +			SHPD(2,I)*F_INV(2,3) - SHPD(3,I)*F_INV(3,3)
	 
			B4_COS(2,COL3) = SHPD(1,I)*F_INV(1,1) + 
     +			SHPD(2,I)*F_INV(2,1) + SHPD(3,I)*F_INV(3,1)
			
			B4_COS(3,COL1) = SHPD(1,I)*F_INV(1,2) + 
     +			SHPD(2,I)*F_INV(2,2) + SHPD(3,I)*F_INV(3,2)
	 
			B4_COS(3,COL2) = -SHPD(1,I)*F_INV(1,1) - 
     +			SHPD(2,I)*F_INV(2,1) - SHPD(3,I)*F_INV(3,1)
			
			
		END DO

		DO I = 1, 24
    	    DO J = 1, 9
	    	B1_COS_TRANS(I,J) = B1_COS(J,I)
	    	B2_COS_TRANS(I,J) = B2_COS(J,I)
	    	B5_COS_TRANS(I,J) = B5_COS(J,I)
	    	B6_COS_TRANS(I,J) = B6_COS(J,I)
		END DO
		END DO
		
		DO I = 1, 24
    	    DO J = 1, 3
	    	B3_COS_TRANS(I,J) = B3_COS(J,I)
	    	B4_COS_TRANS(I,J) = B4_COS(J,I)
		END DO
		END DO
		
		RETURN 
     	END
C
C          
C%*************************************************************************
C% Compute B4_COS_QUAD 3*60 AND TRANSPOSED
C%*************************************************************************
C%%
		SUBROUTINE MAKER_B4_COS_QUAD_HEX20(SHPD_QUAD,F_QUAD,
     +		  B4_COS_QUAD, B4_COS_TRANS_QUAD)   
     		
        IMPLICIT REAL*8 (A-H,O-Z)
        
        INTEGER I, J, K, COL1, COL2, COL3
          
        REAL*8 SHPD_QUAD(3,20), F_QUAD(3,3), F_INV(3,3), Z3,  
     +		B4_COS_QUAD(3,60), B4_COS_TRANS_QUAD(60,3)
        
        PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )
     
     	CALL MATINV(F_QUAD,F_INV,Z3)     	
     		
		DO I = 1, 3
    	    DO J = 1, 60
	    	B4_COS_QUAD(I,J) = ZERO
		END DO
		END DO
		
		DO I = 1, 20
			COL1 = (I-1)*3+1
			COL2 = (I-1)*3+2
			COL3 = (I-1)*3+3
							
			B4_COS_QUAD(1,COL2) = SHPD_QUAD(1,I)*F_INV(1,3) + 
     +			SHPD_QUAD(2,I)*F_INV(2,3) + SHPD_QUAD(3,I)*F_INV(3,3)
	 
			B4_COS_QUAD(1,COL3) = -SHPD_QUAD(1,I)*F_INV(1,2) - 
     +			SHPD_QUAD(2,I)*F_INV(2,2) - SHPD_QUAD(3,I)*F_INV(3,2)
			
			B4_COS_QUAD(2,COL1) = -SHPD_QUAD(1,I)*F_INV(1,3) - 
     +			SHPD_QUAD(2,I)*F_INV(2,3) - SHPD_QUAD(3,I)*F_INV(3,3)
	 
			B4_COS_QUAD(2,COL3) = SHPD_QUAD(1,I)*F_INV(1,1) + 
     +			SHPD_QUAD(2,I)*F_INV(2,1) + SHPD_QUAD(3,I)*F_INV(3,1)
			
			B4_COS_QUAD(3,COL1) = SHPD_QUAD(1,I)*F_INV(1,2) + 
     +			SHPD_QUAD(2,I)*F_INV(2,2) + SHPD_QUAD(3,I)*F_INV(3,2)
	 
			B4_COS_QUAD(3,COL2) = -SHPD_QUAD(1,I)*F_INV(1,1) - 
     +			SHPD_QUAD(2,I)*F_INV(2,1) - SHPD_QUAD(3,I)*F_INV(3,1)
		
		END DO
		
		DO I = 1, 60
    	    DO J = 1, 3
	    	B4_COS_TRANS_QUAD(I,J) = B4_COS_QUAD(J,I)
		END DO
		END DO
		
      RETURN 
     	END
C          
C%*************************************************************************
C% Compute BN AND BG MATRIX FROM DERIVATIVE OF SHAPE FUNCTIONS AND DEFORMATION GRADIENTS
C%*************************************************************************
C%%
        SUBROUTINE MAKER_BNBG_HEX8(SHPD,F,BN,BG)
		
        IMPLICIT REAL*8 (A-H,O-Z)
        
        INTEGER I, J, K, COL1, COL2, COL3
          
        REAL*8 SHPD(3,8), F(3,3), BN(6,24), BG(9,24)
        
        PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )
     
     		DO I = 1, 6
    	    DO J = 1, 24
	    	    BN(I,J) = ZERO
			END DO
			END DO			
C
			DO I = 1, 9
    	    DO J = 1, 24
	    	    BG(I,J) = ZERO
			END DO
			END DO
C
			DO I = 1, 8
				COL1 = (I-1)*3+1
				COL2 = (I-1)*3+2
				COL3 = (I-1)*3+3
				
				BN(1,COL1) = SHPD(1,I)*F(1,1)
				BN(2,COL1) = SHPD(2,I)*F(1,2)
				BN(3,COL1) = SHPD(3,I)*F(1,3)
				BN(4,COL1) = SHPD(1,I)*F(1,2)+SHPD(2,I)*F(1,1)
				BN(5,COL1) = SHPD(2,I)*F(1,3)+SHPD(3,I)*F(1,2)
				BN(6,COL1) = SHPD(1,I)*F(1,3)+SHPD(3,I)*F(1,1)
				
				BN(1,COL2) = SHPD(1,I)*F(2,1)
				BN(2,COL2) = SHPD(2,I)*F(2,2)
				BN(3,COL2) = SHPD(3,I)*F(2,3)
				BN(4,COL2) = SHPD(1,I)*F(2,2)+SHPD(2,I)*F(2,1)
				BN(5,COL2) = SHPD(2,I)*F(2,3)+SHPD(3,I)*F(2,2)
				BN(6,COL2) = SHPD(1,I)*F(2,3)+SHPD(3,I)*F(2,1)
				
				BN(1,COL3) = SHPD(1,I)*F(3,1)
				BN(2,COL3) = SHPD(2,I)*F(3,2)
				BN(3,COL3) = SHPD(3,I)*F(3,3)
				BN(4,COL3) = SHPD(1,I)*F(3,2)+SHPD(2,I)*F(3,1)
				BN(5,COL3) = SHPD(2,I)*F(3,3)+SHPD(3,I)*F(3,2)
				BN(6,COL3) = SHPD(1,I)*F(3,3)+SHPD(3,I)*F(3,1)
C----------------------------------
				BG(1,COL1) = SHPD(1,I)	
				BG(2,COL1) = SHPD(2,I)	
				BG(3,COL1) = SHPD(3,I)		
				
				BG(4,COL2) = SHPD(1,I)	
				BG(5,COL2) = SHPD(2,I)	
				BG(6,COL2) = SHPD(3,I)	
				
				BG(7,COL3) = SHPD(1,I)	
				BG(8,COL3) = SHPD(2,I)	
				BG(9,COL3) = SHPD(3,I)	
				
			END DO     	
     	
     	RETURN 
     	END

C**
C*************************************************
C**
C**  UMAT FOR UPDATED_TUCKER_COSSERAT MODEL
C**   
C*************************************************
C   
      SUBROUTINE UMAT_COSSERAT_TUCKER(PROPERTIES, F, KAPPA_CURVE, 
     +       STRESS_VECT, MATERIAL_MATRIX, COSSERAT_MATRIX, DWDKAPPA)
C
C      INCLUDE 'ABA_PARAM.INC'
C
C     CHARACTER*80 CMNAME
C      
C      DIMENSION STRESS(NTENS),STATEV(NSTATV),
C     +    DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
C     +    STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
C     +    PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
           
C           user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
C           and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT      
C		
C
C**********************************************************************
C             DEFINITION OF PARAMETERS
C**********************************************************************     
C 
C       
		INTEGER I, J, K, K1, K2
C     
		REAL*8   PROPERTIES(14), F(3,3), STRESS_VECT(6), MATERIAL_MATRIX(6,6),
     +      AIDENT(3,3), F_TAU(3,3), F_TAU_TRANS(3,3), KAPPA_CURVE(3,3),
     +		K3_COS, E0(3,1), N0(3,1), EK, EG, EES, PHI_0, PHI_A,
     +		A_S, EG_0, E0_E0(3,3), N0_N0(3,3), C(3,3), C_INV(3,3),
     +		Z3, I1_INVAR, I3_INVAR, J_INVAR, L_INVAR, LAPLAS_INVAR,
     +		DI1DC(3,3), DI3DC(3,3), DJDC(3,3), DLDC(3,3), DW0DLAPLAS,
     +		D2W0DLAPLAS2, LAPLASDC(3,3), D2LAPLASDCDC(3,3,3,3),
     +		DI6DKAPPA(3,3), DWDKAPPA(3,3), DDI6DKAPPADKAPPA(3,3,3,3),
     +		E_VEC(3,1), N_VEC(3,1), DWDC(3,3), F_TAU_INV(3,3),
     +		S_PIOLA(3,3), SIGMA_CAUCHY(3,3), Z1(3,3), Z2(3,3),
     +		QQ_TENSOR(3,3,3,3), C_SPATIAL_TAU(3,3,3,3), C_PRIME(3,3,3,3),
     +		C_ABAQUS(3,3,3,3), D2WDCDC(3,3,3,3), C_MAT_TAU(3,3,3,3),
     +		D2WDKAPPADKAPPA(3,3,3,3), COSSERAT_MATRIX(9,9), DLAPLASDC(3,3)
C         
		PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, SIX=6.D0, ONE_FOURTH=1.D0/4.D0 )
C          
C
C**********************************************************************
C             BEGIN COMPUTATION
C**********************************************************************     
C ----------------------------------------------------------------
C    UMAT FOR MOONEY-RIVLIN HYPERELASTICITY
C    ONLY FOR 3D CASE 
C ----------------------------------------------------------------
C 		INPUTS:
C	  F(3,3) : DEFORMATION GRADIENRT 
C     PROPS(1) - A10
C     PROPS(2) - A01 
C     PROPS(3) - K                   
C
C		 Outputs:
C	%  Stress = 2nd PK stress [S11, S22, S33, S12, S23, S13];
C	%  D = Material stiffness [6x6]    
C  
C**********************************************************************
C             Material properties
C**********************************************************************
C
C
C     E0 :=> THE FIBER DIRECTIONS
C
      	E0(1,1) = PROPERTIES(1)
      	E0(2,1) = PROPERTIES(2)
      	E0(3,1) = PROPERTIES(3)
C
C     N0 :=> UNIT VECTOR NORMAL TO THE PLANE OF THE UNDEFORMED SHEET
C
      	N0(1,1) = PROPERTIES(4)
      	N0(2,1) = PROPERTIES(5)
      	N0(3,1) = PROPERTIES(6)
      
C
C     RELATED TO TRANSVERSLY HYPERELASTIC MODEL
C
      	EK = PROPERTIES(7)
      	EG = PROPERTIES(8)
      	EES = PROPERTIES(9)
      	PHI_0 = PROPERTIES(10)
      	PHI_A = PROPERTIES(11) ! MAXIMUM POSSIBLE FIBER VOLUME FRACTION
      	A_S = PROPERTIES(12)  ! THE SPRING CONSTANT
      	EG_0 = PROPERTIES(13) ! A SMALL PORTION OF EG
      	K3_COS = PROPERTIES(14) ! COSSERAT MATERIAL PARAMETER BY I6
       
          !PRINT *, "IN UMAT  --> ","PROPERTIES = "
          !CALL PRTMAT(PROPERTIES,14,1) 
             
      	CALL ONEM(AIDENT)
      	
      	CALL ZEROM(E0_E0)
      	E0_E0(1,1) = E0(1,1)*E0(1,1)
      	E0_E0(1,2) = E0(1,1)*E0(2,1)
      	E0_E0(1,3) = E0(1,1)*E0(3,1)
      	E0_E0(2,1) = E0(2,1)*E0(1,1)
      	E0_E0(2,2) = E0(2,1)*E0(2,1)
      	E0_E0(2,3) = E0(2,1)*E0(3,1)
      	E0_E0(3,1) = E0(3,1)*E0(1,1)
      	E0_E0(3,2) = E0(3,1)*E0(2,1)
      	E0_E0(3,3) = E0(3,1)*E0(3,1)
      
      	CALL ZEROM(N0_N0)
      	N0_N0(1,1) = N0(1,1)*N0(1,1)
     	    N0_N0(1,2) = N0(1,1)*N0(2,1)
      	N0_N0(1,3) = N0(1,1)*N0(3,1)
      	N0_N0(2,1) = N0(2,1)*N0(1,1)
      	N0_N0(2,2) = N0(2,1)*N0(2,1)
      	N0_N0(2,3) = N0(2,1)*N0(3,1)
      	N0_N0(3,1) = N0(3,1)*N0(1,1)
      	N0_N0(3,2) = N0(3,1)*N0(2,1)
      	N0_N0(3,3) = N0(3,1)*N0(3,1) 
          
          !PRINT *, "E0_E0 = "
          !CALL PRTMAT(E0_E0,3,3) 
          !PRINT *, "N0_N0 = "
          !CALL PRTMAT(N0_N0,3,3) 
C      
C
C**********************************************************************
C             BEGIN MAIN LOOP COMPUTATION
C**********************************************************************         
C
C	INITIALIZE DEFORMATION GRADIENT
C       

          	F_TAU(1,1) = F(1,1)
          	F_TAU(2,2) = F(2,2)
          	F_TAU(3,3) = F(3,3)
          	F_TAU(1,2) = F(1,2)
			F_TAU(2,3) = F(2,3)
			F_TAU(3,1) = F(3,1)
			F_TAU(2,1) = F(2,1)
			F_TAU(3,2) = F(3,2)
			F_TAU(1,3) = F(1,3)
			
			
C
C     CALCULATE MATRIX C
C     
      	CALL MTRANS(F_TAU,F_TAU_TRANS)
      
      	CALL MPROD(F_TAU_TRANS,F_TAU,C) 
      
      	CALL MATINV(C,C_INV,Z3)
          
          !PRINT *, "F_TAU = "
          !CALL PRTMAT(F_TAU,3,3) 
          !PRINT *, "F_TAU_TRANS = "
          !CALL PRTMAT(F_TAU_TRANS,3,3) 
          !PRINT *, "C = "
          !CALL PRTMAT(C,3,3) 
          !PRINT *, "C_INV = "
          !CALL PRTMAT(C_INV,3,3) 
          !PRINT *, "C DETERMINED = " , Z3
      
C
C     CALCULATE DEFORMATION INVARIENTS: I1,I3 OR J, L_INVAR , LAPLAS_INVAR
C 
      	I1_INVAR = C(1,1) + C(2,2) + C(3,3)
      
      	CALL MDET(C,I3_INVAR)
      
      	J_INVAR = DSQRT(I3_INVAR)
      
      	L_INVAR = 
     +    E0_E0(1,1)*C(1,1) + E0_E0(1,2)*C(1,2) + E0_E0(1,3)*C(1,3) +
     +    E0_E0(2,1)*C(2,1) + E0_E0(2,2)*C(2,2) + E0_E0(2,3)*C(2,3) +
     +    E0_E0(3,1)*C(3,1) + E0_E0(3,2)*C(3,2) + E0_E0(3,3)*C(3,3)
      
      	LAPLAS_INVAR = 
     +    N0_N0(1,1)*C_INV(1,1) + N0_N0(1,2)*C_INV(1,2) + 
     +    N0_N0(1,3)*C_INV(1,3) + N0_N0(2,1)*C_INV(2,1) +
     +    N0_N0(2,2)*C_INV(2,2) + N0_N0(2,3)*C_INV(2,3) +
     +    N0_N0(3,1)*C_INV(3,1) + N0_N0(3,2)*C_INV(3,2) + 
     +    N0_N0(3,3)*C_INV(3,3)
          
          !PRINT *, "C = "
          !CALL PRTMAT(C,3,3) 
          !PRINT *, "C_INV = "
          !CALL PRTMAT(C_INV,3,3) 
          !PRINT *, "I1_INVAR = " , I1_INVAR
          !PRINT *, "I3_INVAR = " , I3_INVAR
          !PRINT *, "J_INVAR = " , J_INVAR
          !PRINT *, "L_INVAR = " , L_INVAR
          !PRINT *, "LAPLAS_INVAR = " , LAPLAS_INVAR
      
C
C     CALCULATE DERIVATIVE OF INVARIENTS WITH RESPECT TO C
C 
      	CALL ZEROM(DI1DC)
      	DI1DC(1,1) = AIDENT(1,1)
      	DI1DC(1,2) = AIDENT(1,2)
      	DI1DC(1,3) = AIDENT(1,3)
      	DI1DC(2,1) = AIDENT(2,1)
      	DI1DC(2,2) = AIDENT(2,2)
      	DI1DC(2,3) = AIDENT(2,3)
      	DI1DC(3,1) = AIDENT(3,1)
      	DI1DC(3,2) = AIDENT(3,2)
      	DI1DC(3,3) = AIDENT(3,3)
      
      	CALL ZEROM(DI3DC)
      	DI3DC(1,1) = I3_INVAR*C_INV(1,1)
      	DI3DC(1,2) = I3_INVAR*C_INV(1,2)
      	DI3DC(1,3) = I3_INVAR*C_INV(1,3)
      	DI3DC(2,1) = I3_INVAR*C_INV(2,1)
      	DI3DC(2,2) = I3_INVAR*C_INV(2,2)
      	DI3DC(2,3) = I3_INVAR*C_INV(2,3)
      	DI3DC(3,1) = I3_INVAR*C_INV(3,1)
      	DI3DC(3,2) = I3_INVAR*C_INV(3,2)
      	DI3DC(3,3) = I3_INVAR*C_INV(3,3)
      
      	CALL ZEROM(DJDC)
      	DJDC(1,1) = ONE_HALF*J_INVAR*C_INV(1,1)
      	DJDC(1,2) = ONE_HALF*J_INVAR*C_INV(1,2)
      	DJDC(1,3) = ONE_HALF*J_INVAR*C_INV(1,3)
      	DJDC(2,1) = ONE_HALF*J_INVAR*C_INV(2,1)
      	DJDC(2,2) = ONE_HALF*J_INVAR*C_INV(2,2)
      	DJDC(2,3) = ONE_HALF*J_INVAR*C_INV(2,3)
      	DJDC(3,1) = ONE_HALF*J_INVAR*C_INV(3,1)
      	DJDC(3,2) = ONE_HALF*J_INVAR*C_INV(3,2)
      	DJDC(3,3) = ONE_HALF*J_INVAR*C_INV(3,3)
      
      	CALL ZEROM(DLDC)
      	DLDC(1,1) = E0_E0(1,1)
      	DLDC(1,2) = E0_E0(1,2)
      	DLDC(1,3) = E0_E0(1,3)
      	DLDC(2,1) = E0_E0(2,1)
      	DLDC(2,2) = E0_E0(2,2)
      	DLDC(2,3) = E0_E0(2,3)
      	DLDC(3,1) = E0_E0(3,1)
      	DLDC(3,2) = E0_E0(3,2)
      	DLDC(3,3) = E0_E0(3,3)
          
          !PRINT *, "DI1DC = "
          !CALL PRTMAT(DI1DC,3,3) 
          !PRINT *, "DI3DC = "
          !CALL PRTMAT(DI3DC,3,3) 
          !PRINT *, "DJDC = "
          !CALL PRTMAT(DJDC,3,3)
          !PRINT *, "DLDC = "
          !CALL PRTMAT(DLDC,3,3)


      	DW0DLAPLAS = (A_S/TWO)*(
     +    (LAPLAS_INVAR-DSQRT(LAPLAS_INVAR))/
     +    ((ONE/PHI_0 - DSQRT(LAPLAS_INVAR)/PHI_A)**FOUR) ) + 
     +    ((EG-EG_0)/TWO)*(ONE/(LAPLAS_INVAR**TWO) - ONE/LAPLAS_INVAR) 
C
      	D2W0DLAPLAS2 = (A_S/TWO)*(
     +    ( (ONE-ONE_HALF/DSQRT(LAPLAS_INVAR))*
     +    (ONE/PHI_0 - DSQRT(LAPLAS_INVAR)/PHI_A) - 
     +    (TWO/PHI_A)*(DSQRT(LAPLAS_INVAR)-ONE) )/   
     +    ( (ONE/PHI_0 - DSQRT(LAPLAS_INVAR)/PHI_A)**FIVE ) ) + 
     +    ONE_HALF*(EG-EG_0)*( -TWO/(LAPLAS_INVAR**THREE) + 
     +    ONE/(LAPLAS_INVAR**TWO) )       
C      
          DO 51 P = 1,3
	    DO 51 Q = 1,3
              DLAPLASDC(P,Q) = ZERO
	        DO 52 I = 1,3
	        DO 52 J = 1,3
	        DLAPLASDC(P,Q) = DLAPLASDC(P,Q) - N0(I,1)*N0(J,1)*
     +        (C_INV(I,P)*C_INV(Q,J))        
52	    CONTINUE
51	    CONTINUE

		DO 69 I = 1,3
	    DO 69 J = 1,3
	    DO 69 K = 1,3
	    DO 69 L = 1,3
              D2LAPLASDCDC(I,J,K,L) = ZERO
	        DO 70 P = 1,3
	        DO 70 Q = 1,3
	        D2LAPLASDCDC(I,J,K,L) = D2LAPLASDCDC(I,J,K,L)- N0(P,1)*N0(Q,1)*
     +        (C_INV(P,K)*C_INV(L,I)*C_INV(J,Q) +
     +		   C_INV(P,I)*C_INV(J,K)*C_INV(L,Q))  
70	    CONTINUE
          !PRINT *, "D2LAPLASDCDC(I,J,K,L) " , D2LAPLASDCDC(I,J,K,L)
69	    CONTINUE
          
          !PRINT *, "DW0DLAPLAS = " , DW0DLAPLAS
          !PRINT *, "D2W0DLAPLAS2 = " , D2W0DLAPLAS2
          !
          !PRINT *, "DLAPLASDC = "
          !CALL PRTMAT(DLAPLASDC,3,3) 
          
          !PRINT *, "KAPPA_CURVE = "
          !CALL PRTMAT(KAPPA_CURVE,3,3) 
          
          !PRINT *, "E0 = "
          !CALL PRTMAT(E0,3,1) 

C
C		CALCULATE DI6DKAPPA 
C
		CALL ZEROM(DI6DKAPPA)
		DO I=1,3
          DO J=1,3
        	DI6DKAPPA(I,J) = ZERO
        	    DO K=1,3
        		    DI6DKAPPA(I,J) = DI6DKAPPA(I,J)+
     +                TWO*KAPPA_CURVE(I,K)*E0(K,1)*E0(J,1)
              END DO
          END DO
          END DO
          
C
C		CALCULATE DWDKAPPA 
C
		CALL ZEROM(DWDKAPPA)
		DO I=1,3
          DO J=1,3
        	    DWDKAPPA(I,J) = K3_COS*DI6DKAPPA(I,J)
          END DO
          END DO  
        
          !PRINT *, "DI6DKAPPA = "
          !CALL PRTMAT(DI6DKAPPA,3,3) 
          !PRINT *, "DWDKAPPA = "
          !CALL PRTMAT(DWDKAPPA,3,3) 

C
C		CALCULATE DDI6DKAPPADKAPPA
C

		DO I=1,3
          DO J=1,3
          DO K=1,3
          DO L=1,3
        	    DDI6DKAPPADKAPPA(I,J,K,L) = TWO*AIDENT(I,K)*E0(L,1)*E0(J,1)
          !PRINT *,"DDI6DKAPPADKAPPA(I,J,K,L)=",DDI6DKAPPADKAPPA(I,J,K,L)
          END DO
          END DO
          END DO
          END DO

C
C     CALCULATE E_VEC = E0*F_TAU_TRANS/DSQRT(L_INVAR)
C
      	E_VEC(1,1) = ( E0(1,1)*F_TAU_TRANS(1,1) +E0(2,1)*F_TAU_TRANS(2,1)
     +    + E0(3,1)*F_TAU_TRANS(3,1) )/DSQRT(L_INVAR)
      
      	E_VEC(2,1) = ( E0(1,1)*F_TAU_TRANS(1,2) +E0(2,1)*F_TAU_TRANS(2,2)
     +    + E0(3,1)*F_TAU_TRANS(3,2) )/DSQRT(L_INVAR)

      	E_VEC(3,1) = ( E0(1,1)*F_TAU_TRANS(1,3) +E0(2,1)*F_TAU_TRANS(2,3)
     +    + E0(3,1)*F_TAU_TRANS(3,3) )/DSQRT(L_INVAR)

C
C     CALCULATE N_VEC = N0*F_INV/DSQRT(LAPLAS_INVAR)
C
C
C     N0 :=> UNIT VECTOR NORMAL TO THE PLANE OF THE UNDEFORMED SHEET
C
      	CALL MATINV(F_TAU,F_TAU_INV,Z3)
      
      	N_VEC(1,1) = ( N0(1,1)*F_TAU_INV(1,1) + N0(2,1)*F_TAU_INV(2,1)
     +    + N0(3,1)*F_TAU_INV(3,1) )/DSQRT(LAPLAS_INVAR)
      
      	N_VEC(2,1) = ( N0(1,1)*F_TAU_INV(1,2) + N0(2,1)*F_TAU_INV(2,2)
     +    + N0(3,1)*F_TAU_INV(3,2) )/DSQRT(LAPLAS_INVAR)

      	N_VEC(3,1) = ( N0(1,1)*F_TAU_INV(1,3) + N0(2,1)*F_TAU_INV(2,3)
     +    + N0(3,1)*F_TAU_INV(3,3) )/DSQRT(LAPLAS_INVAR)
          
          !PRINT *, "E_VEC = "
          !CALL PRTMAT(E_VEC,3,1) 
          !PRINT *, "N_VEC = "
          !CALL PRTMAT(N_VEC,3,1) 
          !PRINT *, "F_TAU_INV = "
          !CALL PRTMAT(F_TAU_INV,3,3)
               
C
C     CALCULATE THE DWDC = (1/2)*S_PIOLA
C
      	DO K1=1, 3
      	DO K2=1, 3
          DWDC(K1,K2) = DW0DLAPLAS*DLAPLASDC(K1,K2) +
     +        ONE_FOURTH*EES*PHI_0*(ONE - ONE/L_INVAR)*E0_E0(K1,K2) +
     +        ONE_HALF*EG*(AIDENT(K1,K2)-C_INV(K1,K2)) + 
     +        ONE_HALF*EK*DLOG(J_INVAR)*C_INV(K1,K2)
      	END DO
      	END DO
C
C     CALCULATE THE SECOND PIOLA-KIRCHHOF STRESS: S
C
      	DO K1=1, 3
      	DO K2=1, 3
          S_PIOLA(K1,K2) = TWO*DW0DLAPLAS*DLAPLASDC(K1,K2) +
     +        ONE_HALF*EES*PHI_0*(ONE - ONE/L_INVAR)*E0_E0(K1,K2) +
     +        EG*(AIDENT(K1,K2)-C_INV(K1,K2)) + 
     +        EK*DLOG(J_INVAR)*C_INV(K1,K2)
      	END DO
      	END DO
          
          !PRINT *, "DLAPLASDC = "
          !CALL PRTMAT(DLAPLASDC,3,3) 
          !PRINT *, "E0_E0 = "
          !CALL PRTMAT(E0_E0,3,3)
          !PRINT *, "C_INV = "
          !CALL PRTMAT(C_INV,3,3) 
          !PRINT *, "J_INVAR = " , J_INVAR
          !PRINT *, "L_INVAR = " , L_INVAR
C
C      CACLULATE SIGMA_CAUCHY(I,J) := (1/J)F*S_PIOLA*F_TRANS
C
      	CALL ZEROM(Z1)
      	CALL ZEROM(Z2)
      	CALL ZEROM(SIGMA_CAUCHY)
      	CALL MPROD(S_PIOLA,F_TAU_TRANS,Z1)
      	CALL MPROD(F_TAU,Z1,Z2)
      	DO K1=1, 3
          DO K2=1, 3
              SIGMA_CAUCHY(K1,K2) = (ONE/J_INVAR)*Z2(K1,K2)
          END DO
      	END DO  
          
          !PRINT *, "DWDC = "
          !CALL PRTMAT(DWDC,3,3) 
          !PRINT *, "S_PIOLA = "
          !CALL PRTMAT(S_PIOLA,3,3) 
          !PRINT *, "SIGMA_CAUCHY = "
          !CALL PRTMAT(SIGMA_CAUCHY,3,3)
C
C
C     CALCULATE THE "SPATIAL" ELASTICITY TENSOR ::  C_SPATIAL_TAU
C
C
          DO 61 I = 1,3
	    DO 61 J = 1,3
	    DO 61 K = 1,3
	    DO 61 L = 1,3
	        QQ_TENSOR(I,J,K,L) = N_VEC(I,1)*N_VEC(K,1)*AIDENT(J,L) +
     +            N_VEC(I,1)*N_VEC(L,1)*AIDENT(J,K) +
     +            N_VEC(J,1)*N_VEC(K,1)*AIDENT(I,L) +
     +            N_VEC(J,1)*N_VEC(L,1)*AIDENT(I,K) 
          !PRINT *,"QQ_TENSOR(I,J,K,L)=",QQ_TENSOR(I,J,K,L)
61	    CONTINUE  
C
      	DO 62 I = 1,3
      	DO 62 J = 1,3
      	DO 62 K = 1,3
      	DO 62 L = 1,3
          C_SPATIAL_TAU(I,J,K,L) = 
     +        FOUR*D2W0DLAPLAS2*(LAPLAS_INVAR**TWO)*(ONE/J_INVAR)*
     +        N_VEC(I,1)*N_VEC(J,1)*N_VEC(K,1)*N_VEC(L,1) + 
     +        TWO*DW0DLAPLAS*LAPLAS_INVAR*(ONE/J_INVAR)*
     +        QQ_TENSOR(I,J,K,L) +
     +        EES*PHI_0*(L_INVAR**TWO)*(ONE/J_INVAR)*
     +        E_VEC(I,1)*E_VEC(J,1)*E_VEC(K,1)*E_VEC(L,1) +
     +        EK*(ONE/J_INVAR)*AIDENT(I,J)*AIDENT(K,L) +
     +        ( EG - EK*DLOG(J_INVAR) )*(ONE/J_INVAR)*
     +        ( AIDENT(I,K)*AIDENT(J,L) + AIDENT(I,L)*AIDENT(J,K) )
      !PRINT *,"C_SPATIAL_TAU(I,J,K,L)=",C_SPATIAL_TAU(I,J,K,L)
62    CONTINUE

C
C     CALCULATE THE "ABAQUS" ELASTICITY TENSOR : C_ABAQUS(I,J,K,L) = C_SPATIAL_TAU(I,J,K,L)+C_PRIME(I,J,K,L)
C          

C
C         CACLULATE C_PRIME(I,J,K,L)
C
      	DO 66 I = 1,3
      	DO 66 J = 1,3
      	DO 66 K = 1,3
      	DO 66 L = 1,3
          C_PRIME(I,J,K,L) = ONE_HALF*(
     +    AIDENT(I,K)*SIGMA_CAUCHY(J,L)+ AIDENT(I,L)*SIGMA_CAUCHY(J,K)+
     +    AIDENT(J,K)*SIGMA_CAUCHY(I,L)+ AIDENT(J,L)*SIGMA_CAUCHY(I,K) )
          !PRINT *,"C_PRIME(I,J,K,L)=",C_PRIME(I,J,K,L),I,J,K,L
66    	CONTINUE
C
C         CACLULATE C_ABAQUS(I,J,K,L) = C_SPATIAL_TAU(I,J,K,L)+C_PRIME(I,J,K,L)
C          
          DO 67 I = 1,3
	    DO 67 J = 1,3
	    DO 67 K = 1,3
	    DO 67 L = 1,3
	        C_ABAQUS(I,J,K,L) = C_SPATIAL_TAU(I,J,K,L) + C_PRIME(I,J,K,L) 
              !PRINT *,"C_ABAQUS(I,J,K,L)=",C_ABAQUS(I,J,K,L),I,J,K,L
67	    CONTINUE           

C
C     CALCULATE THE D2WDCDC TENSOR ::  (MATERIAL_MATRIX = 4*D2WDCDC )
C


      	DO 68 I = 1,3
      	DO 68 J = 1,3
      	DO 68 K = 1,3
      	DO 68 L = 1,3
          D2WDCDC(I,J,K,L) = 
     +        D2W0DLAPLAS2*DLAPLASDC(K,L)*DLAPLASDC(I,J)+
     +        DW0DLAPLAS*D2LAPLASDCDC(I,J,K,L)+
     +        ONE_FOURTH*EES*PHI_0*(ONE/L_INVAR**TWO)*
     +        E0(I,1)*E0(J,1)*E0(K,1)*E0(L,1) +
     +        ONE_HALF*EG*C_INV(I,K)*C_INV(L,J) +
     +        ONE_FOURTH*EK*C_INV(I,J)*C_INV(K,L) -
     +        ONE_HALF*EK*DLOG(J_INVAR)*C_INV(I,K)*C_INV(L,J)
          !PRINT *,"D2WDCDC(I,J,K,L)=",D2WDCDC(I,J,K,L)
68    CONTINUE      
			
      	DO I = 1,3
      	DO J = 1,3
      	DO K = 1,3
      	DO L = 1,3
      	    C_MAT_TAU(I,J,K,L) = ZERO
              C_MAT_TAU(I,J,K,L) = FOUR*D2WDCDC(I,J,K,L)
              !PRINT *,"C_MAT_TAU(I,J,K,L)=",C_MAT_TAU(I,J,K,L)
      	END DO
     	    END DO 
		END DO
		END DO
C
		DO I = 1,6
			STRESS_VECT(I) = ZERO
		END DO
		
		STRESS_VECT(1) = S_PIOLA(1,1)
      	STRESS_VECT(2) = S_PIOLA(2,2)
      	STRESS_VECT(3) = S_PIOLA(3,3)
      	STRESS_VECT(4) = S_PIOLA(1,2)
      	STRESS_VECT(5) = S_PIOLA(2,3)
      	STRESS_VECT(6) = S_PIOLA(3,1)
C
C%
C% 		CALCULATE THE Material stiffness
C%
        DO I=1,6
        DO J=1,6
              MATERIAL_MATRIX(I,J)=ZERO
        END DO
        END DO

      	MATERIAL_MATRIX(1,1) = C_MAT_TAU(1,1,1,1)
      	MATERIAL_MATRIX(1,2) = C_MAT_TAU(1,1,2,2)
      	MATERIAL_MATRIX(1,3) = C_MAT_TAU(1,1,3,3)
      	MATERIAL_MATRIX(1,4) = C_MAT_TAU(1,1,1,2)
      	MATERIAL_MATRIX(1,5) = C_MAT_TAU(1,1,2,3)
      	MATERIAL_MATRIX(1,6) = C_MAT_TAU(1,1,3,1)
      
      	MATERIAL_MATRIX(2,1) = C_MAT_TAU(2,2,1,1)
      	MATERIAL_MATRIX(2,2) = C_MAT_TAU(2,2,2,2)
      	MATERIAL_MATRIX(2,3) = C_MAT_TAU(2,2,3,3)
      	MATERIAL_MATRIX(2,4) = C_MAT_TAU(2,2,1,2)
      	MATERIAL_MATRIX(2,5) = C_MAT_TAU(2,2,2,3)
      	MATERIAL_MATRIX(2,6) = C_MAT_TAU(2,2,3,1)
      
      	MATERIAL_MATRIX(3,1) = C_MAT_TAU(3,3,1,1)
      	MATERIAL_MATRIX(3,2) = C_MAT_TAU(3,3,2,2)
      	MATERIAL_MATRIX(3,3) = C_MAT_TAU(3,3,3,3)
      	MATERIAL_MATRIX(3,4) = C_MAT_TAU(3,3,1,2)
      	MATERIAL_MATRIX(3,5) = C_MAT_TAU(3,3,2,3)
      	MATERIAL_MATRIX(3,6) = C_MAT_TAU(3,3,3,1)
      
      	MATERIAL_MATRIX(4,1) = C_MAT_TAU(1,2,1,1)
      	MATERIAL_MATRIX(4,2) = C_MAT_TAU(1,2,2,2)
      	MATERIAL_MATRIX(4,3) = C_MAT_TAU(1,2,3,3)
      	MATERIAL_MATRIX(4,4) = C_MAT_TAU(1,2,1,2)
      	MATERIAL_MATRIX(4,5) = C_MAT_TAU(1,2,2,3)
      	MATERIAL_MATRIX(4,6) = C_MAT_TAU(1,2,3,1)
      
      	MATERIAL_MATRIX(5,1) = C_MAT_TAU(2,3,1,1)
      	MATERIAL_MATRIX(5,2) = C_MAT_TAU(2,3,2,2)
      	MATERIAL_MATRIX(5,3) = C_MAT_TAU(2,3,3,3)
      	MATERIAL_MATRIX(5,4) = C_MAT_TAU(2,3,1,2)
      	MATERIAL_MATRIX(5,5) = C_MAT_TAU(2,3,2,3)
      	MATERIAL_MATRIX(5,6) = C_MAT_TAU(2,3,3,1)
      
      	MATERIAL_MATRIX(6,1) = C_MAT_TAU(3,1,1,1)
      	MATERIAL_MATRIX(6,2) = C_MAT_TAU(3,1,2,2)
      	MATERIAL_MATRIX(6,3) = C_MAT_TAU(3,1,3,3)
      	MATERIAL_MATRIX(6,4) = C_MAT_TAU(3,1,1,2)
      	MATERIAL_MATRIX(6,5) = C_MAT_TAU(3,1,2,3)
      	MATERIAL_MATRIX(6,6) = C_MAT_TAU(3,1,3,1)

          !PRINT *, "STRESS_VECT = "
          !CALL PRTMAT(STRESS_VECT,6,1)
          !PRINT *, "MATERIAL_MATRIX = "
          !CALL PRTMAT(MATERIAL_MATRIX,6,6)
C
C		CALCULATE D2WDKAPPADKAPPA TENSOR
C

		DO I=1,3
          DO J=1,3
          DO K=1,3
          DO L=1,3
        	    D2WDKAPPADKAPPA(I,J,K,L) = K3_COS*DDI6DKAPPADKAPPA(I,J,K,L)
          !PRINT *,"D2WDKAPPADKAPPA(I,J,K,L)=",D2WDKAPPADKAPPA(I,J,K,L)
          END DO
          END DO
          END DO
          END DO
        
C
C%
C% 		CALCULATE THE COSSERAT MATERIAL MATRIX
C%
        DO I=1,9
        DO J=1,9
              COSSERAT_MATRIX(I,J)=ZERO
        END DO
        END DO

      	COSSERAT_MATRIX(1,1) = D2WDKAPPADKAPPA(1,1,1,1)
      	COSSERAT_MATRIX(1,2) = D2WDKAPPADKAPPA(1,1,2,2)
      	COSSERAT_MATRIX(1,3) = D2WDKAPPADKAPPA(1,1,3,3)
      	COSSERAT_MATRIX(1,4) = D2WDKAPPADKAPPA(1,1,1,2)
      	COSSERAT_MATRIX(1,5) = D2WDKAPPADKAPPA(1,1,2,1)
      	COSSERAT_MATRIX(1,6) = D2WDKAPPADKAPPA(1,1,2,3)
      	COSSERAT_MATRIX(1,7) = D2WDKAPPADKAPPA(1,1,3,2)
      	COSSERAT_MATRIX(1,8) = D2WDKAPPADKAPPA(1,1,1,3)
      	COSSERAT_MATRIX(1,9) = D2WDKAPPADKAPPA(1,1,3,1)
      	
      	COSSERAT_MATRIX(2,1) = D2WDKAPPADKAPPA(2,2,1,1)
      	COSSERAT_MATRIX(2,2) = D2WDKAPPADKAPPA(2,2,2,2)
      	COSSERAT_MATRIX(2,3) = D2WDKAPPADKAPPA(2,2,3,3)
      	COSSERAT_MATRIX(2,4) = D2WDKAPPADKAPPA(2,2,1,2)
      	COSSERAT_MATRIX(2,5) = D2WDKAPPADKAPPA(2,2,2,1)
      	COSSERAT_MATRIX(2,6) = D2WDKAPPADKAPPA(2,2,2,3)
      	COSSERAT_MATRIX(2,7) = D2WDKAPPADKAPPA(2,2,3,2)
      	COSSERAT_MATRIX(2,8) = D2WDKAPPADKAPPA(2,2,1,3)
      	COSSERAT_MATRIX(2,9) = D2WDKAPPADKAPPA(2,2,3,1)
      
      	COSSERAT_MATRIX(3,1) = D2WDKAPPADKAPPA(3,3,1,1)
      	COSSERAT_MATRIX(3,2) = D2WDKAPPADKAPPA(3,3,2,2)
      	COSSERAT_MATRIX(3,3) = D2WDKAPPADKAPPA(3,3,3,3)
      	COSSERAT_MATRIX(3,4) = D2WDKAPPADKAPPA(3,3,1,2)
      	COSSERAT_MATRIX(3,5) = D2WDKAPPADKAPPA(3,3,2,1)
      	COSSERAT_MATRIX(3,6) = D2WDKAPPADKAPPA(3,3,2,3)
      	COSSERAT_MATRIX(3,7) = D2WDKAPPADKAPPA(3,3,3,2)
      	COSSERAT_MATRIX(3,8) = D2WDKAPPADKAPPA(3,3,1,3)
      	COSSERAT_MATRIX(3,9) = D2WDKAPPADKAPPA(3,3,3,1)
      
      	COSSERAT_MATRIX(4,1) = D2WDKAPPADKAPPA(1,2,1,1)
      	COSSERAT_MATRIX(4,2) = D2WDKAPPADKAPPA(1,2,2,2)
      	COSSERAT_MATRIX(4,3) = D2WDKAPPADKAPPA(1,2,3,3)
      	COSSERAT_MATRIX(4,4) = D2WDKAPPADKAPPA(1,2,1,2)
      	COSSERAT_MATRIX(4,5) = D2WDKAPPADKAPPA(1,2,2,1)
      	COSSERAT_MATRIX(4,6) = D2WDKAPPADKAPPA(1,2,2,3)
      	COSSERAT_MATRIX(4,7) = D2WDKAPPADKAPPA(1,2,3,2)
      	COSSERAT_MATRIX(4,8) = D2WDKAPPADKAPPA(1,2,1,3)
      	COSSERAT_MATRIX(4,9) = D2WDKAPPADKAPPA(1,2,3,1)
      	
      	COSSERAT_MATRIX(5,1) = D2WDKAPPADKAPPA(2,1,1,1)
      	COSSERAT_MATRIX(5,2) = D2WDKAPPADKAPPA(2,1,2,2)
      	COSSERAT_MATRIX(5,3) = D2WDKAPPADKAPPA(2,1,3,3)
      	COSSERAT_MATRIX(5,4) = D2WDKAPPADKAPPA(2,1,1,2)
      	COSSERAT_MATRIX(5,5) = D2WDKAPPADKAPPA(2,1,2,1)
      	COSSERAT_MATRIX(5,6) = D2WDKAPPADKAPPA(2,1,2,3)
      	COSSERAT_MATRIX(5,7) = D2WDKAPPADKAPPA(2,1,3,2)
      	COSSERAT_MATRIX(5,8) = D2WDKAPPADKAPPA(2,1,1,3)
      	COSSERAT_MATRIX(5,9) = D2WDKAPPADKAPPA(2,1,3,1)
      
      	COSSERAT_MATRIX(6,1) = D2WDKAPPADKAPPA(2,3,1,1)
      	COSSERAT_MATRIX(6,2) = D2WDKAPPADKAPPA(2,3,2,2)
      	COSSERAT_MATRIX(6,3) = D2WDKAPPADKAPPA(2,3,3,3)
      	COSSERAT_MATRIX(6,4) = D2WDKAPPADKAPPA(2,3,1,2)
      	COSSERAT_MATRIX(6,5) = D2WDKAPPADKAPPA(2,3,2,1)
      	COSSERAT_MATRIX(6,6) = D2WDKAPPADKAPPA(2,3,2,3)
      	COSSERAT_MATRIX(6,7) = D2WDKAPPADKAPPA(2,3,3,2)
      	COSSERAT_MATRIX(6,8) = D2WDKAPPADKAPPA(2,3,1,3)
      	COSSERAT_MATRIX(6,9) = D2WDKAPPADKAPPA(2,3,3,1)
      
      	COSSERAT_MATRIX(7,1) = D2WDKAPPADKAPPA(3,2,1,1)
      	COSSERAT_MATRIX(7,2) = D2WDKAPPADKAPPA(3,2,2,2)
      	COSSERAT_MATRIX(7,3) = D2WDKAPPADKAPPA(3,2,3,3)
      	COSSERAT_MATRIX(7,4) = D2WDKAPPADKAPPA(3,2,1,2)
      	COSSERAT_MATRIX(7,5) = D2WDKAPPADKAPPA(3,2,2,1)
      	COSSERAT_MATRIX(7,6) = D2WDKAPPADKAPPA(3,2,2,3)
      	COSSERAT_MATRIX(7,7) = D2WDKAPPADKAPPA(3,2,3,2)
      	COSSERAT_MATRIX(7,8) = D2WDKAPPADKAPPA(3,2,1,3)
      	COSSERAT_MATRIX(7,9) = D2WDKAPPADKAPPA(3,2,3,1)
      	
      	COSSERAT_MATRIX(8,1) = D2WDKAPPADKAPPA(1,3,1,1)
      	COSSERAT_MATRIX(8,2) = D2WDKAPPADKAPPA(1,3,2,2)
      	COSSERAT_MATRIX(8,3) = D2WDKAPPADKAPPA(1,3,3,3)
      	COSSERAT_MATRIX(8,4) = D2WDKAPPADKAPPA(1,3,1,2)
      	COSSERAT_MATRIX(8,5) = D2WDKAPPADKAPPA(1,3,2,1)
      	COSSERAT_MATRIX(8,6) = D2WDKAPPADKAPPA(1,3,2,3)
      	COSSERAT_MATRIX(8,7) = D2WDKAPPADKAPPA(1,3,3,2)
      	COSSERAT_MATRIX(8,8) = D2WDKAPPADKAPPA(1,3,1,3)
      	COSSERAT_MATRIX(8,9) = D2WDKAPPADKAPPA(1,3,3,1)
      
      	COSSERAT_MATRIX(9,1) = D2WDKAPPADKAPPA(3,1,1,1)
      	COSSERAT_MATRIX(9,2) = D2WDKAPPADKAPPA(3,1,2,2)
      	COSSERAT_MATRIX(9,3) = D2WDKAPPADKAPPA(3,1,3,3)
      	COSSERAT_MATRIX(9,4) = D2WDKAPPADKAPPA(3,1,1,2)
      	COSSERAT_MATRIX(9,5) = D2WDKAPPADKAPPA(3,1,2,1)
      	COSSERAT_MATRIX(9,6) = D2WDKAPPADKAPPA(3,1,2,3)
      	COSSERAT_MATRIX(9,7) = D2WDKAPPADKAPPA(3,1,3,2)
      	COSSERAT_MATRIX(9,8) = D2WDKAPPADKAPPA(3,1,1,3)
      	COSSERAT_MATRIX(9,9) = D2WDKAPPADKAPPA(3,1,3,1)
          
          !PRINT *, "COSSERAT_MATRIX = "
          !CALL PRTMAT(COSSERAT_MATRIX,9,9)
           
C
C***********************************************************************
C          End of loop over elements in an NBLOCK. Update the element
C          number pointer.
C***********************************************************************
C
          RETURN
          END
C
C
C          
C%*************************************************************************
C% Compute shape function, derivatives, and determinant of hexahedron element
C%*************************************************************************
C%%
        SUBROUTINE SHAPEL_HEX8(XI, ELXYZ, DETGJ, GJ, GDSF, SF)

        IMPLICIT REAL*8 (A-H,O-Z)
        
        INTEGER I, J, K
          
        REAL*8 XNODE(3,8), QUAR, SF(8,1), DSF(3,8), GDSF(3,8), XI(1,3),
     +     ELXYZ(8,3), XP, YP, ZP, XI0(1,3), GJ(3,3), GJINV(3,3), DETGJ
        
        PARAMETER (ZERO=0.D0 , ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +		  ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +		  THREE_HALF=3.D0/2.D0, THREE=3.D0, PI=4.D0*DATAN(1.D0),
     +          FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0 )
        
        XNODE(1,1) = -ONE
        XNODE(1,2) = ONE
        XNODE(1,3) = ONE
        XNODE(1,4) = -ONE
        XNODE(1,5) = -ONE
        XNODE(1,6) = ONE
        XNODE(1,7) = ONE
        XNODE(1,8) = -ONE
        
        XNODE(2,1) = -ONE
        XNODE(2,2) = -ONE
        XNODE(2,3) = ONE
        XNODE(2,4) = ONE
        XNODE(2,5) = -ONE
        XNODE(2,6) = -ONE
        XNODE(2,7) = ONE
        XNODE(2,8) = ONE
        
        XNODE(3,1) = -ONE
        XNODE(3,2) = -ONE
        XNODE(3,3) = -ONE
        XNODE(3,4) = -ONE
        XNODE(3,5) = ONE
        XNODE(3,6) = ONE
        XNODE(3,7) = ONE
        XNODE(3,8) = ONE
        
        QUAR = 0.125D0

          DO I = 1, 8
	        SF(I,1) = ZERO
	        DSF(1,I) = ZERO
	        DSF(2,I) = ZERO
	        DSF(3,I) = ZERO
          END DO
		
		DO I = 1, 8
			XP = XNODE(1,I)
    		    YP = XNODE(2,I)
    		    ZP = XNODE(3,I)
C
    		    XI0(1,1) = ONE+XI(1,1)*XP   	
    		    XI0(1,2) = ONE+XI(1,2)*YP		
    		    XI0(1,3) = ONE+XI(1,3)*ZP		
C
              SF(I,1) = QUAR*XI0(1,1)*XI0(1,2)*XI0(1,3) ! ---> THIS IS SHAPE FUNCTION AT EACH GAUSS POINTS
    		    DSF(1,I) = QUAR*XP*XI0(1,2)*XI0(1,3)	! ---> THESE ARE THE DERIVATIVE OF SHAPE FUNCTION W.R.T NETURAL COORDS AT EACH GAUSS POINTS
    		    DSF(2,I) = QUAR*YP*XI0(1,1)*XI0(1,3)
    		    DSF(3,I) = QUAR*ZP*XI0(1,1)*XI0(1,2)
          END DO
        
C      CALCULATE GLOBAL JACOUBIAN :=  GJ(3,3) = DSF(3,8)*ELXY(8,3) 
		DO I = 1,3
		DO J = 1,3
			GJ(I,J) = ZERO
			DO K = 1,8
				GJ(I,J) = GJ(I,J) + DSF(I,K)*ELXYZ(K,J)
			END DO
		END DO
		END DO
		
C		CALCULATE DETERMINATE GJ AND JACOUBIAN INVERSE :=  DETGJ = det(GJ) ;  GJINV=inv(GJ)
		CALL MATINV(GJ,GJINV,DETGJ)

C      CALCULATE GLOBAL DERIVATION OF SHAPE FUNCTIONS :=   GDSF(3,8) = GJINV(3,3)*DSF(3,8)
		DO I = 1,3
		DO J = 1,8
			GDSF(I,J) = ZERO
			DO K = 1,3
				GDSF(I,J) = GDSF(I,J) + GJINV(I,K)*DSF(K,J)
			END DO
		END DO
		END DO
C	
		RETURN
      	END

C****************************************************************************    
C    	Print the matrix A of dimension {m x n}
C**************************************************************************** 
      SUBROUTINE PRTMAT(A,M,N)

          IMPLICIT REAL*8 (A-H,O-Z)
          
          DIMENSION A(M,N)
    
          PRINT *, '*****************************'
    
          DO 10 I = 1,M
	        PRINT *, (A(I,J),J=1,N)
10        CONTINUE
    
          PRINT *, '*****************************'
    
      RETURN
      END

C****************************************************************************    
C    	Write the matrix A of dimension {m x n} in file number k
C****************************************************************************     
      SUBROUTINE WRTMAT(A,M,N,K)

	    IMPLICIT REAL*8 (A-H,O-Z)

	    DIMENSION A(M,N)
    
          WRITE(K,*) '*****************************'
    
          DO 11 I = 1,M
	        WRITE(K,*) (A(I,J),J=1,N)
11        CONTINUE
    
          WRITE(K,*) '*****************************'
    
      RETURN
      END   

C****************************************************************************    
C    	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 3 BY 3 MATRIX [A]
C****************************************************************************  
      SUBROUTINE ONEM(A)

          REAL*8 A(3,3)
 
          A(1,1) = 1.
          A(1,2) = 0.
          A(1,3) = 0.	 
          A(2,1) = 0.
          A(2,2) = 1.
          A(2,3) = 0.
          A(3,1) = 0.
          A(3,2) = 0.
          A(3,3) = 1.

      RETURN
      END

C****************************************************************************    
C    	THIS SUBROUTINE STORES THE ZERO MATRIX IN THE 3 BY 3 MATRIX [V]
C****************************************************************************  
      SUBROUTINE ZEROM(V)

          REAL*8 V(3,3)

          V(1,1) = 0.0
          V(1,2) = 0.0
          V(1,3) = 0.0
          V(2,1) = 0.0
          V(2,2) = 0.0
          V(2,3) = 0.0
          V(3,1) = 0.0
          V(3,2) = 0.0
          V(3,3) = 0.0
    
      RETURN
      END

C****************************************************************************    
C    	This subroutine calculates the inverse of a {3 x 3} matrix and the determinant of the inverse
C****************************************************************************  
      SUBROUTINE MATINV(A,A_INV,DET_A)

          IMPLICIT REAL*8 (A-H,O-Z)

          DIMENSION A(3,3), A_INV(3,3)
	
          PARAMETER(ZERO=0., ONE=1.)

	    OPEN(UNIT=80)

          DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) - 
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) + 
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3)) 

	    IF (DET_A .LE. ZERO) THEN
	        WRITE(80,*) 'WARNING: SUBROUTINE MATINV:'
	        WRITE(80,*) 'WARNING: DET of MAT is zero/negative!!'
          ENDIF

          DET_A_INV = ONE/DET_A

          A_INV(1,1) = DET_A_INV*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
	    A_INV(1,2) = DET_A_INV*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
	    A_INV(1,3) = DET_A_INV*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
	    A_INV(2,1) = DET_A_INV*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
	    A_INV(2,2) = DET_A_INV*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
	    A_INV(2,3) = DET_A_INV*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
	    A_INV(3,1) = DET_A_INV*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
	    A_INV(3,2) = DET_A_INV*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
	    A_INV(3,3) = DET_A_INV*(A(1,1)*A(2,2)-A(2,1)*A(1,2))

	RETURN
      END

C****************************************************************************    
C    	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C****************************************************************************  
	SUBROUTINE MPROD(A,B,C)

      	IMPLICIT REAL*8(A-H,O-Z)

          REAL*8 A(3,3), B(3,3), C(3,3)

          DO 2 I = 1, 3
	    DO 2 J = 1, 3
		    C(I,J) = 0.0
		    DO 1 K = 1, 3
			    C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1		    CONTINUE
2         CONTINUE
	
      RETURN
      END
      
C****************************************************************************    
C    	THIS SUBROUTINE MULTIPLIES TWO N BY M MATRICES [A] AND [B],	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C****************************************************************************  
	SUBROUTINE MPROD_GEN(A,ISIZE1,JSIZE1,B,ISIZE2,JSIZE2,C)

      	IMPLICIT REAL*8(A-H,O-Z)
      	

          REAL*8 A(ISIZE1,JSIZE1), B(ISIZE2,JSIZE2), C(ISIZE1,JSIZE2)

        DO I = 1, ISIZE1
	    DO J = 1, JSIZE2
		    C(I,J) = 0.0D0
		    DO K = 1, JSIZE1
			    C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
		    END DO
		END DO
		END DO
	
      RETURN
      END

C****************************************************************************    
C    	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 MATRIX [A], AND PLACES THE RESULT IN ATRANS.
C****************************************************************************  
	SUBROUTINE MTRANS(A,ATRANS)

      	IMPLICIT REAL*8(A-H,O-Z)

          REAL*8  A(3,3), ATRANS(3,3)

          CALL ONEM(ATRANS)
    	        DO 1 I = 1, 3
    	        DO 1 J = 1, 3
	    	    ATRANS(I,J) = A(J,I)
1             CONTINUE
	
      RETURN
	END

C****************************************************************************    
C    	THIS SUBROUTINE CALCULATES THE DETERMINANT OF A 3 BY 3 MATRIX [A].
C****************************************************************************
	SUBROUTINE MDET(A,DET)

          IMPLICIT REAL*8(A-H,O-Z)
    
          DIMENSION A(3,3)

          DET =	  A(1,1)*A(2,2)*A(3,3) 
     +        + A(1,2)*A(2,3)*A(3,1)  
     +        + A(1,3)*A(2,1)*A(3,2)  
     +        - A(3,1)*A(2,2)*A(1,3)  
     +        - A(3,2)*A(2,3)*A(1,1)  
     +        - A(3,3)*A(2,1)*A(1,2)  

	RETURN
      END		
			
		
!=========================== ABAQUS format user element subroutine ===================	
!
!       Variables that must be computed in this routine
!       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
!       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
!       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
!       ENERGY(1:8)
!                                  Energy(1) Kinetic Energy
!                                  Energy(2) Elastic Strain Energy
!                                  Energy(3) Creep Dissipation
!                                  Energy(4) Plastic Dissipation
!                                  Energy(5) Viscous Dissipation
!                                  Energy(6) Artificial strain energy
!                                  Energy(7) Electrostatic energy
!                                  Energy(8) Incremental work done by loads applied to the element
!       PNEWDT                     Allows user to control ABAQUS time increments.
!                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
!                                  a time increment equal to PNEWDT*DTIME
!                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
!
!       Variables provided for information
!       NDOFEL                     Total # DOF for the element
!       NRHS                       Dimension variable
!       NSVARS                     Total # element state variables
!       PROPS(1:NPROPS)            User-specified properties of the element
!       NPROPS                     No. properties
!       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
!       NJPROPS                    No. integer valued properties
!       COORDS(i,N)                ith coordinate of Nth node on element
!       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
!       U                          Vector of DOF at the end of the increment
!       DU                         Vector of DOF increments
!       V                          Vector of velocities (defined only for implicit dynamics)
!       A                          Vector of accelerations (defined only for implicit dynamics)
!       TIME(1:2)                  TIME(1)   Current value of step time
!                                  TIME(2)   Total time
!       DTIME                      Time increment
!       KSTEP                      Current step number (always 1 in EN234_FEA)
!       KINC                       Increment number
!       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
!       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
!       NDLOAD                     Number of user-defined distributed loads defined for this element
!       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
!       ADLMAG(1:NDLOAD)           Distributed load magnitudes
!       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
!       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
!       PREDEF(1,...)              Value of predefined field
!       PREDEF(2,...)              Increment in predefined field
!       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
!       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
!       NPREDF                     Number of predefined fields (1 for en234FEA)
!       LFLAGS                     Control variable
!       LFLAGS(1)                  Defines procedure type
!       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
!       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
!                                   2 => Subroutine must return stiffness AMATRX = -dF/du
!                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
!                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
!                                   5 => Define the RHS only
!                                   6 => Define the mass matrix for the initial acceleration calculation
!                                   100 => Define perturbation quantities for output
!       LFLAGS(4)                   0 => General step   1 => linear perturbation step
!       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
!       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
!       PERIOD                      Time period of the current step
!		
!
		