C     THIS IS THE CREEP MODEL FOR POWER LAW
C     ISOTROPIC IN 3D
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (NS_YU = 20)
C 
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C		
      DIMENSION YU_S(NTENS),DEL_CR_STRAIN(NTENS),YU_DELTA_STRESS(NTENS)
      DIMENSION YU_STRAIN(NTENS),YU_EN(NS_YU,NTENS,NTENS),YU_IN(NS_YU,NTENS,NTENS)
      DIMENSION XN(NTENS,NTENS), XMATD(NTENS,NTENS)
      DIMENSION YU_BETA_MU(NS_YU), YU_E_MU(NS_YU,NTENS,NTENS),YU_LAM_MU(NS_YU)
      DIMENSION YU_T_MU(NS_YU), YU_DELTA_Y_MU(NS_YU)
C	  
      DIMENSION TEMP_D(NTENS,NTENS), TEMP_E(NTENS,NTENS)
	  DIMENSION YU_E(NTENS,NTENS),YU_C(NTENS,NTENS)
	  DIMENSION XCOSL(3),XCOSM(3),XCOSN(3)
	  DIMENSION TRANSFORM_T(NTENS,NTENS)
C===========================================		
C     NS_YU STANDS FOR THE NUMBER OF KELVIN UNITS
C===========================================	
	  NPRNY = 15		
C===========================================					
C     IN THIS CODE, THE Z AXIS IS SET TO 
C        BE PERPENDICULAR TO THE BEDDING PLANE.	
C===========================================	
C	  XMATD STORED B_IJ	
C     ATTENTION! VOIGT NOTATION OF STRESS AND STRAIN
C       IS SLIGHTLY DIFFERENT FROM ABAQUS CONVENTIONAL 
C       ONE.
C===========================================
		XMATD = 0.D0
		XMATD(1,1) = PROPS(1)
		XMATD(2,2) = PROPS(2)
		XMATD(3,3) = PROPS(3)
		XMATD(4,4) = PROPS(6)
		XMATD(5,5) = PROPS(5)
		XMATD(6,6) = PROPS(4)
		XMATD(1,2) = PROPS(7)
		XMATD(1,3) = PROPS(8)
		XMATD(2,3) = PROPS(9)
		XMATD(2,1) = XMATD(1,2)
		XMATD(3,1) = XMATD(1,3)
		XMATD(3,2) = XMATD(2,3)
C===========================================		
C       XN STORED THE POWER INDICES N_IJ	
C===========================================
		XN = 0.D0
		XN(1,1) = PROPS(1+9)
		XN(2,2) = PROPS(2+9)
		XN(3,3) = PROPS(3+9)
		XN(4,4) = PROPS(6+9)
		XN(5,5) = PROPS(5+9)
		XN(6,6) = PROPS(4+9)
		XN(1,2) = PROPS(7+9)
		XN(1,3) = PROPS(8+9)
		XN(2,3) = PROPS(9+9)
		XN(2,1) = XN(1,2)
		XN(3,1) = XN(1,3)
		XN(3,2) = XN(2,3)
C===========================================		
		ROT = 0.D0
		PI = 4.D0*DATAN(1.D0) 
	    ROT = 0.D0*PI/6.D0	
C===========================================		
C		ROT, ANGLE FOR ROTATION,
C       	IS POSITIVE IF MEASURED FROM X TO X'
C 			IS CLOCKWISE.
C       IN THIS CASE, Y IS IDENTICAL TO Y'.
C===========================================
C===========================================		
C		CALCULATE TRANSFORMATION MATRIX
C===========================================
		XCOSL = 0.D0
		XCOSM = 0.D0
		XCOSN = 0.D0
C		
		XCOSL(1) = COS(ROT)
C
		XCOSL(2) = 0.D0
C		
		XCOSL(3) = -1.D0*SIN(ROT)
C		
		XCOSM(1) = 0.D0
C		
		XCOSM(2) = 1.D0
C		
		XCOSM(3) = 0.D0
C		
		XCOSN(1) = SIN(ROT)
C		
		XCOSN(2) = 0.D0
C		
		XCOSN(3) = COS(ROT)
C		
		TRANSFORM_T = 0.D0
		TRANSFORM_T(1,1)=XCOSL(1)**2
		TRANSFORM_T(2,1)=XCOSL(2)**2
		TRANSFORM_T(3,1)=XCOSL(3)**2
		TRANSFORM_T(4,1)=XCOSL(1)*XCOSL(2)
		TRANSFORM_T(5,1)=XCOSL(1)*XCOSL(3)
		TRANSFORM_T(6,1)=XCOSL(2)*XCOSL(3)
		TRANSFORM_T(1,2)=XCOSM(1)**2
		TRANSFORM_T(2,2)=XCOSM(2)**2
		TRANSFORM_T(3,2)=XCOSM(3)**2
		TRANSFORM_T(4,2)=XCOSM(1)*XCOSM(2)
		TRANSFORM_T(5,2)=XCOSM(1)*XCOSM(3)
		TRANSFORM_T(6,2)=XCOSM(2)*XCOSM(3)
		TRANSFORM_T(1,3)=XCOSN(1)**2
		TRANSFORM_T(2,3)=XCOSN(2)**2
		TRANSFORM_T(3,3)=XCOSN(3)**2
		TRANSFORM_T(4,3)=XCOSN(1)*XCOSN(2)
		TRANSFORM_T(5,3)=XCOSN(1)*XCOSN(3)
		TRANSFORM_T(6,3)=XCOSN(2)*XCOSN(3)
		TRANSFORM_T(1,4)=2.D0*XCOSL(1)*XCOSM(1)
		TRANSFORM_T(2,4)=2.D0*XCOSL(2)*XCOSM(2)
		TRANSFORM_T(3,4)=2.D0*XCOSL(3)*XCOSM(3)
		TRANSFORM_T(4,4)=XCOSL(1)*XCOSM(2)+XCOSL(2)*XCOSM(1)
		TRANSFORM_T(5,4)=XCOSL(1)*XCOSM(3)+XCOSL(3)*XCOSM(1)
		TRANSFORM_T(6,4)=XCOSL(2)*XCOSM(3)+XCOSL(3)*XCOSM(2)
		TRANSFORM_T(1,5)=2.D0*XCOSM(1)*XCOSN(1)
		TRANSFORM_T(2,5)=2.D0*XCOSM(2)*XCOSN(2)
		TRANSFORM_T(3,5)=2.D0*XCOSM(3)*XCOSN(3)
		TRANSFORM_T(4,5)=XCOSM(1)*XCOSN(2)+XCOSM(2)*XCOSN(1)
		TRANSFORM_T(5,5)=XCOSM(1)*XCOSN(3)+XCOSM(3)*XCOSN(1)
		TRANSFORM_T(6,5)=XCOSM(2)*XCOSN(3)+XCOSM(3)*XCOSN(2)
		TRANSFORM_T(1,6)=2.D0*XCOSL(1)*XCOSN(1)
		TRANSFORM_T(2,6)=2.D0*XCOSL(2)*XCOSN(2)
		TRANSFORM_T(3,6)=2.D0*XCOSL(3)*XCOSN(3)
		TRANSFORM_T(4,6)=XCOSL(1)*XCOSN(2)+XCOSL(2)*XCOSN(1)
		TRANSFORM_T(5,6)=XCOSL(1)*XCOSN(3)+XCOSL(3)*XCOSN(1)
		TRANSFORM_T(6,6)=XCOSL(2)*XCOSN(3)+XCOSL(3)*XCOSN(2)
C		
C===========================================
C      ASSIGN INITIAL VALUES TO STATEV
C===========================================		
		IF (TIME(2) .LE. 0) THEN
C
		DO II=1, NSTATV
			STATEV(II)=0.D0
		END DO
C
        END IF
C===========================================		
C       GIVE VALUE TO RETARDATION TIME,
C       INCREASING BY 10
C       YU_T_MU STORED THE RETARDATION TIMES TAU_MU
C       YU_DELTA_Y_MU STORED DELTA_T/TAU_MU
C===========================================
		YU_T_MU(1) = 1.0D-6
        DO I = 2, NS_YU
           YU_T_MU(I) = YU_T_MU(I-1)*10.0D0
        END DO
        DO I = 1, NS_YU
	      YU_DELTA_Y_MU(I) = DTIME/YU_T_MU(I)
        END DO
C===========================================		
C     	NPRNY STANDS FOR K IN THE KTH ORDER APPROXIMATION 
C		FCT STANDS FOR FACTORIAL(K-1)
C       DVL STANDS FOR GAMMANK,N*(N-1)*(N-2)*...*(N-K+1)
C===========================================		
		FCT = 1.D0
		DO I = 1, NPRNY-1
			FCT=FCT*I
		END DO
C===========================================
C       CALCULATE BETA_MU AND LAMBDA_MU
C===========================================

        DO I = 1, NS_YU
	      TEMP_YU = EXP( 0.0D0-YU_DELTA_Y_MU(I) )
		  YU_BETA_MU(I)=TEMP_YU
	      YU_LAM_MU(I) = (1.0D0 - TEMP_YU)/YU_DELTA_Y_MU(I)
        END DO
C===========================================		
C		READING THE STRESS HISTORY TERM FROM 
C			THE LAST STEP, STORE IT IN YU_S		
C===========================================		
		DO K1=1,NTENS
	      YU_S(K1) = STRESS(K1)
        END DO 
C		
C		
		YU_E = 0.D0
		YU_C = 0.D0
		DDSDDE = 0.D0
		YU_E_MU= 0.D0
		TEMP_D = 0.D0
		TEMP_E = 0.D0
C===========================================
C     IN THE FOLLOWING LOOPS, 
C        THE INDICES NEED TO BE SPECIFIED AS
C        :KJ1=K, KJ2=L, KJ3=I, KJ4=J 		
C===========================================
        DO KJ1 = 1,NTENS
           DO KJ2 = 1,NTENS
C===========================================			
C      YU_E: THE INVERSE OF EMU, COMPLIANCE MATRIX.   
C	   YU_C: 1/EMU_IJ IN EQ.(16)	
C===========================================
				YU_E_ZERO = 0.D0 
C				
                DO I=1, NS_YU
                    DO KJ3 = 1, NTENS
						DO KJ4 = 1, NTENS
C						
						DVL = 1.D0
						DO II = 1, NPRNY-1
							DVL=DVL*(XN(KJ3,KJ4)-II*1.D0)
						END DO	
C						
					    IF (XMATD(KJ3,KJ4) .NE. 0.D0) THEN
C
C 						
				 YU_E_MU(I,KJ1,KJ2)=YU_E_MU(I,KJ1,KJ2)+TRANSFORM_T(KJ3,KJ1)*TRANSFORM_T(KJ4,KJ2)
     &			 *XMATD(KJ3,KJ4)*(LOG(10.D0)
     &			 *(-1.D0)**(NPRNY+1)/FCT*NPRNY**(XN(KJ3,KJ4))*DVL*(XN(KJ3,KJ4))
     &           *YU_T_MU(I)**XN(KJ3,KJ4))
C	 
                    IF (I .EQ. 1) THEN
					XQ = 0.D0
				    XQ = (-1.D0)**(NPRNY+1)/FCT*NPRNY**(XN(KJ3,KJ4))*DVL
     &			*(YU_T_MU(1)/SQRT(10.D0))**XN(KJ3,KJ4)
C					
					YU_E_ZERO = YU_E_ZERO + TRANSFORM_T(KJ3,KJ1)*TRANSFORM_T(KJ4,KJ2)
     &				*XMATD(KJ3,KJ4)*XQ
C
                    END IF
                    END IF
C	 
                        ENDDO
                    ENDDO
C				
					YU_E(KJ1,KJ2) = YU_E(KJ1,KJ2) + ( 1.0D0-YU_LAM_MU(I) )*YU_E_MU(I,KJ1,KJ2)
					YU_C(KJ1,KJ2) = YU_C(KJ1,KJ2) + ( 1.0D0-YU_BETA_MU(I) )*YU_E_MU(I,KJ1,KJ2)
C				
				END DO
C				
				TEMP_D(KJ1,KJ2) = YU_E_ZERO + YU_E(KJ1,KJ2)
C
			END DO
		END DO
C===========================================
C       EQ.(17B)
C===========================================		
        CALL INVERT_MATRIX(TEMP_D, TEMP_E)
C		
		YU_E = TEMP_E
		DDSDDE = YU_E
C===========================================
C      READ INELASTIC CREEP STRAIN AT LAST STEP,
C         STORE THEM IN YU_IN
C      IN TOTAL NSTATEV=NS_YU*6*6=720,
C         1st KELVIN UNIT HAS 1~6.
C===========================================
        K1 = 1
        DO KMU = 1, NS_YU
           DO I = 1, 6
			   DO J = 1,6 
					YU_IN(KMU,I,J) = STATEV(K1)
					K1 = K1 + 1
			   END DO
           ENDDO
        ENDDO		
C===========================================		
C    CALCULATE THE INELASTIC STRAIN INCREMENT,
C       EQUATION(16)
C===========================================
        DO I = 1,NTENS
	      YU_STRAIN(I) = 0.D0
		  DEL_CR_STRAIN(I) = 0.D0
           DO J = 1, NTENS
			   DEL_CR_STRAIN(I)=DEL_CR_STRAIN(I)+YU_C(I,J)*YU_S(J)
			   DO KMU = 1,NS_YU
					DEL_CR_STRAIN(I) = DEL_CR_STRAIN(I) 
     +		           + (-1.D0) * YU_IN(KMU,I,J) * (1.D0-YU_BETA_MU(KMU))
			   END DO
           END DO
        END DO 	
C===========================================		
C      ELASTIC STRAIN INCREMENT
C      DELTA EPSILON = DELTA TOTAL EPSILON - DELTA EPSILON''
C===========================================
        DO I = 1, NTENS
	      YU_STRAIN(I) = DSTRAN(I) - DEL_CR_STRAIN(I)
        END DO
         YU_DELTA_STRESS = 0.0D0
C===========================================
C		CALCULATE THE STRERSS INCREMENT
C===========================================
         DO K1=1,NTENS
             DO K2=1,NTENS
	          YU_DELTA_STRESS(K1) = YU_DELTA_STRESS(K1) + 
     +		                DDSDDE(K1,K2)*YU_STRAIN(K2)
             END DO
          STRESS(K1) = STRESS(K1) + YU_DELTA_STRESS(K1)
        END DO 
C===========================================
C      CALCULATE PARTIAL STRAINS ACCORDING TO EQ(15)
C===========================================
		YU_EN = 0.0D0 
        DO KMU = 1, NS_YU
           DO I = 1, NTENS
			   DO J = 1, NTENS
				   YU_EN(KMU,I,J)=YU_IN(KMU,I,J)*YU_BETA_MU(KMU)+(1.D0-YU_BETA_MU(KMU))
     &              *YU_E_MU(KMU,I,J)*YU_S(J)+(1.D0-YU_LAM_MU(KMU))*YU_E_MU(KMU,I,J)
     &              *YU_DELTA_STRESS(J)
			   END DO
           END DO
        END DO		
C===========================================	   
C       STORE THE UPDATED GAMA FOR NEXT STEP
C===========================================
        K1 = 1
		DO KMU = 1, NS_YU
           DO I = 1, 6
			   DO J = 1,6 
					STATEV(K1) = YU_EN(KMU,I,J)
						K1 = K1 + 1
			   END DO
           ENDDO
        ENDDO	
C		
C		
      RETURN
      END
C
C     
C	  	  
      SUBROUTINE INVERT_MATRIX(MATRIX, INV_MATRIX)
C     Inverts a square matrix using Gauss-Jordan elimination
C     Inputs:
C       MATRIX - Input matrix to be inverted
C     Outputs:
C       INV_MATRIX - Inverted matrix
C=================================================
C       Size of the matrix
C=================================================
      PARAMETER (N = 6)     
C=================================================
C     Input matrix
C=================================================
      REAL*8 MATRIX(N, N)  
C=================================================
C	  Inverted matrix
C=================================================
      REAL*8 INV_MATRIX(N, N)  
C=================================================
      INTEGER I, J, K
      REAL*8 TEMP
C=================================================
C     Create the augmented matrix [MATRIX | IDENTITY]
C=================================================
      DO I = 1, N
         DO J = 1, N
            IF (J == I) THEN
               INV_MATRIX(I, J) = 1.D0
            ELSE
               INV_MATRIX(I, J) = 0.D0
            ENDIF
         ENDDO
      ENDDO
C=================================================
C     Perform Gauss-Jordan elimination
C=================================================
      DO I = 1, N
C=================================================
C        Divide row I by the pivot element
C=================================================
         TEMP = MATRIX(I, I)
         DO JJ = 1, N
            MATRIX(I, JJ) = MATRIX(I, JJ) / TEMP
            INV_MATRIX(I, JJ) = INV_MATRIX(I, JJ) / TEMP
         ENDDO
C=================================================
C        Eliminate other elements in the column
C=================================================
       DO K = 1, N
          IF (K /= I) THEN
C		  
             TEMP= MATRIX(K,I)
             DO J = 1, N
C               
                MATRIX(K, J) = MATRIX(K, J) - TEMP * MATRIX(I,J)
                INV_MATRIX(K,J)= INV_MATRIX(K, J)-TEMP * INV_MATRIX(I,J)
             ENDDO
          ENDIF
       ENDDO
      ENDDO
C
      RETURN
      END