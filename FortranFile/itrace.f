
SUBROUTINE ITRACE (IAPX)
C          Follow a geomagnetic field line until passing its apex
C
C          INPUTS:
C            (all are in common blocks)
C          OUTPUTS:
C            IAPX = 2 (when apex passed) or 1 (not)
C
C          This uses the 4-point Adams formula after initialization.
C          First 7 iterations advance point by 3 steps.
C
C          COMMON BLOCKS:
C            COMMON /APXIN/   YAPX(3,3)
C            COMMON /FLDCOMD/ BX, BY, BZ, BB
C            COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
C
C          APXIN has step locations determined in ITRACE:
C            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
C                    three points about the apex.  Set in subroutine ITRACE.
C
C          FLDCOMD has geomagnetic field at current trace point:
C            BX    = X component (Gauss)
C            BY    = Y component (Gauss)
C            BZ    = Z component (Gauss)
C            BB    = Magnitude   (Gauss)
C
C          ITRA has field line tracing variables determined in LINAPX:
C            NSTP  = Step count.
C            Y     = Array containing current tracing point cartesian coordinates.
C            YOLD  = Array containing previous tracing point cartesian coordinates.
C            SGN   = Determines direction of trace.
C            DS    = Step size (arc length in km).
C
C          REFERENCES:
C            Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
C            Greenbelt, Maryland
C------------------------------------------------------------------------------
C          HISTORY:
C          Oct 1973: Initial version completed on the 29th by W. Clark, NOAA ERL
C                    Laboratory.
C          Feb 1988: Revised by H. Passi, NCAR.
C          Apr 2004: Replace computed GO TO with IF blocks because some compilers
C                    are threatening to remove this old feature
C
      COMMON /APXIN/   YAPX(3,3)
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
      DIMENSION YP(3,4)
      SAVE
C          Statement function
      RDUS(D,E,F) = SQRT (D**2 + E**2 + F**2)
      IAPX = 1
C          Cartesian component magnetic field (partial) derivitives steer the trace
      YP(1,4) = SGN*BX/BB
      YP(2,4) = SGN*BY/BB
      YP(3,4) = SGN*BZ/BB
      IF (NSTP .LE. 7) THEN
	DO 10 I=1,3
	IF (NSTP .EQ. 1) THEN
	  D2        = DS/2.
	  D6        = DS/6.
	  D12       = DS/12.
	  D24       = DS/24.
	  YP(I,1)   = YP(I,4)
	  YOLD(I)   = Y(I)
	  YAPX(I,1) = Y(I)
	  Y(I)      = YOLD(I) + DS*YP(I,1)
	ELSE IF (NSTP .EQ. 2) THEN
	  YP(I,2) = YP(I,4)
	  Y(I)    = YOLD(I) + D2*(YP(I,2)+YP(I,1))
	ELSE IF (NSTP .EQ. 3) THEN
	  Y(I) = YOLD(I) + D6*(2.*YP(I,4)+YP(I,2)+3.*YP(I,1))
	ELSE IF (NSTP .EQ. 4) THEN
	  YP(I,2)   = YP(I,4)
	  YAPX(I,2) = Y(I)
	  YOLD(I)   = Y(I)
	  Y(I)      = YOLD(I) + D2*(3.*YP(I,2)-YP(I,1))
	ELSE IF (NSTP .EQ. 5) THEN
	  Y(I) = YOLD(I) + D12*(5.*YP(I,4)+8.*YP(I,2)-YP(I,1))
	ELSE IF (NSTP .EQ. 6) THEN
	  YP(I,3)   = YP(I,4)
	  YOLD(I)   = Y(I)
	  YAPX(I,3) = Y(I)
	  Y(I)      = YOLD(I) + D12*(23.*YP(I,3)-16.*YP(I,2)+5.*YP(I,1))
	ELSE IF (NSTP .EQ. 7) THEN
	  YAPX(I,1) = YAPX(I, 2)
	  YAPX(I,2) = YAPX(I, 3)
	  Y(I)      = YOLD(I) + D24*(9.*YP(I,4) + 19.*YP(I,3) -
     +                               5.*YP(I,2) +     YP(I,1))
	  YAPX(I,3) = Y(I)
	ENDIF
   10   CONTINUE
	IF (NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN        ! signal if apex passed
	  RC = RDUS (YAPX(1,3), YAPX(2,3), YAPX(3,3))
	  RP = RDUS (YAPX(1,2), YAPX(2,2), YAPX(3,2))
	  IF (RC .LT. RP) IAPX = 2
	ENDIF
      ELSE                 ! NSTP > 7
	DO 30 I=1,3
	YAPX(I,1) = YAPX(I,2)
	YAPX(I,2) = Y(I)
	YOLD(I)   = Y(I)
	Y(I)      = YOLD(I) + D24*(55.*YP(I,4) - 59.*YP(I,3) +
     +                             37.*YP(I,2) -  9.*YP(I,1))
	YAPX(I,3) = Y(I)
	DO 20 J=1,3
   20   YP(I,J) = YP(I,J+1)
   30   CONTINUE
	RC = RDUS (   Y(1),    Y(2),    Y(3))
	RP = RDUS (YOLD(1), YOLD(2), YOLD(3))
	IF (RC .LT. RP) IAPX = 2
      ENDIF
      RETURN
      END