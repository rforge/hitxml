*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*                                                                      *
      SUBROUTINE MGDRAW ( ICODE, MREG, NEWREG)
*      
*     main commons
*      
*     defined variables
      INCLUDE '(DBLPRC)'
*     array dimensions
      INCLUDE '(DIMPAR)'
*     logical in- and output numbers
      INCLUDE '(IOUNIT)'
*
*     standard commons
      INCLUDE '(CASLIM)'
      INCLUDE '(COMPUT)'
      INCLUDE '(SOURCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(MGDDCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(QUEMGD)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'
      INCLUDE '(FLKMAT)'
      INCLUDE '(USRBIN)'
      INCLUDE '(SCOHLP)'

*
      DIMENSION DTQUEN ( MXTRCK, MAXQMG )
*	  
      CHARACTER*20 FILNAM
*	  necessary for region names
      CHARACTER*8 MRGNAM, NRGNAM
*         kinetic energy, weight in amu
      DOUBLE PRECISION EKIN, WAMU, EKINAMU
*         flag for first call
      LOGICAL FIRSTCAL
      SAVE FIRSTCAL
      DATA FIRSTCAL / .TRUE. /
      
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /
  

      LOGICAL LFIRST
      SAVE LFIRST
      DATA LFIRST /.TRUE./
*
      DOUBLE PRECISION MY_E, MY_P, MY_M, LET_DELTA, MY_DELTA, LET_UNRES
      INTEGER MY_MAT
      RETURN
*
*======================================================================*
*     Variables for BXDRAW
*
*     STEP   : number of routine call
*     ICODE  : Fluka physical compartment originating the call
*              = 1: call from subroutine KASKAD (hadrons and muons)
*              = 2: call from subroutine EMFSCO (e − , e + and photons)
*              = 3: call from subroutine KASNEU (low-energy neutrons)
*              = 4: call from subroutine KASHEA (heavy ions)
*              = 5: call from subroutine KASOPH (optical photons)
*     LTRACK : unique particle identifier
*     JTRACK : identity number of the particle 
*     PRNAME(i) : alphabetical name of the i_th particle      
*     MREG   : current geometry region
*     NEWREG : region the particle is entering
*                 
*     ETRACK : total energy of the particle
*     AM(i)  : mass energy of the i_th particle
*     PTRACK : momentum of the particle
*     ICHRGE(i) : electric charge of the i_th particle
*
*     XSCO, YSCO, ZSCO : point where the boundary crossing occurs
*     CXTRCK, CYTRCK, CZTRCK : direction cosines of the current
*                              particle 
*======================================================================*
*      Material index for water - change if needed

*      Write message the first time the routine is called
 
*         WRITE(LUNOUT,*) '****************************************'
*         WRITE(LUNOUT,*) ' '
         

  
      ENTRY BXDRAW  ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )
*
*     check if call from ion compartment
      IF( ICODE .NE. 19 .AND. ICODE .NE. 49) THEN
         RETURN
      ENDIF

*     get corresponding region name for region number
      CALL GEOR2N ( MREG, MRGNAM, IERR1 )
      CALL GEOR2N ( NEWREG, NRGNAM, IERR2 )
      IF(IERR1 .NE. 0 .OR. IERR2 .NE. 0) STOP "Error in name conversion"

*     check if at selected BOUNDARY
      if( MRGNAM(1:6) .NE. "TARGET" .OR. NRGNAM(1:6) .NE. "TARGET") THEN
         RETURN
      ENDIF

*     get kinetic energy per amu
      EKIN    = ETRACK - AM(JTRACK)
      WAMU    = AM(JTRACK) / AmugeV
      IF ( WAMU .GT. ZERZER) THEN
          EKINAMU = EKIN / WAMU * GEVMEV
      END IF
*  
*     Getting LET
      MY_MAT = 27
      LET_DELTA = ZERZER
      LET_UNRES= ZERZER
      MY_DELTA= 1E-3
      MY_M= AM(JTRACK)
      MY_P = PTRACK
      MY_E = SQRT(PTRACK**2+MY_M**2)-MY_M
      IF ((JTRACK.GT.-7).AND.(JTRACK.LE.2))THEN
*        WRITE(LUNOUT,*) 'Shirin lets see  '
         LET_DELTA =GETLET(JTRACK, MY_E, MY_P,MY_DELTA, MY_MAT)
         LET_UNRES= GETLET(JTRACK, MY_E, MY_P,-1, MY_MAT)
*        WRITE(LUNOUT,*) 'ME_P  ' , MY_P
*        WRITE(LUNOUT,*) 'MY_E', MY_E
*        WRITE(LUNOUT,*) 'EKIN', EKIN
*        WRITE(LUNOUT,*) 'LET_DELTA', LET_DELTA
      ENDIF	
*     readout files
*      
*     output file

*     output file
*
      OPEN ( UNIT = 99, FILE = "phase_space", STATUS ='UNKNOWN')
*     header 
      IF ( FIRSTCAL ) THEN
         WRITE(99,*) "ICODE OLDREG NEWREG JTRACK PRNAME XSCO YSCO  
     & ZSCO CXTRCK CYTRCK CZTRCK EKIN ICHRGE WAMU LET.UNR"
         FIRSTCAL = .FALSE.
      END IF
*     writing	
999   FORMAT ( I7, 2(A10), I7, A10, 7(F10.4), I7 , 2(F10.4))
      WRITE(99,999) ICODE, MRGNAM, NRGNAM, JTRACK, PRNAME(JTRACK),
     &     XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK,
     &     EKINAMU, ICHRGE(JTRACK), WAMU, LET_UNRES

      RETURN
*


*
*======================================================================*
*                                                                      *
      ENTRY EEDRAW ( ICODE )
      RETURN
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO )
      RETURN
*
*======================================================================*
*                                                                      *
      ENTRY SODRAW
      RETURN
*
*======================================================================*
*                                                                      *
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      RETURN
*=== End of subroutine Mgdraw ==========================================*
      END
