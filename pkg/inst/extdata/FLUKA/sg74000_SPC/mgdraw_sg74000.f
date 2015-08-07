*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*
*=== modified to score phase space at a number of boundaries within ===*
*=== a target. 2015-07, DKFZ E040-8, SG, SR, AN                     ===*
*
* Uses GETLET function to obtain the LET in a material
* Syntax: 
*
*   DOUBLE PRECISION FUNCTION GETLET ( IJ, EKIN, PLA, TDELTA, MATLET ) 
*
* with: 
* IJ = particle index 
* EKIN = particle kinetic energy (GeV) 
* PLA = particle momentum (GeV/c) 
* TDELTA = maximum secondary electron energy (GeV) (unrestricted if =< 0) 
* MATLET = material index for which LET is requested 
*
*


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
*     region names
      CHARACTER*8 MRGNAM, NRGNAM

*     Variables for phase space scoring
      DOUBLE PRECISION EKIN, WAMU
      DOUBLE PRECISION MY_E, MY_P, MY_M, LET_UNRES
      INTEGER MY_MAT

      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /
  
*     flag for first call
      LOGICAL LFIRST
      SAVE LFIRST
      DATA LFIRST /.TRUE./

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

      ENTRY BXDRAW  ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )


*     check if call from ion compartment
      IF( ICODE .NE. 19 .AND. ICODE .NE. 49) THEN
         RETURN
      ENDIF

*     get corresponding region name for region number
      CALL GEOR2N ( MREG, MRGNAM, IERR1 )
      CALL GEOR2N ( NEWREG, NRGNAM, IERR2 )
      IF(IERR1 .NE. 0 .OR. IERR2 .NE. 0) STOP "Error in name conversion"

*     check whether particle passes selected boundary
      if( MRGNAM(1:6) .NE. "TARGET" .OR. NRGNAM(1:6) .NE. "TARGET") THEN
         RETURN
      ENDIF

*     Message the first time the routine is called
      IF( LFIRST ) THEN
         WRITE(LUNOUT,*) '****************************************'
         WRITE(LUNOUT,*) 'Customized MGDRAW routine sg74000 called'
         WRITE(LUNOUT,*) '****************************************'
      ENDIF

*     Material index for water - change if needed (check log file to
*     see material indices)
      MY_MAT     = 27

*     Get kinetic energy per amu
      EKIN    = ETRACK - AM(JTRACK)
      WAMU    = AM(JTRACK) / AmugeV
      IF ( WAMU .GT. ZERZER) THEN
          EKINAMU = EKIN / WAMU * GEVMEV
      END IF

*     Get LET
      LET_UNRES  = ZERZER
      MY_M       = AM(JTRACK)
      MY_P       = PTRACK
      MY_E       = SQRT(PTRACK**2+MY_M**2)-MY_M
      LET_UNRES  = GETLET(JTRACK, MY_E, MY_P,-1, MY_MAT)
      
*     output to phase space file
      OPEN ( UNIT = 99, FILE = "phase_space", STATUS ='UNKNOWN')

*     header (first call)
      IF ( LFIRST ) THEN
         WRITE(99,*) "ICODE OLDREG NEWREG JTRACK PRNAME XSCO YSCO ZSCO 
     &CXTRCK CYTRCK CZTRCK EKIN ICHRGE WAMU LET.UNR"
         LFIRST = .FALSE.
      END IF

*     writing data
999   FORMAT ( I7, 2(A10), I7, A10, 7(F10.4), I7 , 2(F10.4))
      WRITE(99,999) ICODE, MRGNAM, NRGNAM, JTRACK, PRNAME(JTRACK),
     &     XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK,
     &     EKINAMU, ICHRGE(JTRACK), WAMU, LET_UNRES

      RETURN

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
