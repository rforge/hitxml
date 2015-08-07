*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*
*=== modified to writeout dEdx                                      ===*
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

*     kinetic energy, weight in amu, kinetic energy per amu
      DOUBLE PRECISION EKIN, WAMU, EKINAMU
      DOUBLE PRECISION MY_E, MY_P, MY_M, LET_DELTA, MY_DELTA, LET
      INTEGER MY_MAT

      INTEGER MY_IJS(3)
      DATA MY_IJS/1, -6, -2/

      DOUBLE PRECISION MY_ES(11)
      DATA MY_ES/1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0/

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

*     Material index for water - change if needed (check log file to
*     see material indices)
      MY_MAT     = 27

*     Writeout dEdx tables during first call
      IF( LFIRST ) THEN
       OPEN ( UNIT = 95, FILE = "dEdx", STATUS ='UNKNOWN')
       WRITE(95,*) "PRNAME Z EKIN.MEV.U LET.UNR.KEV.UM"

        DO 998 I = 1,3
*        Current particle id
         MY_IJ   = MY_IJS(I)

         DO 997 J = 1,6
          DO 996 K = 1,11

*         Current kinetic energy (MeV/amu)
          MY_E = 10.0**(J-3) * MY_ES(K)

          MY_M = AM(MY_IJ)
          WAMU = MY_M / AmugeV
          EKIN = MY_E * WAMU / GEVMEV
          MY_P = SQRT( (EKIN + MY_M) * (EKIN + MY_M) - MY_M**2)
          LET  = GETLET(MY_IJ, EKIN, MY_P, -1, MY_MAT) 
995       FORMAT ( A10, I5, 2(F10.4))
          WRITE(95,995) PRNAME(MY_IJ), ICHRGE(MY_IJ), MY_E, LET

996       CONTINUE
997     CONTINUE
998    CONTINUE
      END IF

      LFIRST = .FALSE.

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
