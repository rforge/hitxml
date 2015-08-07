* MODIFIED SOURCE.F USER ROUTINE TO SAMPLE FROM A ENERGY DISTRIBUTION
* Eduardo G. Yukihara
* modified S Rahmanian
* 10 July 2015

      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'

      LOGICAL LFIRST
      SAVE LFIRST
      DATA LFIRST / .TRUE. /

* The following variables are need to sample from the
* available distribution of particles per position.
* Current limit is 100,000 beamlets (total positions).
* If needed, this needs to be increased.

* USREN = energy of the beam (GeV/u)
* USRFOC = focus size of the beam (cm)
* USRX = x position of the beam (cm)
* USRY = y position of the beam (cm)
* USRMIN = beginning of particle interval (included)
* USRMAX = end of particle interval (excluded)
* note that USRMAX(I-1) = USRMIN(I)
* NPOINTS = number of positions (beamlets)
* MAX_EN = maximum energy in the sampled plan
* TOTALP = total number of particles
* USRNPART = number of particle in particular spot

      DOUBLE PRECISION USREN(100000), USRFOC(100000),
     &     USRX(100000), USRY(100000), USRMIN(100000), USRMAX(100000)
      DOUBLE PRECISION USRNPART, TOTALP
      INTEGER NPOINTS,  I
      DOUBLE PRECISION RND, XPOS, YPOS, EN, FOC
      DOUBLE PRECISION MAX_EN
      
      NOMORE = 0

* ---------------------------------------------------------------------
* First call initializations:

      IF ( LFIRST ) THEN
         WRITE(LUNOUT,*) ' ' 
         WRITE(LUNOUT,*) 'Calling SOURCE.F user-routine by yukihara'
         WRITE(LUNOUT,*) 'version from 24 June 2015'
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.

* This will read the file to get the positions and particle distributions.
* Data format should be text with space separated number in the sequence:
* energy (GeV/u), focus (cm), x(cm), y(cm), 

* Initializes values with zeros
         NPOINTS = 0
         MAX_EN = ZERZER

* read file to be sampled (OPEN card, 98 ASC, status = OLD)
         DO I = 1, 100000
            READ(98,*, END = 999) USREN(I), USRFOC(I),
     &      USRX(I), USRY(I), USRNPART
            USRMIN(I) = TOTALP
            USRMAX(I) = TOTALP + USRNPART
            TOTALP = TOTALP + USRNPART
            IF (USREN(I) .GT. MAX_EN) THEN 
               MAX_EN = USREN(I) 
            END IF
            NPOINTS = NPOINTS + 1
         END DO
 999     WRITE(LUNOUT,*) ' '
         WRITE(LUNOUT,*) 'total number of beamlets = ', NPOINTS
         WRITE(LUNOUT,*) 'maximum energy (MeV/u) = ', MAX_EN
         WRITE(LUNOUT,*) ' '
         WRITE(LUNOUT,*) 'Example of particle distribution created'
         WRITE(LUNOUT,*) '(first 10 and last 10 positions only)'

* if gaussian beam distribution is selected, use file FWHM
* if rectangular distribtuion is selected, use BEAM width
         IF (LDXGSS) THEN
            WRITE(LUNOUT,*) 'Sampling gaussian dist. in x'
         ELSE
            WRITE(LUNOUT,*) 'Sampling rectangular dist. in x'
         END IF
         
         IF (LDYGSS) THEN
            WRITE(LUNOUT,*) 'Sampling gaussian dist. in y'
         ELSE
            WRITE(LUNOUT,*) 'Sampling rectangular dist. in y'
         END IF

* normalize distribution and prints sample
         DO I = 1, NPOINTS
            USRMIN(I) = USRMIN(I) / USRMAX(NPOINTS)
            USRMAX(I) = USRMAX(I) / USRMAX(NPOINTS)
            IF ((I .LT. 10) .OR. (I .GT. (NPOINTS-10))) THEN
               WRITE(LUNOUT,*) USREN(I), USRFOC(I), USRX(I),
     &         USRY(I), USRMIN(I), USRMAX(I)
            END IF
         END DO
 
* CHECK IF ENERGY IS LARGER THAN SPECIFIED IN BEAM CARD
* ======================================================        
* currently only valid for C-12 and protons
* IJBEAM = -2 is assumed to be C-12 and not other HEAVYION
* IJBEAM = 1: proton

* for C-12
         IF (IJBEAM.EQ.-2) THEN
            MAX_EN = MAX_EN*TWELVE

* for protons
         ELSE IF (IJBEAM.EQ.1) THEN
            MAX_EN = MAX_EN*AMPRMU
* for other particles
         ELSE

         END IF

* check and message
         IF (PBEAM.LT.ZERZER) THEN
            IF (MAX_EN.GT.-PBEAM) THEN
               WRITE(LUNOUT,*) 'Warning: energy of the sample particle'
               WRITE(LUNOUT,*) 'is larger than the energy in beam card!'
               WRITE(LUNOUT,*) 'Increase energy in beam card!'
               WRITE(LUNOUT,*) '(Yukihara)'
            END IF
         END IF

*  |  *** User initialization ***
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
*  Npflka is the stack counter: of course any time source is called it
*  must be =0
      NPFLKA = NPFLKA + 1
*  Wt is the weight of the particle
      WTFLK  (NPFLKA) = ONEONE
      WEIPRI = WEIPRI + WTFLK (NPFLKA)
*  Particle type (1=proton.....). Ijbeam is the type set by the BEAM
*  card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         ILOFLK (NPFLKA) = IJBEAM
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
      END IF
*  |
*  +-------------------------------------------------------------------*
*  From this point .....
*  Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
*  User dependent flag:
      LOUSE  (NPFLKA) = 0
*  No channeling:
      LCHFLK (NPFLKA) = .FALSE.
      DCHFLK (NPFLKA) = ZERZER
*  User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
*  User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
*  Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
*  ... to this point: don't change anything


* SAMPLE FROM THE USER-PROVIDED DISTRIBUTION
* ============================================================
* The following routine picks a particle based on the data in
* the source file provided by the user, specified in the 
* OPEN card with logical unit 98.

* Important: the particle in the BEAM card needs to match the
* information in user-provided file.


      EN = ZERZER
      RND = FLRNDM(XDUMMY)
      DO I = 1, NPOINTS
         IF ((RND .GE. USRMIN(I)) .AND. (RND .LT. USRMAX(I))) THEN
            XPOS = USRX(I)
            YPOS = USRY(I)
            EN = USREN(I)
            FOC = USRFOC(I)
         END IF
      END DO

      IF (IJBEAM.EQ.-2) THEN
         EN = EN*TWELVE

      ELSE IF (IJBEAM.EQ.1) THEN
         EN = EN*AMPRMU

* other particles are not handled by the code at the moment
      ELSE
      
      END IF

* indicate problems when no particle is found and energy is zero
      IF (EN .LE. 0) THEN
         WRITE(LUNOUT,*) 'Problem with particle sampling!!!'
         WRITE(LUNOUT,*) '(Yukihara)'
      END IF

* =====================================================================

*  Particle age (s)
      AGESTK (NPFLKA) = +ZERZER
      AKNSHR (NPFLKA) = -TWOTWO

*  Kinetic energy of the particle 

*      TKEFLK (NPFLKA) = SQRT ( PBEAM**2 + AM (IONID)**2 ) - AM (IONID)
      TKEFLK(NPFLKA) = EN

*  Particle momentum
*      PMOFLK (NPFLKA) = PBEAM
*     PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
*    &                       + TWOTWO * AM (IONID) ) )
      PMOFLK(NPFLKA) = SQRT(EN*(EN + TWOTWO * AM(IONID)))

*  Cosines (tx,ty,tz)
      TXFLK  (NPFLKA) = UBEAM
      TYFLK  (NPFLKA) = VBEAM
      TZFLK  (NPFLKA) = WBEAM
*     TZFLK  (NPFLKA) = SQRT ( ONEONE - TXFLK (NPFLKA)**2
*    &                       - TYFLK (NPFLKA)**2 )

*  Polarization cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER

* If BEAM card request Gaussian Beam, then the width is taken
* from the user-provided file; if BEAM card requests Rectangular
* beam, the distribution is as in BEAM card

* for x position
      CALL FLNRR2(RGAUS1,RGAUS2)
      IF (LDXGSS) THEN
         FWHMX = FOC
         SIGX = FWHMX/S2FWHM
         XFLK(NPFLKA) = XPOS + RGAUS1*SIGX
      ELSE
         XFLK(NPFLKA) = XBEAM + XSPOT*(FLRNDM(XDUMMY) - HLFHLF)
      END IF

* for y position
      IF (LDYGSS) THEN
         FWHMY = FOC
         SIGY = FWHMY/S2FWHM
         YFLK(NPFLKA) = YPOS + RGAUS2*SIGY
      ELSE
         YFLK(NPFLKA) = YBEAM + YSPOT*(FLRNDM(XDUMMY) - HLFHLF)
      END IF

* for z position
      ZFLK(NPFLKA) = ZBEAM





* DON'T CHANGE THIS PART
* =================================================================
*  Calculate the total kinetic energy of the primaries: don't change
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
      RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV

      RETURN
      END

