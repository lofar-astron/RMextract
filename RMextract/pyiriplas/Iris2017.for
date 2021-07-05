C IRIS2017.FOR_____________________________ May.2017 T.L.Gulyaeva
C IRIS2015c.FOR_____________________________ Feb.2016 T.L.Gulyaeva
C IRIS2015.FOR_____________________________ Jun.2015 T.L.Gulyaeva
C IRIS2014.FOR_____________________________ Apr.2014 T.L.Gulyaeva
C IRIS2013.FOR_____________________________ Oct.2013 T.L.Gulyaeva
C IRIS2011.FOR_____________________________ Jan.2012 T.L.Gulyaeva
C IRIS2010.FOR_____________________________ Aug 2010 T.L.Gulyaeva
C IRIS2009.FOR_____________________________ Feb 2010 T.L.Gulyaeva
C IRIS2007.FOR_____________________________ Mar 2008 T.L.Gulyaeva
C IRIS2006.FOR_____________________________ Oct 2006 T.L.Gulyaeva
C IRIS2003.FOR _____________________________Nov 2004 T.L.Gulyaeva
C IRIS16.FOR ---------- Oct 95 (D Bilitza) :Sep 2001 T L Gulyaeva
C
C
C*****************************************************************
C                        ISO_IRI-Plas PROJECT
C
C********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
C***RUSSIAN STANDARD MODEL OF IONOSPHERE AND PLASMASPHERE (SMI)***
C*****************************************************************
C**************** ALL-IN-ONE SUBROUTINE  *************************
C*****************************************************************
C includes IRIS2017.FOR (computes parameters for specified
C altitude range, for specified altitude, latitude, longitude, 
C year, month, day,day of year, or hour range) 
C IRI2001 STORM-MODEL INCLUDED from IRI2001 (D.Bilitza, June 2001)
C
C IRI2001 functions B0_98 (instead of B0_TAB), xe3_1, xe4_1 included
C B1-diurnal change included with B0_GUL option.(T.Gulyaeva, Sep. 2005)
C
C SMI: HPL,XNEPL,NE(H) AND TECPL FOR 1336 KM : 20000 KM INCLUDED
C Changed Feb. 2008 : SMI extension at [h05top : hpl]. hpl<=20000 km !!!
C
C SUBROUTINE GKPM PRODUCING KP AND KPM INDICES INCLUDED 
C (T.Gulyaeva, Sep. 2001)
C*****************************************************************
C
C*05/05/2017 Solar activity smoothed proxy index specified by jind:
C*  jind =1 ssn1_12.dat  DEFAULT: SSN1  =>  rz(3); F10.7  =>  sf(3) 
C*  jind =2 ssn2_12.dat  SSN2 converted to SSN1  =>  rz(3); F10.7  => sf(3)
C*  jind =3 f107_12.dat  F10.7    =>  sf(3) converted to SSN1  =>  rz(3) 
C*  jind =4 geci_12.dat  GEC => rz(3); F10.7  => sf(3)
C*  jind =5 teci_12.dat  TEC => rz(3); F10.7  => sf(3)
C*  jind =6 igin_12.dat  IG  => rz(3); F10.7  => sf(3)
C*  jind =7 mgii_12.dat  MgII => rz(3); F10.7  => sf(3)
C*  Jind =8 lyma_12.dat  Lyman-alpha => rz(3); F10.7  => sf(3)
C* Ref. T.L.Gulyaeva et al., 2017. TEC proxy index of solar activity for the
C*      International Reference Ionosphere IRI and its extension to 
C*      Plasmasphere IRI-PLAS model. Int. J. Sci. Eng. Applied Sci., 3, 5, 144-150,
C*      http://ijseas.com/index.php/issue-archive-2/volume3/issue-5/
C
      SUBROUTINE IRIS2017(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,    
     &    HEIBEG,HEIEND,NUMHEI,OUTF,JIND)                          
C-----------------------------------------------------------------
C
C     INPUT: ALATI,ALONG  LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C         IYYYY         Year as YYYY, e.g. 1985
C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL
C                          HOURS
C         HEIBEG,       HEIGHT RANGE IN KM; maximal 50 heights, i.e.
C         HEIEND,HEISTP  int((heiend-heibeg)/heistp)+1.le.50
C     UP TO 120 HEIGHTS ARE ALLOWED FOR [80 KM : 20,000 KM]
C     if jf(8)=false OARR(1)=input value for foF2 or NmF2
C     if jf(9)=false  OARR(2)=input value for hmF2 or M(3000)F2
C     if jf(10)=false   OARR(3)=input value for Ne(300km)
C     if jf(10)=false   OARR(4)=input value for Ne(400km)
C     if jf(10)=false   OARR(5)=input value for Ne(600km)
C     if jf(17)=false   OARR(33)=input value for sunspot number Rzs
C     if jf(18)=false   NO OUTPUT Ne(h) PROFILE
C     if jf(26)=false   NO STORM foF2 calculation
C
C         JF(1:30)      TRUE/FALSE FLAGS FOR SEVERAL OPTIONS
C          JF(1)=.TRUE.[.FALSE.]   ELECTRON DENSITY IS [NOT] CALCULATED
C          JF(2)=T[F]    TEMPERATURES ARE [NOT] CALCULATED
C          JF(3)=T[F]    ION COMPOSITION IS [NOT] CALCULATED
C          JF(4)=T[F]    B0 FROM TABLE [FROM GULYEAVA 1987]
C          JF(5)=T[F]    F2 PEAK FROM CCIR [FROM URSI]
C WITHDRAWN: JF(6)=T[F]    ION COMP. STANDARD [DANILOV-YAICHNIKOV-1985]
C          JF(7)=T[F]    STAND. IRI TOPSIDE [IRI-79]
C          JF(8)=T[F]    NMF2 PEAK MODEL [INPUT VALUES]
C          JF(9)=T[F]    HMF2 PEAK MODEL [INPUT VALUES]
C          JF(10)=T[F]   TE MODEL [TE-NE MODEL WITH NE INPUT]
C          JF(11)=T[F]   NE STANDARD [LAY-FUNCTIONS VERSION]
C          JF(12)=T[F]   MESSAGE ARE WRITTEN TO UNIT=6 [=12]
C          JF(17)=T[F]   sunspot number from IG_RZ file [Rz input]
C          JF(18)=T[F]   OUTPUT NE(h)&fN(h) PROFILES [NO PROFILES]
C          JF(26)=T[F]   OUTPUT foF2 STORM model [NO STORM MODEL]
C+++++++++++++++++++++^^^^
C          jf(21) =.true.      OARR(41)=user input for daily F10.7 index
C          jf(23) =.false.     OARR(41)=user input for daily F10.7 index
C          optional for jf(21:24); default is F10.7D=COV
C          jf(25) =.false.     OARR(41)=user input for daily F10.7 index
C          if oarr(41).le.0 then 12-month running mean is 
C          taken from internal file]
C*         jf(27) =.false.     OARR(39)=user input for selected solar proxy 
C+++++++++++++++++++++vvvv
C
C     JF(1:11)=.TRUE. GENERATES THE STANDARD IRI-90 PARAMETERS.
C     IF YOU SET JF(8)=.FALSE., THAN YOU HAVE TO PROVIDE THE F2 PEAK
C     NMF2/M-3 OR FOF2/MHZ IN OARR(1). SIMILARLY, IF YOU SET JF(9)=
C     .FALSE., THAN YOU HAVE TO PROVIDE THE F2 PEAK HEIGHT HMF2/KM IN
C     OARR(2). IF YOU SET JF(10)=.FALSE., THAN YOU HAVE TO PROVIDE THE
C     ELECTRON DENSITY IN M-3 AT 300KM AND/OR 400KM AND/OR 600KM IN
C     OARR(3), OARR(4), AND OARR(5). IF YOU WANT TO USE THIS OPTION AT
C     ONLY ONE OF THE THREE ALTITUDES, THAN SET THE DENSITIES AT THE
C     OTHER TWO TO ZERO.
C     JF(17)=.FALSE. REQUIRES INPUT OF SUNSPOT NUMBER RZS
C     JF(18)=.TRUE. GENERATES OUTPUT OF NE(h) AND fN(h) PROFILES
C     JF(26)=.TRUE. GENERATES STORM foF2 CALCULATION
C     JF(30)=.TRUE.[.FALSE.] TEC MODEL [INPUT VALUES]
C
C     OUTPUT:  OUTF(0:11,0:120)
C        OUTF(0,*)  HEIGHT, KM
C        OUTF(1,*)  ELECTRON DENSITY/M-3
C        OUTF(2,*)  NEUTRAL TEMPERATURE/K
C        OUTF(3,*)  ION TEMPERATURE/K
C        OUTF(4,*)  ELECTRON TEMPERATURE/K
C        OUTF(5,*)  O+ ION DENSITY/M-3
C        OUTF(6,*)  H+ ION DENSITY/M-3
C        OUTF(7,*)  HE+ ION DENSITY/M-3
C        OUTF(8,*)  O2+ ION DENSITY/M-3
C        OUTF(9,*)  NO+ ION DENSITY/M-3
C        OUTF(10,*)  PLASMA FREQUENCY/MHz 
C        OUTF(11,*)  LOG NE(M-3) // TECeff in step of 200 km
C !!! MODIFIED FOR OUTPUT RESULTS  OUTF(10,*) & OUTF(11,*) !!!
C =============================================================
C
C            OARR(1:50)   ADDITIONAL OUTPUT PARAMETERS
C
C     OARR(1) = NMF2/M-3            OARR(2) = HMF2/KM
C     OARR(3) = NMF1/M-3            OARR(4) = HMF1/KM
C     OARR(5) = NME/M-3             OARR(6) = HME/KM
C     OARR(7) = NMD/M-3             OARR(8) = HMD/KM
C     OARR(9) = HHALF/KM            OARR(10) = B0/KM
C     OARR(11) =VALLEY-BASE/M-3     OARR(12) = VALLEY-TOP/KM
C     OARR(13) = TE-PEAK/K          OARR(14) = TE-PEAK HEIGHT/KM
C     OARR(15) = TE-MOD(300KM)      OARR(16) = TE-MOD(400KM)/K
C     OARR(17) = TE-MOD(600KM)      OARR(18) = TE-MOD(1400KM)/K
C     OARR(19) = TE-MOD(3000KM)     OARR(20) = TE(120KM)=TN=TI/K
C     OARR(21) = TI-MOD(430KM)      OARR(22) = X/KM, WHERE TE=TI
C     OARR(23) = SOL ZENITH ANG/DEG OARR(24) = SUN DECLINATION/DEG
C     OARR(25) = DIP/deg            OARR(26) = DIP LATITUDE/deg
C     OARR(27) = MODIFIED DIP LAT.  OARR(28) = DELA
C     OARR(29) = sunrise/dec. hours OARR(30) = sunset/dec. hours
C     OARR(31) = ISEASON (1=spring) OARR(32) = RL (L value)
C     OARR(33) = RZS (SUNSPOT INPUT)
C     OARR(35) = B1             OARR(36) = M(3000)F2
C     OARR(37) = LATI  OARR(38) = LONGI
C     OARR(39) = gind  OARR(40) = MLONG
C     OARR(41)=f107d   OARR(42) = MLAT  
C     OARR(43)=daynr   OARR(44) = gmlt  
C     OARR(45)=stormcor OARR(46)= Hscale_top
C-------------------------------------------------------------------
C******************************************************************
C***  THE ALTITUDE LIMITS ARE: LOWER (DAY/NIGHT)  UPPER         ***
C***     ELECTRON DENSITY       60/80 KM     20,000 KM ... hpl  ***
C***     TEMPERATURES           120 KM       10 000 KM          ***
C***     ION DENSITIES          100 KM        1000 KM           ***
C******************************************************************
C*****************************************************************
C*********            INTERNALLY                    **************
C*********       ALL ANGLES ARE IN DEGREE           **************
C*********       ALL DENSITIES ARE IN M-3           **************
C*********       ALL ALTITUDES ARE IN KM            **************
C*********     ELECTRON CONTENT*1E-16 IN M-2        ************** 
C*********     ALL TEMPERATURES ARE IN KELVIN       **************
C*********     ALL TIMES ARE IN DECIMAL HOURS       **************
C*****************************************************************
C*****************************************************************
C*****************************************************************
      INTEGER           DAYNR,DDO,DO2,SEASON,SEADAY,NUMHEI
      REAL              LATI,MO2,MO,MODIP,NMF2,MAGBR
      REAL              NMF1,NME,NMD,MM,MLAT,MLONG,NOBO2,NEI,RZS,COV
      real:: tetitn(3)                         
      CHARACTER FILNAM*12
      CHARACTER*4 xxpr,yypr                 
      DIMENSION  F(3),RIF(4),E(4),XDELS(4),DNDS(4)
      DIMENSION  FF0(988),XM0(441),F2(13,76,2),FM3(9,49,2)
      DIMENSION  FF0N(988),XM0N(441),F2N(13,76,2),FM3N(9,49,2)
      DIMENSION  AMP(4),HXL(4),SCL(4)                  
      DIMENSION  XSM(4),MM(5),DTI(4)
      DIMENSION  AHH(7),STTE(6),DTE(5),ATE(7),TEA(6),XNAR(3)
      DIMENSION  PG1O(80),PG2O(32),PG3O(80),PF1O(12),PF2O(4),PF3O(12)
      DIMENSION  HO(4),MO(5),DDO(4),HO2(2),MO2(3),DO2(2),DION(7)
      DIMENSION  OUTF(0:11,0:120),OARR(50),ARIG(3),RZAR(3)   
      DIMENSION  INDAP(13),AKP(13),icurkp(0:7),iprekp(0:7)
      LOGICAL           EXT,SCHALT,NIGHT,TECON(3),sam_mon,sam_yea,sam_ut
      LOGICAL           F1REG,FOF2IN,HMF2IN,URSIF2,LAYVER,DY,GULB0
      LOGICAL           NODEN,NOTEM,NOION,TENEOP,sam_doy,TECIN
      LOGICAL           OLD79,JF(30),URSIFO,RZIN,INKP,F1_OCPRO,F1_L_COND
      LOGICAL           FOF1IN,HMF1IN,FOEIN,HMEIN,DREG,sam_jind
     &,IGIN,igino,rzino             
     &,sam_date
      COMMON    /BLOCK1/HMF2,NMF2,HMF1,F1REG           /CONST/UMR,PI
     &          /BLOCK2/B0,B1,C1      /BLOCK3/HZ,T,HST,STR
     &          /BLOCK4/HME,NME,HEF   /BLOCK5/NIGHT,E
     &          /BLOCK6/HMD,NMD,HDX   /BLOCK7/D1,XKK,FP30,FP3U,FP1,FP2
     &          /BLOCK8/HS,TNHS,XSM,MM,DTI,MXSM
     &          /BLOTN/XSM1,TEXOS,TLBDH,SIGMA /BLOTE/AHH,ATE1,STTE,DTE
     &          /BLO10/BETA,ETA,DELTA,ZETA       /ARGEXP/ARGMAX
     &          /BLO11/XHI/B6/CGLAT,H0NE,GLT,FG1,FG2
     &          /CONST1/HUMR,DUMR,XFAP,INDAP,AKP,COV,RZS,HMLON
     &          /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU,HTOP,QF0
     &/PLAS/HPL,XNEPL,XKP,ICALLS,RUT,RLT,OARR,TECI,TCB,TCT,TCPL,TEC
     &/dem/W,ylongi,hour,daynr,cn1000,alon,blon,enre,HSC
      common /temod/month
      COMMON /FIKP/icurkp,iprekp,iyear_cur,imn_cur,idy_cur
      COMMON /BXY/xxpr,yypr                     
      EXTERNAL          XE1,XE2,XE3,XE3_1,XE4,XE5,XE6,XXE6,TEDER,ABLON
      DATA     XNAR       /3*0.0/,
     &      XDELS   /3*5.,10./,      DNDS   /.016,.01,2*.016/,
     &      DDO   /9,5,5,25/,        DO2        /5,5/
      COMMON /path/datapath
      character datapath*200

c  save  icalls,rssn,rsat,cov
      DO 7397 KI=0,11
      do 7397 kk=0,120
7397  OUTF(KI,kk)=-1.
C
      do 8398 kind=6,32,1
8398  oarr(kind)=-1.

C
C PROGRAM CONSTANTS
C
        icalls=icalls+1
        ARGMAX=88.0
        PI=ATAN(1.0)*4.
        UMR=PI/180.
        HUMR=PI/12.
        DUMR=PI/182.5
        ALOG2=ALOG(2.)
        ALG100=ALOG(100.)
C# NEW CONSTANTS:
      Y05=.6931473
      QF=1.
      h05top=0.
      stormcor=1.0
      DH2=9999.
c      numhei=int(abs(heiend-heibeg)/abs(heistp))+1
c      if(numhei.gt.50) numhei=50
C
C Code inserted to aleviate block data problem for PC version.
C Thus avoiding DATA statement with parameters from COMMON block.
C
        AHH(1)=120.
        AHH(2)=0.
        AHH(3)=300.
        AHH(4)=400.
        AHH(5)=600.
        AHH(6)=1400.
        AHH(7)=3000.
        DTE(1)=5.
        DTE(2)=5.
        DTE(3)=10.
        DTE(4)=20.
        DTE(5)=20.
        DTI(1)=10.
        DTI(2)=10.
        DTI(3)=20.
        DTI(4)=20.
C
C
C FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
C AGNR=OUTPUT (OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...
C IUCCIR=UNIT NUMBER FOR CCIR COEFFICIENTS ........................
C
      MONITO=6
      IUCCIR=12
      KONSOL=6
c     if(.not.jf(12)) konsol=12
c
c selection of density and ion composition options ..................
c
      NODEN=(.not.jf(1))
      NOTEM=(.not.jf(2))
      NOION=(.not.jf(3))
      DY=(.not.jf(6))
      LAYVER=(.not.jf(11))
      OLD79=(.not.jf(7))
      GULB0=(.not.jf(4))
      TECIN=(.NOT.JF(30))
      F1_OCPRO=(jf(19))
      f1_l_cond=.false.
      if(F1_OCPRO) F1_L_COND=(.not.jf(20))
      DREG=jf(24)
C
      IF (XKP.GE.0.0) THEN 
        INKP=.TRUE.
      ELSE
        INKP=.FALSE.
      ENDIF
c
c rz12 input option ....................................................
c
      RZIN=(.not.jf(17))
       IF(RZIN) THEN
          OARR(33)=RZS
          ARZIN=OARR(33)
      rssn=arzin
        else
          oarr(33)=-1.
          ENDIF
      IGIN=(.not.jf(27))       ! Add
      IF(IGIN) THEN         !
          AIGIN=OARR(39)        !      
               else
          oarr(39)=-1.
               ENDIF              !
       IF(.not.jf(25)) THEN
         f107d=OARR(41)
        else
          oarr(41)=-1.
          ENDIF
c
c F2 peak density ....................................................
c
      AFOF2=0.
      AHMF2=0.
      FOF2IN=(.not.jf(8))
       IF(FOF2IN) THEN
          AFOF2=OARR(1)
          IF(AFOF2.GT.100.) AFOF2=SQRT(AFOF2/1.24E10)
      else
       oarr(1)=-1.
          ENDIF
      URSIF2=(.not.jf(5))
c
c F2 peak altitude ..................................................
c
      HMF2IN=(.not.jf(9))
       IF(HMF2IN) then
      AHMF2=OARR(2)
      else
      oarr(2)=-1.
      endif
c
c F1 peak density ....................................................
c
      FOF1IN=(.not.jf(13))
       IF(FOF1IN) THEN
          AFOF1=OARR(3)
          IF(AFOF1.GT.100.) AFOF1=SQRT(AFOF1/1.24E10)
        else
          oarr(3)=-1.
          ENDIF
c
c F1 peak altitude ..................................................
c
      HMF1IN=(.not.jf(14))
       IF(HMF1IN) then
                AHMF1=OARR(4)
                if(.not.layver.and.(konsol.gt.1)) write(konsol,1939)
1939  format(' *Ne* User input of hmF1 is only possible for the LAY-',
     &          'version')
        else
                oarr(4)=-1.
        endif
c
c E peak density ....................................................
c
      FOEIN=(.not.jf(15))
       IF(FOEIN) THEN
          AFOE=OARR(5)
          IF(AFOE.GT.100.) AFOE=SQRT(AFOE/1.24E10)
        else
          oarr(5)=-1.
          ENDIF
c
c E peak altitude ..................................................
c
      HMEIN=(.not.jf(16))
       IF(HMEIN) then
                AHME=OARR(6)
        else
                oarr(6)=-1.
        endif
c
C TE-NE MODEL OPTION ..............................................
C
      TENEOP=(.not.jf(10))
        IF(TENEOP) THEN
           DO 8154 JXNAR=1,3
              XNAR(JXNAR)=OARR(JXNAR+2)
              TECON(JXNAR)=.FALSE.
8154          IF(XNAR(JXNAR).GT.0.) TECON(JXNAR)=.TRUE.
      else
               oarr(15)=-1.  !
              oarr(16)=-1.  !  
C!      oarr(3)=-1.
C!      oarr(4)=-1.
C!      oarr(5)=-1.
           ENDIF
c
c lists the selected options before starting the table
c
      if(icalls.gt.1) goto 8201
      write(*,*) '*** ISO_IRI parameters are being calculated ***'
c      write(*,*) '*** ISO_IRI parameters are being calculated ***'
      if(NODEN) goto 2889
c         if(LAYVER) write(*,*) 'Ne, E-F: The LAY-Version is ',
c     &      'prelimenary. Erroneous profile features can occur.'
c         if(GULB0) write(*,*) 'Ne, B0: Bottomside thickness is ',
c     &      'obtained with Gulyaeva-1987 model.'
      if(tecin) then 
      write(*,*) 'TEC: input values are used.'
      write(*,*) '    hmF2: CCIR model is fitted to TEC.' 
      write(*,*) 'Ne, foF2: CCIR model is fitted to TEC.'
      goto 8201
      endif
         if(OLD79) write(*,*) 'Ne: Using IRI-79. Correction',
     &      ' of equatorial topside is not included.'
      if (.not.(HMF2IN)) THEN
      write(*,*) '    hmF2: CCIR model is used.' 
       goto 700
                      ELSE
      if(hmf2.GT.0.) then
        write(*,*) 'hmF2: Input values are used.'
      else
        write(*,*) 'M3000F2: Input values are used.'
      endif
                      ENDIF
 700  continue
         if(FOF2IN) then
          write(*,*) 'Ne, foF2: Input values are used.'
            goto 2889
            endif
         if(URSIF2) then
          write(*,*) 'Ne, foF2: URSI model is used.'
         else
          write(*,*) 'Ne, foF2: CCIR model is used.'
         endif
      if(jf(26)) then
       write(*,*) 'Ne, foF2:  STORM  model included.'
       write(*,*) 'hmF2(foF2) STORM2 model included.'
      endif
2889     continue
      if((.not.NOION).and.(DY))
     &    write(*,*) 'Ion Com.: Using Danilov-Yaichnikov-1985.'
         if((.not.NOTEM).and.(TENEOP))
     &    write(*,*) 'Te: Temperature-density correlation is used.'
8201    continue
C
C****  
C CALCULATION OF DAY OF YEAR AND SUN DECLINATION......................
C CALCULATION OF UT/LT AND RELATED YEAR, MONTH, DAYNRs ...............
C CALCULATION OF (UT-)SEASON (SUMMER=2, WINTER=4).....................
C
        iyear=iyyyy
       if(iyear.lt.100) iyear=iyyyy+1900
      idayy=365
        yds=365.
        if(iyear/4*4.eq.iyear) then
        idayy=366
        yds=366.
        endif
      ryear=iyear*1.0
        if(MMDD.lt.0) then
                DAYNR=-MMDD
                call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
        else
                MONTH=MMDD/100
                IDAY=MMDD-MONTH*100
                call MODA(0,iyear,MONTH,IDAY,DAYNR,nrdaym)
        endif
C ADD+++
        ryear = iyear + (daynr-1.0)/idayy
c
c lyear,lmonth,lday,ldaynre,lnrday related to LT
c
c
        lyear=iyear
        lmonth=month
        lday=iday
        ldaynr=daynr
        lnrday=nrdaym
C
C
C CALCULATION OF GEOG. OR GEOM. COORDINATES IN DEG....................
C CALCULATION OF MAGNETIC INCLINATION (DIP), DECLINATION (DEC)........
C   DIP LATITUDE (MAGBR) AND MODIFIED DIP (MODIP). ALL IN DEGREE......
C
       if(along.lt.0.0) along = along + 360.
        IF(JMAG.GT.0) THEN
           MLAT=ALATI
           MLONG=ALONG
        ELSE
           LATI=ALATI
           YLONGI=ALONG
        ENDIF
cc
      CALL GEODIP(IYEAR,LATI,YLONGI,MLAT,MLONG,JMAG)
C
C_ADD
      call igrf_dip(lati,ylongi,ryear,300.0,dec,dip,magbr,modip) ! NEW IGRF2015
        IF (MAGBR.GT.9998.) THEN          
        CALL FIELDG(LATI,LONGI,300.0,XMA,YMA,ZMA,BET,DIP,DEC,MODIP) !OLD POGO-75
        MAGBR=ATAN(0.5*TAN(DIP*UMR))/UMR
        ENDIF
C_END_ADD
        ABSLAT=ABS(LATI)       !
        ABSMLT=ABS(MLAT)
        ABSMDP=ABS(MODIP)
        ABSMBR=ABS(MAGBR)
C+  NEW: SUBR. CONVER TO PRODUCE CORRECTED GEOMAG. LATITUDE
       CALL CONVER (LATI,YLONGI,CGLAT)
C
C !!!!NEW : SMI SUBR. ABLON TO PREPARE FOR PLASMASPHERE MODEL!!!!
      alon=365.
      blon=365.
      if (((YLONGI.ge.130.).and.(YLONGI.le.178.)).or.((YLONGI.ge.330.).
     &and.(YLONGI.le.340.))) goto 209
      CALL ABLON(LATI,YLONGI,ALON,BLON)
  209 continue
C L-value in Earth's radii RL, height over magnetic equator r0km:
      RL=1./(cos(mlat*umr))**2.
      r0km=RL*6370.
      continue
C ========================================
C
      IF(DHOUR.LT.24.1) goto 2619
        UT=DHOUR-25.
        iytmp=iyear
        idtmp=daynr
        call ut_lt(0,ut,hour,YLONGI,iytmp,idtmp)
       if(idtmp.ne.ldaynr) then
                        lyear=iytmp
                        ldaynr=idtmp
                        call MODA(1,lyear,LMONTH,LDAY,LDAYNR,lnrday)
                     endif
        goto 2629
2619  HOUR=DHOUR
        iytmp=lyear
        idtmp=ldaynr
                call ut_lt(1,ut,hour,YLONGI,iytmp,idtmp)
                if(idtmp.ne.daynr) then
                        iyear=iytmp
                        daynr=idtmp
      idayinp=iday           ! input day 
                        call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
C 
      if (iday.ne.idayinp) mmdd=month*100+iday  ! corrected day according to UT
                        endif
c
c zmonth is decimal month (Jan 1 = 1.0 and Dec 31 = 12.97) 
c
2629  zmonth = lmonth + (lday*1.)/lnrday
      RUT=UT
      RLT=HOUR
C+ Geomagnetic local time
      UTS=ut*3600.        ! UT /sec.
      xday=float(daynr)
      gmlt=fmlt(uts,lati,YLONGI,xday)
C Geomagnetic Local Time, GLT:
C>>      GMLT=UT+(MLONG-69.8)/15.
C>>      IF ((GMLT.LT.0.).OR.(GMLT.GE.24.)) GMLT=ABS(ABS(GMLT)-24.)
C>>      IF (GMLT.EQ.24.) GMLT=0.
C
      SEASON=INT((DAYNR+45.0)/92.0)
      IF(SEASON.LT.1) SEASON=4
      NSEASN=SEASON
      seaday=daynr
      seamon=zmonth      !
C+ TLG Day-of-year reduced to 360.=> sday=season_day
      sday=daynr/yds*360.         ! north hemisphere
      IF(LATI.GT.0.0) GOTO 5592
        SEASON=SEASON-2
        IF(SEASON.LT.1) SEASON=SEASON+4
      seamon=zmonth+6.              
         if(seamon.ge.13) seamon=seamon-12.
         seaday=daynr+idayy/2.                     !
         if(seaday.gt.idayy) seaday=seaday-idayy     !
C TLG : sday for the sothern hemisphere:
      sday=sday+180.                ! south hemisphere
      if (sday.gt.360.) sday=sday-360. ! south hemisphere
C
C MEAN F10.7CM SOLAR RADIO FLUX (COV), GLOBAL IONOSPHERIC INDEX
C (GIND)
C
C 12-month running mean sunspot number (rssn) and Ionospheric Global 
C index (gind), daily F10.7 cm solar radio flux (f107d) and monthly 
C F10.7 (cov) index   
C
5592    sam_mon=(month.eq.montho)
      sam_yea=(iyear.eq.iyearo)
      sam_doy=(daynr.eq.idaynro)
      sam_jind = (jind.eq.jindo)
      sam_date=(sam_yea.and.sam_doy.and.sam_jind)   !
        sam_ut=(ut.eq.ut0)         !
        if(sam_date.and..not.rzino.and..not.rzin.       !
     &                   and..not.igin.and..not.igino) goto 2910  !
C     if you want to use the ig_rz.dat version of the coefficients, please use 
C the statements commented                                   <! ig_rz.dat> 
C      call tcon(iyear,month,iday,daynr,rzar,arig,ttt,nmonth) !  ig_rz.dat
C
C NEW TLG: the following statements commented <! gec_rz.dat> are using gec_rz.dat file:
      IF (jind.eq.0) THEN                    
      call tcongec(iyear,month,iday,daynr,rzar,arig,ttt,nmonth)   ! gec_rz.dat
                     ELSE                     
C* NEW TLG May 2017: the following statement for different solar proxies input 
C* R12 proxy - rsar(3), F10.7 proxy - arig(3); jind - option for proxy input: 
      call tconind(iyear,month,iday,daynr,rzar,arig,ttt,nmonth,jind)
                     ENDIF
C
      if(nmonth.lt.0) goto 3330
                if(RZIN) then                 !+
                  rrr = arzin
                  rzar(1) = rrr
                  rzar(2) = rrr
                  rzar(3) = rrr
C                 zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr  ! ig_rz.dat
         IF (jind.eq.0) THEN                                     
                  zi=((0.0194*rrr+1.0906)-1.0)*50.0              
                  if(zi.gt.174.0) zi=174.0
                  if(zi.lt.0.) zi=10.      ! check low limit of IG-index
                  arig(1) = zi
                  arig(2) = zi
                  arig(3) = zi                      
         ELSE
                  cov=62.6645+rzar(3)*0.9066             !* F10.7 proxy
                  arig(1) = cov                     
                  arig(2) = cov                       !*
                  arig(3) = cov                       !*
               gind= 0.97236*rzar(3)+0.18836          !*
         ENDIF                                           !*
      endif                   ! 
C
C!++++++
C        if(IGIN) then           ! ig_rz.dat
C           zi = aigin              ! ig_rz.dat
C           arig(1) = zi            ! ig_rz.dat
C           arig(2) = zi            ! ig_rz.dat
C           arig(3) = zi            ! ig_rz.dat
C           endif                ! ig_rz.dat
        rssn=rzar(3)
      IF (jind.eq.0) THEN   !*
      gind=arig(3)      
      cov=63.75+RSSN*(0.728+RSSN*0.00089)
      ELSE                 !*
        cov=arig(3)
      ENDIF                !*
      COVSAT=cov
        if(covsat.gt.188.) covsat=188.
         f107d=cov
         if(.not.jf(25)) then       !* 
           f107d=oarr(41)           !* 
C*           RSSN=1.076*f107d-65.7817  !*
         endif                      !* 
C
C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG).........................
C NOON VALUE (XHINON).................................................
C
2910  continue
      if(.not.rzin) RZS=rssn     ! 
         rg=rzs               ! 
        CALL SOCO(ldaynr,HOUR,LATI,YLONGI,SUNDEC,XHI,SAX,SUX)
        CALL SOCO(ldaynr,12.0,LATI,YLONGI,SUNDE1,XHINON,SAXNON,SUXNON)
C New: call APF for COV81, RZ81
            call subapf(iyear,month,iday,ut,indap,covin,rzc) !
      nmono=month
C New: call GKPM for calculation of kp_indices and Kpm_index
C
      IF (INKP) THEN 
        GOTO 4
      ELSE 
        CALL GKPM(indap,akp,xkp)
      ENDIF 
   4  CONTINUE
C
C
              NIGHT=.FALSE.
        if(abs(sax).gt.25.0) then
         if(sax.lt.0.0) NIGHT=.TRUE.
         goto 1334
         endif
        if(SAX.le.SUX) goto 1386
         if((hour.gt.sux).and.(hour.lt.sax)) night=.true.
        goto 1334
1386     IF((HOUR.GT.SUX).OR.(HOUR.LT.SAX)) NIGHT=.TRUE.
C
C CALCULATION OF ELECTRON DENSITY PARAMETERS................
C lower height boundary (HNEA), upper boundary (HNEE)
C
1334  HNEA=65.
      IF(NIGHT) HNEA=80.
      HNEE=2000.
      IF(NODEN) GOTO 4933
      DELA=4.32
      IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)
      DELL=1+EXP(-(ABSLAT-20.)/10.)
C     
C E peak critical frequency (foE), density (NmE), and peak height (hmE)
c
        IF(FOEIN) THEN
          FOE=AFOE
          IF(AFOE.GT.100.) AFOE=SQRT(AFOE/1.24E10)
        ELSE
            FOE=FOEEDI(COV,XHI,XHINON,ABSLAT)
      ENDIF
       NME=1.24E10*FOE*FOE
c
c E peak altitude ..................................................
        IF(HMEIN) THEN
          HME=AHME
        ELSE
          HME=110.0
        ENDIF
C
c F2 peak critical frequency foF2, density NmF2, and height hmF2
c
C READ CCIR AND URSI COEFFICIENT SET FOR CHOSEN MONTH ............
C
 499  CONTINUE
C
      IF(URSIF2.NEQV.URSIFO) GOTO 7797
       if(.not.rzin.and..not.rzino.and..not.igin.and..not.igino) then
      IF(sam_mon.AND.(nmonth.EQ.nmono).and.sam_yea) GOTO 4292
      IF(sam_mon) GOTO 4293
       endif
c
c the program expects the coefficients files in ASCII format; if you
C want to use the binary version of the coefficients, please use the
C the statements that are commented-out below and comment-out the
C ASCII-related statements.
c
7797    URSIFO=URSIF2
      WRITE(FILNAM,104) MONTH+10
 104   FORMAT('ccir',I2,'.asc')
c-binary- if binary files than use:
c-binary-104   FORMAT('ccir',I2,'.bin')
c-web- special for web-version:
c-web-104   FORMAT('/usr/local/etc/httpd/cgi-bin/models/IRI/ccir',I2,'.asc')
c
c344    lread=1
        OPEN(IUCCIR,FILE=TRIM(ADJUSTL(datapath))//FILNAM,STATUS='OLD',
     &          ERR=8448,FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &          FORM='UNFORMATTED')
        READ(IUCCIR,4689) F2,FM3
 4689    FORMAT(1X,4E15.8)
c-binary- if binary files than use:
c-binary-        READ(IUCCIR) F2,FM3
        CLOSE(IUCCIR)
C
C then URSI if chosen ....................................
C
        if(URSIF2) then
          WRITE(FILNAM,1144) MONTH+10
 1144  FORMAT('ursi',I2,'.ASC')
c-web- special for web-version:
c-web-1144       FORMAT('/usr/local/etc/httpd/cgi-bin/models/IRI/ursi',I2,'.asc')
c-binary- if binary files than use:
c-binary-1144       FORMAT('ursi',I2,'.bin')
          OPEN(IUCCIR,FILE=TRIM(ADJUSTL(datapath))//FILNAM,STATUS='OLD',
     &         ERR=8448,FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &         FORM='UNFORMATTED')
          READ(IUCCIR,4689) F2
c-binary- if binary files than use:
c-binary-          READ(IUCCIR) F2
          CLOSE(IUCCIR)
        endif
C
C READ CCIR AND URSI COEFFICIENT SET FOR NMONTH, i.e. previous
c month if day is less than 15 and following month otherwise
C
c      MONTHO=MONTH
c      GOTO 4291
C
4293    continue
c
c first CCIR ..............................................
c
        WRITE(FILNAM,104) NMONTH+10
        OPEN(IUCCIR,FILE=TRIM(ADJUSTL(datapath))//FILNAM,STATUS='OLD',
     &          ERR=8448,FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &          FORM='unFORMATTED')
        READ(IUCCIR,4689) F2N,FM3N
c-binary- if binary files than use:
c-binary-        READ(IUCCIR) F2N,FM3N
        CLOSE(IUCCIR)
C
C then URSI if chosen .....................................
C
        if(URSIF2) then
          WRITE(FILNAM,1144) NMONTH+10
          OPEN(IUCCIR,FILE=TRIM(ADJUSTL(datapath))//FILNAM,STATUS='OLD',
     &         ERR=8448,FORM='FORMATTED')
c-binary- if binary files than use:
c-binary-     &         FORM='unFORMATTED')
          READ(IUCCIR,4689) F2N
c-binary- if binary files than use:
c-binary-          READ(IUCCIR) F2N
          CLOSE(IUCCIR)
          endif
        GOTO 4291
8448    WRITE(MONITO,8449) FILNAM
8449    FORMAT(1X////,
     &    ' The file ',A30,'is not in your directory.')
        GOTO 3330
C
C* LINEAR INTERPOLATION IN SOLAR ACTIVITY. solar proxy used for foF2
C
4291   CONTINUE
      IF (JIND.eq.0) THEN      !*TEST
      RR2=ARIG(1)/100.
        RR2N=ARIG(2)/100.
                      ELSE      !*TEST                           
           RR2=RZAR(1)/100.     !*TEST
            RR2N=RZAR(2)/100.   !*TEST
                      ENDIF     !*TEST
        RR1=1.-RR2
        RR1N=1.-RR2N
        DO 20 I=1,76
        DO 20 J=1,13
        K=J+13*(I-1)
        FF0N(K)=F2N(J,I,1)*RR1N+F2N(J,I,2)*RR2N
20      FF0(K)=F2(J,I,1)*RR1+F2(J,I,2)*RR2
        RR2=RZAR(1)/100.
        RR2N=RZAR(2)/100.
        RR1=1.-RR2
        RR1N=1.-RR2N
        DO 30 I=1,49
        DO 30 J=1,9
        K=J+9*(I-1)
        XM0N(K)=FM3N(J,I,1)*RR1N+FM3N(J,I,2)*RR2N
30      XM0(K)=FM3(J,I,1)*RR1+FM3(J,I,2)*RR2
4292  zfof2  =  FOUT(MODIP,LATI,YLONGI,UT,FF0)
      fof2n  =  FOUT(MODIP,LATI,YLONGI,UT,FF0N)
      zm3000 = XMOUT(MODIP,LATI,YLONGI,UT,XM0)
      xm300n = XMOUT(MODIP,LATI,YLONGI,UT,XM0N)
      midm=15
      if(month.eq.2) midm=14
      if (iday.lt.midm) then
                yfof2 = fof2n + ttt * (zfof2-fof2n)
      xm3000= xm300n+ ttt * (zm3000-xm300n)
        else
                yfof2 = zfof2 + ttt * (fof2n-zfof2)
      xm3000 = zm3000+ ttt * (xm300n-zm3000)
      endif                                       
501   ratfc=1.  
      IF(FOF2IN) THEN
           FOF2=AFOF2
        ELSE
          FOF2=YFOF2
        ENDIF
C new : flag for one-step-iteration with TECI input:
      iflagtc=0
C ************** NEW FROM IRI-2001********************
C STORM model updating
C
1502  stormcor=-1.
      icoord=1
      if (xfap.lt.200.) XFAP=0.0
      KUT=INT(UT)
      if (jf(26)) then
      xfap=0.
      rlat=0.
      call STORM(indap,lati,YLONGI,icoord,rlat,kut,daynr,stormcor,xfap)
CREM      FOF2=FOF2*stormcor
      FOF2storm=FOF2*stormcor  !NEW++++++++++ for STORM2
      ratfc=stormcor
      endif
C ****************************************************
C
C *** SMI PLASMAPASUE ALTITUDE 
         W=RSSN
      EH=20200.
      hprr=cos(mlat*umr)
      HHH=((5.7*hprr*hprr-1.)*6370.+500.)*(1.-0.0825*XKP)    ! for hpl(kp)
      HPL=20200.  !  for GPS-TEC 
      heiend=hpl
C
C F2 REGION PEAK HEIGHT
C
      IF (HMF2IN) then
        IF(AHMF2.LT.50.0) THEN 
         XM3000=AHMF2
COLD    HMF2=HMF2EDS(RZS,FOF2/FOE,XM3000,HOUR,DAYNR,ABSMLT) ! SMI hmF2
      HMF2=PEAKH(FOE,FOF2,XM3000)
        ENDIF
      ELSE  
cIRI    HMF2=HMF2ED(MAGBR,RSSN,FOF2/FOE,XM3000)
      IF (DH2.GT.9998.) then
COLD      HMF2=HMF2EDS(RZS,FOF2/FOE,XM3000,HOUR,DAYNR,ABSMLT) ! SMI
      HMF2=PEAKH(FOE,FOF2,XM3000)  ! CCIR prediction
C Add  Include hmF2 STORM2 correction:
      IF (ratfc.ne.1.0) THEN
C Use model dh(dne/dnm):
C
      dlogne=2.0*alog10(ratfc)  !new
      call subdlogh(dlogne,mlat,sday,rzs,dlgh) !new
      rathhm=10.**dlgh !new
      if (rathhm.lt.0.8) rathhm=0.8
      HMF2=HMF2*rathhm !new
      FOF2=FOF2storm
C+++==========================
          call ROGUL(SEADAY,XHI,SEAX,GRAT)
      if(NIGHT) GRAT=0.91-HMF2/4000.
          B0CNEW=HMF2*(1.-GRAT)
         BCOEF=B1*(B1*(0.0046*B1-0.0548)+0.2546)+0.3606
      B0=B0CNEW/BCOEF        ! NEW coefficient BCOEF
      ENDIF
      endif
      ENDIF
C ****************************************************
        NMF2=1.24E10*FOF2*FOF2
C Include from IRI2000:
        nmono=nmonth
        MONTHO=MONTH
        iyearo=iyear
        idaynro=daynr
        rzino=rzin
        igino=igin
        ut0=ut
        jindo=jind
c
c topside profile parameters .............................
c
  836    COS2=COS(MLAT*UMR)
      COS2=COS2*COS2
      FLU=(COVSAT-40.0)/30.0
c option to use unlimiited F10.7M for the topside
       IF(OLD79) FLU=(COV-40.0)/30.0
        EX=EXP(-MLAT/15.)
        EX1=EX+1
        EPIN=4.*EX/(EX1*EX1)
        ETA1=-0.02*EPIN
      ETA=0.058798+ETA1+FLU*(-0.014065+0.0069724*COS2)+
     &(0.0024287+0.0042810*COS2-0.0001528*FOF2)*FOF2
      ZETA=0.078922-0.0046702*COS2-0.019132*FLU+0.0076545*FLU*COS2+
     &(0.0032513+0.0060290*COS2-0.00020872*FOF2)*FOF2
      BETA=-128.03+20.253*COS2+FLU*(-8.0755-0.65896*COS2)+(0.44041
     &+0.71458*COS2-0.042966*FOF2)*FOF2
      Z=EXP(94.5/BETA)
      Z1=Z+1
      Z2=Z/(BETA*Z1*Z1)
      DELTA=(ETA/Z1-ZETA/2.0)/(ETA*Z2+ZETA/400.0)
C# NEW!!! T.L.Gulyaeva CORRECTION OF BENT MODEL FOR h05top HEIGHT
      hei05=0.
         CALL TOPH05(RSSN,MLAT,HOUR,HMF2,HEI05,sday)
c
      h05top=hei05
C# NEXT CALCULATION OF XNE05 IS REQUIRED TO ESTIMATE QF-FACTOR
      HTOP=H05TOP             ! transition from IRItop to SMIpl
      xnetop=XE(H05TOP)
C NEW: find topside scale height hsc=hsip(nmf2/e)-hmf2
      call topscale1(hsip,scalh,xsc)
c     
CADD+
      scalh0=scalh
C NEW++++++++++++ 1st call to obtain xnepl:
      if (hpl.gt.1336.) then
      XNEPL=0.
      XNEPL=XXE6(hpl)
      endif
c bottomside profile parameters .............................
C
1501    HMF1=HMF2
        HZ=HMF2
        HEF=HME
c
c F layer - bottomside thickness parameter B0 and shape parameters B1
c
C        B1=3.0           ! Former B1=const with BO_GUL option 
      B1=hpol(HOUR,1.9,2.6,SAX,SUX,1.,1.) ! NEW ! Diurnal B1 Formula: 1.9<B1<2.6
C New coefficient BCOEF is introduced for B0_GUL option
C deduced from 3nd order approximation to coefficients x05bot=(hmF2-h05bot/B0) => replace former B0B1 array:
      BCOEF=B1*(B1*(0.0046*B1-0.0548)+0.2546)+0.3606
C
        if(GULB0) then
          call ROGUL(SEADAY,XHI,SEAX,GRAT)
CNEW!      if(NIGHT) GRAT=0.91-HMF2/4000.
          B0CNEW=HMF2*(1.-GRAT)
      B0=B0CNEW/BCOEF        ! NEW coefficient BCOEF
        else
      B0 = B0_98(HOUR,SAX,SUX,NSEASN,RSSN,YLONGI,MODIP) ! IRI-2001
        endif
      HHALF = GRAT * HMF2
C!!!!!!! F1-REGION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      F1REG=.FALSE.
      HMF1=0.
      PNMF1=0.
      C1=0.
      IF(NIGHT) GOTO 150
        FOF1=FOF1ED(ABSMBR,RSSN,XHI)
C        IF(FOF1.LT.1.E-3) GOTO 150
        IF((FOF1.LT.1.E-3).OR.(FOF1.GE.FOF2)) GOTO 150
C
          F1REG=.TRUE.
          C1=.09+.11/DELA
          PNMF1=1.24E10*FOF1*FOF1
150   NMF1=PNMF1
C!!!!!!! PARAMETER FOR E AND VALLEY-REGION !!!!!!!!!!!!!!!!!!!!!
c E-valley: depth, width, height of deepest point (HDEEP),
c height of valley top (HEF)
c
      XDEL=XDELS(SEASON)/DELA
      DNDHBR=DNDS(SEASON)/DELA
      HDEEP=HPOL(HOUR,10.5/DELA,28.,SAX,SUX,1.,1.)
      WIDTH=HPOL(HOUR,17.8/DELA,45.+22./DELA,SAX,SUX,1.,1.)
      DEPTH=HPOL(HOUR,XDEL,81.,SAX,SUX,1.,1.)
      DLNDH=HPOL(HOUR,DNDHBR,.06,SAX,SUX,1.,1.)
      IF(DEPTH.LT.1.0) GOTO 600
        IF(NIGHT) DEPTH=-DEPTH
        CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
        IF(.NOT.EXT) GOTO 667
c        WRITE(KONSOL,650)
600   WIDTH=.0
667   HEF=HME+WIDTH
       hefold=hef
      VNER = (1. - ABS(DEPTH) / 100.) * NME
c
c Parameters below E  .............................
c
2727  NMD=XMDED(XHI,RSSN,4.0E8)
      HMD=HPOL(HOUR,81.0,88.0,SAX,SUX,1.,1.)
      F(1)=HPOL(HOUR,0.02+0.03/DELA,0.05,SAX,SUX,1.,1.)
      F(2)=HPOL(HOUR,4.6,4.5,SAX,SUX,1.,1.)
      F(3)=HPOL(HOUR,-11.5,-4.0,SAX,SUX,1.,1.)
      FP1=F(1)
      FP2=-FP1*FP1/2.0
      FP30=(-F(2)*FP2-FP1+1.0/F(2))/(F(2)*F(2))
      FP3U=(-F(3)*FP2-FP1-1.0/F(3))/(F(3)*F(3))
      HDX=HMD+F(2)
c
c indermediate region between D and E region; parameters xkk
c and d1 are found such that the function reaches hdx/xdx/dxdh
c
      X=HDX-HMD
      XDX=NMD*EXP(X*(FP1+X*(FP2+X*FP30)))
      DXDX=XDX*(FP1+X*(2.0*FP2+X*3.0*FP30))
      X=HME-HDX
      XKK=-DXDX*X/(XDX*ALOG(XDX/NME))
c
c if exponent xkk is larger than xkkmax, then xkk will be set to 
c xkkmax and d1 will be determined such that the point hdx/xdx is 
c reached; derivative is no longer continuous.
c
        xkkmax=5.
        if(xkk.gt.xkkmax) then
                xkk=xkkmax
                d1=-alog(xdx/nme)/(x**xkk)
        else
      D1=DXDX/(XDX*XKK*X**(XKK-1.0))
       endif
C
C SEARCH FOR HMF1 ..................................................
C
2726  if(LAYVER) goto 6153
C
      hmf1=0.
924   IF(.not.F1REG) GOTO 380
c omit F1 feature if nmf1*0.9 is smaller than nme
      bnmf1=0.9*nmf1
      if(nme.ge.bnmf1) goto 9427
 9245       XE2H=XE2(HEF)
      if(xe2h.gt.bnmf1) then
            hef=hef-1.
            if(hef.le.hme) goto 9427
            goto 9245
            endif      
         CALL REGFA1(HEF,HMF2,XE2H,NMF2,0.001,NMF1,XE2,SCHALT,HMF1)
        IF(.not.SCHALT) GOTO 3801
 9427 continue
C        WRITE(KONSOL,11)
        IREGFA=1
       HMF1=0.
C        NMF1=0.
C        C1=0.0
        F1REG=.FALSE.
c
c Determine E-valley parameters if HEF was changed
c
3801      continue
      if(hef.ne.hefold) then
            width=hef-hme
              IF(NIGHT) DEPTH=-DEPTH
              CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
              IF(.NOT.EXT) GOTO 380
c              if(konsol.gt.1) WRITE(KONSOL,650)
               WIDTH=.0
            hef=hme
            hefold=hef
            goto 9245
            endif
C
C SEARCH FOR HST [NE3(HST)=NME] ..........................................
C
380        IF(F1REG) then
                hf1=hmf1
                xf1=nmf1
            else
            hf1=(hmf2+hef)/2.
            xf1=xe2(hf1)
                ENDIF
       hf2=hef
      xf2=xe3_1(hf2)
      if(xf2.gt.nme) goto 3885
      CALL REGFA1(hf1,HF2,XF1,XF2,0.001,NME,XE3_1,SCHALT,HST)
      if(schalt) goto 3885
        HZ=(HST+HF1)/2.0
        D=HZ-HST
        T=D*D/(HZ-HEF-D)
        GOTO 4933
 3885   continue
C     if(konsol.gt.1) WRITE(KONSOL,100)
C 100   FORMAT(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
       HZ=(HEF+HF1)/2.
      xnehz=xe3_1(hz)
       T=(XNEHZ-NME)/(HZ-HEF)
        HST=-333.
        GOTO 4933
C
C LAY-functions for middle ionosphere
C
6153    HMF1M=165.+0.6428*XHI
        HHALF = GRAT * HMF2
        HV1R = HME + WIDTH
        HV2R = HME + HDEEP
        HHMF2 = HMF2
        CALL INILAY(NIGHT,NMF2,NMF1,NME,VNER,HHMF2,HMF1M,HME,
     &                  HV1R,HV2R,HHALF,HXL,SCL,AMP,IIQU)
C        IF(IIQU.EQ.1) WRITE(KONSOL,7733)
C 7733    FORMAT('*NE* LAY amplitudes found with 2nd choice of HXL(1).')
C        IF(IIQU.EQ.2) WRITE(KONSOL,7722)
C 7722    FORMAT('*NE* LAY amplitudes could not be found.')
C
C---------- CALCULATION OF NEUTRAL TEMPERATURE PARAMETER-------
C
4933  HTA=120.0
      HTE=3000.0
        IF(NOTEM) GOTO 240
      SEC=UT*3600.
      CALL CIRA86(DAYNR,SEC,LATI,YLONGI,HOUR,COV,TEXOS,TN120,SIGMA)
        IF(HOUR.NE.0.0) THEN
      iyz=iyear
      idz=daynr
      call ut_lt(1,utni,0.0,YLONGI,iyz,idz)
      SECNI=utni*3600.
      CALL CIRA86(DAYNR,SECNI,LATI,YLONGI,0.,COV,TEXNI,TN1NI,SIGNI)
        ELSE
      TEXNI=TEXOS
      TN1NI=TN120
      SIGNI=SIGMA
        ENDIF
      TLBDH=TEXOS-TN120
      TLBDN=TEXNI-TN1NI
C
C--------- CALCULATION OF ELECTRON TEMPERATURE PARAMETER--------
C
881   CONTINUE
C !!!!!!!!!! TE(120KM)=TN(120KM) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ATE(1)=TN120
C !!!!!!!!!! TE-MAXIMUM (JICAMARCA,ARECIBO) !!!!!!!!!!!!!!!!!!!!
      HMAXD=60.*EXP(-(MLAT/22.41)**2)+210.
      HMAXN=150.
      AHH(2)=HPOL(HOUR,HMAXD,HMAXN,SAX,SUX,1.,1.)
      TMAXD=800.*EXP(-(MLAT/33.)**2)+1500.
      TMAXN=TN(HMAXN,TEXNI,TLBDN,SIGNI)+20
      ATE(2)=HPOL(HOUR,TMAXD,TMAXN,SAX,SUX,1.,1.)
C !!!!!!!!!! TE(300,400KM)=TE-AE-C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! TE(1400,3000KM)=TE-ISIS !!!!!!!!!!!!!!!!!!!!!!!!!!! replaced by
C NEW !!!!!! Te-model by J. Titheridge above 400 km
        DIPLAT=MAGBR
      CALL TEBA(DIPLAT,HOUR,NSEASN,TEA)
      ATE(3)=TEA(1)        !Te(300km)
c      ATE(4)=TEA(2)       !Te(400km)
c      ATE(6)=TEA(3)       !Te(1400km)
c      ATE(7)=TEA(4)       !Te(3000km)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C NEW !!!!!! TE(400km),TE(600), TE(1400km), Te(3000km) USING Titheridge-2006 model
      rmon=float(month)
      do jt=1,4
      goto (701,702,703,704) jt
  701 ht=400.
      goto 705
  702 ht=600.
      goto 705
  703 ht=1400.
      goto 705
  704 ht=3000.
  705 tetitn(1)=-1.
       tetitn(2)=-1.
      tetitn(3)=-1.
      call temodel(tetitn,ht,rlt,mlat,lati,rmon,cov,xkp)
      ATE(jt+3)=TETITN(1)
      enddo
      TNAHH2=TN(AHH(2),TEXOS,TLBDH,SIGMA)
      IF(ATE(2).LT.TNAHH2) ATE(2)=TNAHH2
      STTE1=(ATE(2)-ATE(1))/(AHH(2)-AHH(1))
C !!!!!!!!!! GRADIENTS ARE CALCULATED WITH !!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! CORRECTED REGION BOUNDARIES !!!!!!!!!!!!!!!!!!!!!!
      DO 1902 I=1,6
1902  STTE(I)=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
      ATE1=ATE(1)
887   CONTINUE
C
C------------ CALCULATION OF ION TEMPERATURE PARAMETERS--------
C
C !!!!!!!!!! TI(430KM,DAY)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!!!
      XSM1=430.0
      XSM(1)=XSM1
      Z1=EXP(-0.09*MLAT)
      Z2=Z1+1.
      TID1 = 1240.0 - 1400.0 * Z1 / ( Z2 * Z2 )
      MM(2)=HPOL(HOUR,3.0,0.0,SAX,SUX,1.,1.)
C !!!!!!!!!  TI < TE   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         TED1=TEA(6)+30.
        IF(TID1.GT.TED1) TID1=TED1
C !!!!!!!!!! TI(430KM,NIGHT)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!
      Z1=ABSMLT
      Z2=Z1*(0.47+Z1*0.024)*UMR
      Z3=COS(Z2)
      TIN1=1200.0-300.0*SIGN(1.0,Z3)*SQRT(ABS(Z3))
C !!!!!!!!!! TN < TI < TE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TEN1=TEA(5)
        TNN1=TN(XSM1,TEXNI,TLBDN,SIGNI)
        IF(TEN1.LT.TNN1) TEN1=TNN1
        IF(TIN1.GT.TEN1) TIN1=TEN1
        IF(TIN1.LT.TNN1) TIN1=TNN1
C !!!!!!!!!! TI(430KM,LT) FROM STEP FUNCTION !!!!!!!!!!!!!!!!!!
        TI1=TIN1
        IF(TID1.GT.TIN1) TI1=HPOL(HOUR,TID1,TIN1,SAX,SUX,1.,1.)
C !!!!!!!!!! TANGENT ON TN DETERMINES HS !!!!!!!!!!!!!!!!!!!!!!
        TI13=TEDER(130.)
        TI50=TEDER(500.)
      CALL REGFA1(130.0,500.0,TI13,TI50,0.01,TI1,TEDER,SCHALT,HS)
      IF(SCHALT) HS=200.
      TNHS=TN(HS,TEXOS,TLBDH,SIGMA)
C#      MM(1)=DTNDH(HS,TEXOS,TLBDH,SIGMA)
         MM(1)=DTNDH(HS,TLBDH,SIGMA)
      IF(SCHALT) MM(1)=(TI1-TNHS)/(XSM1-HS)
      MXSM=2
C !!!!!!!!!! XTETI ALTITTUDE WHERE TE=TI !!!!!!!!!!!!!!!!!!!!!!
2391    XTTS=500.
        X=500.
2390    X=X+XTTS
        IF(X.GE.AHH(7)) GOTO 240
        TEX=ELTE(X)
        TIX=TI(X)
        IF(TIX.LT.TEX) GOTO 2390
        X=X-XTTS
        XTTS=XTTS/10.
        IF(XTTS.GT.0.1) GOTO 2390
        XTETI=X+XTTS*5.
C !!!!!!!!!! TI=TE ABOVE XTETI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MXSM=3
        MM(3)=STTE(6)
        XSM(2)=XTETI
        IF(XTETI.GT.AHH(6)) GOTO 240
        MXSM=4
        MM(3)=STTE(5)
        MM(4)=STTE(6)
        XSM(3)=AHH(6)
        IF(XTETI.GT.AHH(5)) GOTO 240
        MXSM=5
        DTI(1)=5.
        DTI(2)=5.
        MM(3)=STTE(4)
        MM(4)=STTE(5)
        MM(5)=STTE(6)
        XSM(3)=AHH(5)
        XSM(4)=AHH(6)
C
C CALCULATION OF ION DENSITY PARAMETER..................
C
240   IF(NOION) GOTO 141
      HNIA=100.
      HNIE=2000.
C
C INPUT OF THE ION DENSITY PARAMETER ARRAYS PF1O,PF2O AND PF3O......
C
      RIF(1)=2.
      IF(ABSLAT.LT.30.0) RIF(1)=1.
      RIF(2)=2.
      IF(COV.LT.100.0) RIF(2)=1.
      RIF(3)=SEASON
      IF(SEASON.EQ.1) RIF(3)=3.
      RIF(4)=1.
      IF(NIGHT) RIF(4)=2.
      CALL KOEFP1(PG1O)
      CALL KOEFP2(PG2O)
      CALL KOEFP3(PG3O)
      CALL SUFE(PG1O,RIF,12,PF1O)
      CALL SUFE(PG2O,RIF, 4,PF2O)
      CALL SUFE(PG3O,RIF,12,PF3O)
c
c calculate O+ profile parameters
c
      IF(ABS(XHI).LE.90.0) THEN
        ZZZ1=COS(XHI*UMR)
      ELSE
        ZZZ1=0.0
      ENDIF
      msumo=4
      RDOMAX=100.0
      MO(1)=EPSTEP(PF1O(1),PF1O(2),PF1O(3),PF1O(4),ZZZ1)
      MO(2)=EPSTEP(PF1O(5),PF1O(6),PF1O(7),PF1O(8),ZZZ1)
      MO(3)=0.0
      HO(1)=EPSTEP(PF1O(9),PF1O(10),PF1O(11),PF1O(12),ZZZ1)
      HO(2)=290.0
      IF((RIF(2).EQ.2.).AND.(RIF(3).EQ.2.)) HO(2)=237.0
      HO(4)=PF2O(1)
        ho05=pf2o(4)
      MO(4)=PF2O(2)
      MO(5)=PF2O(3)
c
c adjust gradient MO(4) of O+ profile segment above F peak
c
7100    HO(3)=(ALG100-MO(5)*(HO(4)-ho05))/MO(4)+HO(4)
        IF(HO(3).LE.HO(2)+20.) THEN
                MO(4)=MO(4)-0.001
                GOTO 7100
                endif
        hfixo=(ho(2)+ho(3))/2.
c
c find height H0O of maximum O+ relative density
c
      DELX=5.0
      X=HO(2)
      YMAXX=0.0
7102  X=X+DELX
      Y=RPID(X,HFIXO,RDOMAX,msumo,MO,DDO,HO)
      IF(Y.LE.YMAXX) then
        if(delx.le.0.1) GOTO 7104
        x=x-delx
        delx=delx/5.
      ELSE
        YMAXX=Y
      ENDIF
      GOTO 7102
7104  H0O=X-DELX/2.
7101    if(y.lt.100.0) goto 7103
          rdomax=rdomax-0.01
        y=rpid(h0o,hfixo,rdomax,msumo,mo,ddo,ho)
        goto 7101
7103    yo2h0o=100.-y
        yoh0o=y
c
c calculate parameters for O2+ profile
c
        hfixo2  = pf3o(1)
        rdo2mx = pf3o(2)
      DO 7105 L=1,2
                I = L * 2
                HO2(L)=PF3O(1+I)+PF3O(2+I)*ZZZ1
7105            MO2(L+1)=PF3O(7+I)+PF3O(8+I)*ZZZ1
      MO2(1)=PF3O(7)+PF3O(8)*ZZZ1
        if(hfixo2.gt.ho2(1)) then
           ymo2z=mo2(2)
        else
           ymo2z=mo2(1)
        endif
        aldo21=alog(rdo2mx)+ymo2z*(ho2(1)-hfixo2)
        hfixo2=(ho2(2)+ho2(1))/2.
        rdo2mx=exp(aldo21+mo2(2)*(hfixo2-ho2(1)))
c
c make sure that rd(O2+) is less or equal 100-rd(O+) at O+ maximum
c
7106  Y=RPID(H0O,hfixo2,rdo2mx,2,MO2,DO2,HO2)
      IF(Y.GT.yo2h0o) then
        MO2(3)=MO2(3)-0.02
        GOTO 7106
        endif
c
C use ratio of NO+ to O2+ density at O+ maximum to calculate
c NO+ density above the O+ maximum (H0O)
c
      IF(y.LT.1.) then
        NOBO2=0.0
      ELSE
        NOBO2= (yo2h0o-y)/y
      ENDIF
C
C CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
C
141   IF(.NOT.F1REG) HMF1=HZ
      dtec0=0.
C==== TECpl BETWEEN 1336km TO PLASMAPAUSE (=<20,200 KM)
      CN1000=xe(1336.)
c-      CN1000=xnetop
c
c TEC-plasmasphere:
c
      IF (HPL.GT.1336.) THEN
       call dinteg(HPL,1336.,CN1000,TCPL)
c-       call dinteg(HPL,h05top,CN1000,TCPL)
      else
      TCPL=0.
       XNEPL=xe(1336.)
      endif
      TCPL=TCPL/1.0E+16         ! TECU
c
c  calculate total electron content (TEC) in m-2 using the
c  stepsize selection 1 (IRI_TEC)
      call iri_tec (heibeg,1336.,1,tec,tect,tecb)
c-      call iri_tec(heibeg,h05top,1,tec,tect,tecb)
             TEC=TEC/1.0E+16
          TCB=TEC*TECB/100.0
       TCT=TEC*TECT/100.0      ! TCT of IRI-ISO in [hmF2, h05top]
      IF (jf(30)) goto 149   ! None TECI input 
C
C NEW!!! Fitting h05top to TECI input: REPLACED by change of foF2:
C
      TECR=TECI-TCB-TCPL    !  TECtop from TECI
      DTEC0=TEC+TCPL-TECI   ! BIAS between TECm-TECgps
      rattec=TECI/(TEC+TCPL) ! Ratio TECgps/TECmod
      rattop=(TECI-TCB-TCPL)/TCT   ! Ratio for Topside TEC
C
C Calculation of New h05top:
C
      dw=1.
C
C Allow correction of foF2 in one step to fit TECI:
      iflagtc=iflagtc+1
      if (iflagtc.lt.2) then
         fof2prec=FOF2
      if (iflagtc.eq.1) fof2rem=FOF2
         IF (.not.FOF2IN) THEN
      FOF2=fof2prec*sqrt(rattec)
C Add control :
      if (FOF2.lt.FOE) then
      FOF2=fof2prec
      endif
              ENDIF
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Search h05top new in the plasmasphere model
       xnm_pre=NMF2
      NMF2=1.24E10*FOF2*FOF2
       xnetop=0.5*NMF2  ! 0.5 NmF2  
      if (fof2IN) goto 1504
      if (hmf2IN) goto 1504
C
      rtop_pre=(h05top-hmF2)/hmF2
C Add recalculation of hmF2:+++++++++++++++++++++++++++++++++++
C Use model dh(dne/dnm):
      xnef=fof2*fof2
      xnefm=fof2prec*fof2prec
C
      dlogne=alog10(xnef)-alog10(xnefm)  !new
      call subdlogh(dlogne,mlat,sday,rzs,dlgh) !new
      rathhm=10.**dlgh !new
      if (rathhm.lt.0.8) rathhm=0.8
      HMF2=HMF2*rathhm !new
      AHMF2=HMF2
      HMF2IN=.TRUE.
      endif
C ADD Return to main PROGRAM with updated foF2 and hmF2
      OARR(1)=NMF2
      OARR(2)=HMF2
      IF (TECIN)  RETURN 
C *********************************************************
C
C CHANGE bottomsode HHALF according new hmF2
 1504      HHALF = GRAT * HMF2
      if ((TECIN).and.(abs(DTEC0).gt.2.)) then   !allow change of Hsc to reduce biasTECtop
      h05top_pre=h05top
            endif
C
1503  continue
C
C NEW: find topside scale height hsc=hsip(nmf2/e)-hmf2
      call topscale1(hsip,scalh,xsc)
      CN1000=xe(1336.)     !
c-     CN1000=XNETOP
c
c new TEC-plasmasphere:
c
      IF (HPL.GT.1336.) THEN
       call dinteg(HPL,1336.,CN1000,TCPL)
      TCPL=TCPL/1.0E+16         ! TECU
      else
        TCPL=0.0
       XNEPL=xe(1336.)
      endif
c
c  new calculate total electron content (TEC) in m-2 using the
c  stepsize selection 1 (IRI_TEC)
      call iri_tec (heibeg,1336.,1,tec,tect,tecb)
             TEC=TEC/1.0E+16
          TCB=TEC*TECB/100.0
       TCT=TEC*TECT/100.0      ! TCT of IRI-ISO
C
 149     kk=0
C**********TEST
      IF (.NOT.(JF(18))) GOTO 330
         H0NE=ifix(HPL/1000.+0.0005)*1000.
         HEIGHT=H0NE
 300   CONTINUE
C
 230   X=HEIGHT

       KK=KK+1
      OUTF(0,KK)=X
      NEI=XE(X)

      OUTF(1,KK)=NEI

      OUTF(10,KK)=8.98E-06*SQRT(NEI)
      OUTF(11,KK)=ALOG10(NEI)
330   IF(NOTEM) GOTO 7108
C
C+ Add Call of Temodel (J.E.Titheridge)
      tetitn(1)=-1.
      tetitn(2)=-1.
      tetitn(3)=-1.

      if (height.ge.400.) then
      call temodel(tetitn,height,rlt,mlat,lati,rmon,cov,xkp)
      outf(4,kk)=tetitn(1)
      outf(3,kk)=tetitn(2)
      outf(2,kk)=abs(tetitn(3))
      kk400=kk
      endif
      kkti=kk400
C      IF((HEIGHT.GT.HTE).OR.(HEIGHT.LT.HTA)) GOTO 7108
      IF((HEIGHT.GE.400.).OR.(HEIGHT.LT.HTA)) GOTO 7108
        TNH=TN(HEIGHT,TEXOS,TLBDH,SIGMA)
        TIH=TNH
        IF(HEIGHT.GE.HS) TIH=TI(HEIGHT)
        TEH=ELTE(HEIGHT)
        IF(TIH.LT.TNH) TIH=TNH
        IF(TEH.LT.TIH) TEH=TIH
C !!! ADD CHECK OF Ti(h<400km)<ti(400km) BELOW 400 km:
      if (abs(height-200.).lt.1.0E-05) kk200=kk
      if (tih.gt.outf(3,kk400)) then
c   tih=teh-outf(4,kk-1)+outf(3,kk-1)
      tih=-tih
      kkti=kk+1
      endif
C     
        OUTF(2,kk)=TNH
        OUTF(3,kk)=TIH
        OUTF(4,kk)=TEH
7108  IF(NOION) GOTO 7118
      IF((HEIGHT.GT.HNIE).OR.(HEIGHT.LT.HNIA)) GOTO 7118
        if(DY) then
      call IONCOM(HEIGHT,XHI*UMR,LATI*UMR,COV,ZMONTH,DION)
      ROX=DION(1)
      RHX=DION(2)
      RNX=DION(3)
      RHEX=DION(4)
      RNOX=DION(5)
      RO2X=DION(6)
      RCLUST=DION(7)
        else
      ROX=RPID(HEIGHT,HFIXO,RDOMAX,msumo,MO,DDO,HO)
      RO2X=RPID(HEIGHT,HFIXO2,rdo2mx,2,MO2,DO2,HO2)
      CALL RDHHE(HEIGHT,H0O,ROX,RO2X,NOBO2,10.,RHX,RHEX)
      RNOX=RDNO(HEIGHT,H0O,RO2X,ROX,NOBO2)
      RNX=-1.
      RCLUST=-1.
      endif

c ion densities are given in percent of total electron density (xnorm=1.);
c to get ion densities in cm-3 use the following statement
        elede=outf(1,kk)
        xnorm=elede/100.
c      xnorm=1.
      OUTF(5,kk)=ROX*xnorm
      OUTF(6,kk)=RHX*xnorm
      OUTF(7,kk)=RHEX*xnorm
      OUTF(8,kk)=RO2X*xnorm
      OUTF(9,kk)=RNOX*xnorm
C      OUTF(10,kk)=RNX*xnorm
C      OUTF(11,kk)=RCLUST*xnorm

7118  CONTINUE
7777  if (HEIGHT.ge.10000.) sh=2000.
      if ((HEIGHT.gt.3000.).and.(HEIGHT.le.10000.)) sh=1000.
      if ((HEIGHT.gt.2000.).and.(HEIGHT.le.3000.)) sh=500.
      if ((HEIGHT.gt.1000.).and.(HEIGHT.le.2000.)) sh=200.
      if ((HEIGHT.gt.500.).and.(HEIGHT.le.1000.)) sh=50.
      if ((HEIGHT.le.500.)) sh=20.
      HEIGHT=HEIGHT-SH
C     KK <=120 - LIMIT NUMBER OF HEIGHTS FOR OUTPUT ARRAY OUTF(11,120)
      if ((HEIGHT.GE.heibeg).and.(KK.le.120)) goto 230
C
C Number of heights:
C
      NUMHEI=KK
C Add output line for hpl, Nepl:
      outf(0,0)=hpl
      outf(1,0)=XXE6(hpl)
      OUTF(10,0)=8.98E-06*SQRT(xnepl)
      outf(11,0)=ALOG10(xnepl)
C+ Add Call of Temodel (J.E.Titheridge)
      tetitn(1)=-1.
      tetitn(2)=-1.
      tetitn(3)=-1.
      call temodel(tetitn,hpl,rlt,mlat,lati,rmon,cov,xkp)
      outf(4,0)=tetitn(1)
      outf(3,0)=tetitn(2)
      outf(2,0)=abs(tetitn(3))
C Check Ti(h<400km) for interpolation:
      if (kkti.gt.kk-2) then
      kkti=kk-2
      outf(3,kkti)=outf(3,kk400)-100.
      endif
      if (outf(0,kkti).gt.200.) kkti=kk200   
C   Fit Ti-profile at [200km:400km] to Ti(kkti):
      do i=kk400+1,kk
      if ((i.lt.kkti).and.(i.gt.kk400)) then
       outf(3,i)=outf(3,kkti)+(outf(3,kk400)-outf(3,kkti))/
     &(outf(0,kk400)-outf(0,kkti))*(outf(0,i)-outf(0,kkti))
            endif
      if (outf(2,i).gt.outf(3,i)) then 
       outf(2,i)=outf(3,i)/outf(3,i-1)*outf(2,i-1)  ! correct: Tn < Ti
       endif
      enddo
C
C ADDITIONAL PARAMETER FIELD OARR
C
5555     IF(NODEN) GOTO 6192
      OARR(1)=NMF2
      OARR(2)=HMF2
      OARR(3)=NMF1
      OARR(4)=HMF1
      OARR(5)=NME
      OARR(6)=HME
      OARR(7)=NMD
      OARR(8)=HMD
      OARR(9)=HHALF
      OARR(10)=B0
      OARR(11)=VNER
      OARR(12)=HEF
6192    IF(NOTEM) GOTO 6092
      OARR(13)=ATE(2)
      OARR(14)=AHH(2)
      OARR(15)=ATE(3)
      OARR(16)=ATE(4)
      OARR(17)=ATE(5)
      OARR(18)=ATE(6)
      OARR(19)=ATE(7)
      OARR(20)=ATE(1)
      OARR(21)=TI1
      OARR(22)=XTETI
6092  OARR(23)=XHI
      OARR(24)=SUNDEC
      OARR(25)=DIP
      OARR(26)=MAGBR
      OARR(27)=MODIP
      OARR(28)=DELA
      OARR(29)=SAX
      OARR(30)=SUX
      OARR(31)=SEASON
      OARR(32)=RL        ! L value
      OARR(33)=RSSN
      OARR(34)=COV
      OARR(35)=B1
      OARR(36)=xm3000      !
      OARR(37)=LATI
      OARR(38)=YLONGI
      OARR(39)=gind     !
      OARR(40)=MLONG     !
      OARR(41)=f107d     !
      OARR(42)=MLAT     !
      OARR(43)=daynr    !
      OARR(44)=gmlt     !
      OARR(45)=stormcor  !
      OARR(46)=scalh    !
      OARR(47)=xsc      !
CADD+
      OARR(48)=scalh0    ! 
      OARR(49)=DTEC0      !
      OARR(50)=fof2rem          !
3330  CONTINUE
      RETURN
      END
