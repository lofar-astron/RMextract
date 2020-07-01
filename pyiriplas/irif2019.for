C IRIF2019.FOR________________________May 2019 (T.L.Gulyaeva)
C
C 05.05.2019 subroutine tcongec corrected for dimension of 
C       real ionoindx(890),indrz(890)
C
C IRIF2017.FOR________________________May 2017 (T.L.Gulyaeva)
C
C IRIF2015c.FOR________________________Feb 2016 (T.L.Gulyaeva)
C
C IRIF2015.FOR________________________May 2015 (T.L.Gulyaeva)
C
C IRIF2014.FOR________________________Apr 2014 (T.L.Gulyaeva)
C
C IRIF2013.FOR________________________Oct 2013 (T.L.Gulyaeva)
C
C IRIF2012.FOR________________________Apr 2012 (T.L.Gulyaeva)
C
C IRIF2011.FOR________________________Dec 2011 (T.L.Gulyaeva)
C
C IRIF2007.FOR________________________Mar 2008 (T.L.Gulyaeva)
C 
C Changes of SMI plasmasphere model: one function in [h05top:hpl]
C
C IRIF2006.FOR________________________May 2006 (T.L.Gulyaeva)
C
C Changes from IRIF15.FOR                          July 2003
C
C Topside IRI-Bent model is corrected with QF-factor by fitting
C Ne(h) profile to h05top height at Ne05=0.5*NmF2 deduced from
C ISIS and IK19 topside electron density profiles
C
C IRIF15.FOR
C**************************************************************
C changes from IRIF13 to IRIF14:                   Apr. 2001
C        XE(H) INCLUDING XXE6(H)
C        XXE6(H) NE(H) PROFILE ABOVE h05top km
C        ABLON PROCEDURE USED BY XXE(6)
C**************************************************************
C********** INTERNATIONAL REFERENCE IONOSPHERE ****************
C**************************************************************
C****************  FUNCTIONS,SUBROUTINES  *********************
C**************************************************************
C** NE:         XE1,DXE1N,XE2,XE3,XE4,XE5,XE6,XE
C** TE/TI:      TEBA,SPHARM,ELTE,TEDE,TI,TEDER,TN,DTNDH
C** NI:         RPID,RDHHE,RDNO,KOEFP1,KOEFP2,KOEFP3,SUFE
C** PEAKS:      FOUT,XMOUT,HMF2EDS,FOF1ED,FOEEDI,XMDED,GAMMA1
C** MAG. FIELD: GGM,FIELDG
C** FUNCTIONS:  REGFA1,TAL
C** TIME:       SOCO,HPOL,MODA,UT_LT
C** INTERPOL.:  B0POL,B0_TAB
C** EPSTEIN:    RLAY,D1LAY,D2LAY,EPTR,EPST,EPSTEP,EPLA
C** LAY:        XE2TO5,XEN,ROGUL,VALGUL,LNGLSN,LSKNM,INILAY
C** NI-new:     IONCOM,IONCO1,IONCO2,RPDA
C** INDICES:    TCON, TCONGEC,TCONIND,SUBPROXY              !*
C** Updating:   LSTID
C**************************************************************
C
C**************************************************************
C***  -------------------ADDRESSES------------------------  ***
C***  I  PROF. K. RAWER             DR. D. BILITZA       I  ***
C***  I  HERRENSTR. 43              GSFC CODE 933        I  ***
C***  I  7801 MARCH 1               GREENBELT MD 20771   I  ***
C***  I  F.R.G.                     USA                  I  ***
C***  ----------------------------------------------------  ***
C***  I              DR. T.L. GULYAEVA                   I  ***
C***  I              IZMIRAN, Troitsk,                   I  ***
C***  I              10848 MOSCOW                        I  ***
C***  I              RUSSIA                              I  ***
C***  I              gulyaeva@izmiran.ru                 I  ***
C***  ----------------------------------------------------  ***
C**************************************************************
C*************** ELECTRON DENSITY ****************************
C*************************************************************
C
C
      FUNCTION XE1(H)
c----------------------------------------------------------------
C REPRESENTING ELECTRON DENSITY(M-3) IN THE TOPSIDE IONOSPHERE
C (H=HMF2....h05top) BY HARMONIZED BENT-MODEL ADMITTING
C VARIABILITY OFGLOBAL PARAMETER ETA,ZETA,BETA,DELTA WITH
C GEOM. LATITUDE, SMOOTHED SOLAR FLUX AND CRITICAL FREQUENCY
C (SEE MAIN PROGRAM).
C [REF.:K.RAWER,S.RAMAKRISHNAN,1978]
c----------------------------------------------------------------
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1,F1REG
     &          /BLO10/         BETA,ETA,DELTA,ZETA
     &          /ARGEXP/        ARGMAX
C# new common block:
     &          /QTOP/Y05,H05TOP,QF,XNETOP,xm3000,hhalf,tau,HTOP,QF0
      logical    f1reg

      DXDH = (1000.-HMF2)/700.
      x0 = 300. - delta
      xmx0 = (H-HMF2)/DXDH
      x = xmx0 + x0
      eptr1 = eptr(x,beta,394.5) - eptr(x0,beta,394.5)
      eptr2 = eptr(x,100.,300.0) - eptr(x0,100.,300.0)

      yy = BETA * ETA * eptr1 + ZETA * (100. * eptr2 - xmx0)

      Y = yy * dxdh
      if(abs(Y).gt.argmax) Y = sign(argmax,Y)
C++
      if (y.lt.0.) y=abs(y) ! extra check ++++++++++
C++
C# NEW CORRECTING FACTOR QFAC:
      IF((QF.EQ.1.).AND.(ABS(H-H05TOP).LT.1.)) then 
      if (y.gt.0.) QF=Y05/Y
      if (QF.gt.50.0) QF=50.           !++++++++++
      endif
C ADD :

C# Former IRI-Bent expression for Ne :
C# XE1 = XNMF2 * EXP(-Y)

C# IS REPLACED BY CORRECTING FACTOR QFAC (T.L.GULYAEVA, 2003):
      XE1 = XNMF2 * EXP(-Y*QF)
      RETURN
      END
C
C
      FUNCTION DXE1N(H)
C LOGARITHMIC DERIVATIVE OF FUNCTION XE1 (KM-1).
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1,F1REG
     &          /BLO10/         BETA,ETA,DELTA,ZETA
      logical    f1reg
      x0 = 300. - delta
      X=(H-HMF2)/(1000.0-HMF2)*700.0 + x0
      epst2 = epst(x,100.0,300.0)
      epst1 = epst(x,beta ,394.5)
      DXE1N = - ETA * epst1 + ZETA * (1. - epst2)
      RETURN
      END
C
C
      REAL FUNCTION XE2(H)
C ELECTRON DENSITY FOR THE BOTTOMSIDE F-REGION (HMF1...HMF2).
      COMMON    /BLOCK1/HMF2,XNMF2,HMF1,F1REG
     &          /BLOCK2/B0,B1,C1        /ARGEXP/ARGMAX
      logical    f1reg
      X=(HMF2-H)/B0
      z=x**b1
      if(z.gt.argmax) z=argmax
      XE2=XNMF2*EXP(-z)/COSH(X)

      RETURN
      END
C
C
      REAL FUNCTION XE3(H)
C ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1).
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1,F1REG
     &          /BLOCK2/        B0,B1,C1
      logical    f1reg
      XE3=XE2(H)+XNMF2*C1*SQRT(ABS(HMF1-H)/B0)
      RETURN
      END
C
        REAL FUNCTION XE3_1(H)
C ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1)
C USING THE NEW DEFINED F1-LAYER FUNCTION (Reinisch and Huang, Advances 
C in Space Research, Volume 25, Number 1, 81-88, 2000)
        COMMON   /BLOCK1/ HMF2,XNMF2,HMF1,F1REG
     &     /BLOCK2/ B0,B1,D1F1
      logical f1reg
C
      h1bar=h
        if (f1reg) H1BAR=HMF1*(1.0-((HMF1-H)/HMF1)**(1.0+D1F1))
        XE3_1=XE2(H1BAR)
        RETURN
        END
C
C
        REAL FUNCTION XE4_1(H)
C ELECTRON DENSITY FOR THE INTERMEDIATE REGION (HEF...HZ)
C USING THE NEW DEFINED FUNCTION
        COMMON   /BLOCK3/ HZ,T,HST,STR
     &     /BLOCK4/ HME,XNME,HEF
C
      if(hst.lt.0.0) then
      xe4_1=xnme+t*(h-hef)
      return
      endif
        IF(HST.EQ.HEF) THEN
           H1BAR=H
        ELSE
           H1BAR=HZ+0.5*T-SIGN(1.0,T)*SQRT(T*(0.25*T+HZ-H))
        ENDIF
        XE4_1=XE3_1(H1BAR)
        RETURN
        END

C
      REAL FUNCTION XE4(H)
C ELECTRON DENSITY FOR THE INDERMEDIUM REGION (HEF..HZ).
      COMMON    /BLOCK3/        HZ,T,HST,STR
     &          /BLOCK4/        HME,XNME,HEF
      IF(HST.LT.0.) GOTO 100
      XE4=XE3(HZ+T/2.0-SIGN(1.0,T)*SQRT(T*(HZ-H+T/4.0)))
      RETURN
100   XE4=XNME+T*(H-HEF)
      RETURN
      END
C
C
      REAL FUNCTION XE5(H)
C ELECTRON DENSITY FOR THE E AND VALLEY REGION (HME..HEF).
      LOGICAL NIGHT
      COMMON    /BLOCK4/        HME,XNME,HEF
     &          /BLOCK5/        NIGHT,E(4)
      T3=H-HME
      T1=T3*T3*(E(1)+T3*(E(2)+T3*(E(3)+T3*E(4))))
      IF(NIGHT) GOTO 100
      XE5=XNME*(1+T1)
      RETURN
100   XE5=XNME*EXP(T1)
      RETURN
      END
C
C
      REAL FUNCTION XE6(H)
C ELECTRON DENSITY FOR THE D REGION (HA...HME).
      COMMON    /BLOCK4/        HME,XNME,HEF
     &          /BLOCK6/        HMD,XNMD,HDX
     &          /BLOCK7/        D1,XKK,FP30,FP3U,FP1,FP2
      IF(H.GT.HDX) GOTO 100
      Z=H-HMD
      FP3=FP3U
      IF(Z.GT.0.0) FP3=FP30
      XE6=XNMD*EXP(Z*(FP1+Z*(FP2+Z*FP3)))
      RETURN
100   Z=HME-H
      XE6=XNME*EXP(-D1*Z**XKK)
      RETURN
      END
C
C
      REAL FUNCTION XE(H)
C ELECTRON DENSITY BETWEEN HA(KM) AND 1336 KM
C SUMMARIZING PROCEDURES  NE1....6; 
C SMI XXE6 PROCEDURE INCLUDED FOR 1336 TO 20,000 km
      DIMENSION OARR(50)
      COMMON    /BLOCK1/HMF2,XNMF2,XHMF1,f1reg
     &          /BLOCK3/HZ,T,HST,STR
     &          /BLOCK4/HME,XNME,HEF
      COMMON /DEM/W,GLON,HOUR,DAYNR,CN1000,ALON,BLON,ENRE,HSC
     &          /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU,HTOP,QF0
     &/PLAS/HPL,XNEPL,XKP,ICALLS,RUT,RLT,OARR,TECI,TCB,TCT,TCPL,TEC

      logical    f1reg
      integer daynr
      EXTERNAL ABLON
      if(f1reg) then
        hmf1=xhmf1
      else
        hmf1=hmf2
      endif
      hup=hpl
      hs=HTOP
      if (H.GT.hs) GOTO 10 ! Start of plasmasphere
C
      IF(H.LT.HMF2) GOTO 100
      XE=XE1(H)
      RETURN
100   IF(H.LT.HMF1) GOTO 300
      XE=XE2(H)
      RETURN
300   IF(H.LT.HZ) GOTO 400
      XE=XE3_1(H)
      RETURN
400   IF(H.LT.HEF) GOTO 500
      XE=XE4_1(H)
      RETURN
500   IF(H.LT.HME) GOTO 600
      XE=XE5(H)
      RETURN
600   XE=XE6(H)
      RETURN
   10 XE=XXE6(H)
      RETURN
      END
C
C**********************************************************                     
      FUNCTION XXE6(X)
C==== SMI PLASMASPHERE ELECTRON DENSITY PROFILE ABOVE h05top
C
C T.Gulyaeva, July 2004; Corrected: August 2014.
C !! Formulation is changed including Topside Boundary of ionosphere:
C Height HEITOP, Electron Density ELTOP.
C=NEW= SMI PLASMASPHERE ELECTRON DENSITY PROFILE ABOVE HEITOP km
C     HOUR:  Local Time
C
      INTEGER  ND,DAYNR
      real N6370,L,lt
      COMMON /DEM/W,GLON,HOUR,DAYNR,CN1000,ALON,BLON,enre,hsc 
     &   /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU,HTOP,QF0

c     HEITOP,ELTOP,ABSMLT
      PARAMETER ( PI=3.14159265, RE=6370., RAD=PI/180.)
CREM  heitop=htop
      heitop=H05TOP

      LT=HOUR
      ALT=X
C New:
       HHUP=HEITOP
CREM  ELTOP=XE(HHUP)
      ELTOP=XNETOP
      ND=DAYNR
      L= 1.0 + ALT/RE
c
c>>   A1 = 1.-0.1*COS(PI*(LT-4.)/12.)   !Old
C A1 new - corrected A1(L,LT)
      xn = 0.9                 ! night value
      xd = 1.-0.1*cos(PI/L**5)    ! noon value
      A1 = xn+(xd-xn)/(1.+exp(11.-LT))+(xn-xd)/(1.+exp(21.-LT))  ! LT function
         A2 = 1.+0.2*SQRT(W)
      A3 = 1.+0.001*W
      AL = 1.-EXP(-0.04*L**4)
      AL2 = 1.-EXP(-0.04*2**4)
C
C !!New!!    CORRECTION: C1=3E+09 AS UPPER LIMIT OF COEFFICIENT
C Old    C1=3.0E+09
C New expression for C1(LT):
      c1=(1.0+cos(-pi*(LT-12.)/24.)/2.)*1.0E+09

      if (alon.le.360.) goto 10
      A0 = (1.0+0.7*COS(RAD*(GLON+70.)))*(1.+COS(2*PI*(ND+16)/365.25))
      BL = 0.5*A0*AL-0.7*L*EXP(0.1*A0*AL)
      BL2 = 0.5*A0*AL2-1.4*EXP(0.1*A0*AL2)
C-      N6370 = 3.0E9*A1*A2*EXP(A3*BL2)
      N6370 = c1*A1*A2*EXP(A3*BL2)
C+!!! old formula:      H = 1342.5/ALOG(CN1000/N6370) IS REPLACED BY :
      H=(RE-HHUP)*.25/ALOG(ELTOP/N6370)
      if (ALT.le.6370.) then
C!!! old formula:    DN = CN1000*EXP((1000.-ALT)/H/L/L) IS REPLACED BY :
      DN=ELTOP*EXP((HHUP-ALT)/H/L/L)
      else
C!!! old formula:    DN = 3.0E9*A1*A2*EXP(A3*BL) IS REPLACED BY :
      DN=C1*A1*A2*EXP(A3*BL)
      endif
      XXE6=DN
      RETURN
   10 a0a= (1.0+0.7*COS(RAD*(ALON+70.)))*(1.+COS(2*PI*(ND+16)/365.25))
      b0b= (1.0+0.7*COS(RAD*(BLON+70.)))*(1.+COS(2*PI*(ND+16)/365.25))
      bla= 0.5*a0a*AL-0.7*L*EXP(0.1*a0a*AL)
      blb= 0.5*b0b*AL-0.7*L*EXP(0.1*b0b*AL)
      bl2a= 0.5*a0a*AL2-0.7*2*EXP(0.1*a0a*AL2)
      bl2b= 0.5*b0b*AL2-0.7*2*EXP(0.1*b0b*AL2)
C!!! Old formula      drea = 3.0E9*A1*A2*EXP(A3*bl2a)
      drea = C1*A1*A2*EXP(A3*bl2a)
C!!! Old formula      dreb = 3.0E9*A1*A2*EXP(A3*bl2b)
      dreb = C1*A1*A2*EXP(A3*bl2b)
      drea=alog10(drea)
      dreb=alog10(dreb)
      delon=blon-alon
      if (delon.le.0.) delon=delon+360.
      deglon=glon-alon
      if (deglon.le.0.) deglon=deglon+360.
      dne=drea+(dreb-drea)*deglon/delon
      N6370=10.**dne
C!!!      H = 1342.5/ALOG(CN1000/N6370)
      H=(RE-HHUP)*.25/ALOG(ELTOP/N6370)         ! 0.25 instead of 0.25
      if (ALT.le.6370.) then
C!!!         DN = CN1000*EXP((1000.-ALT)/H/L/L)
      DN=ELTOP*EXP((HHUP-ALT)/H/L/L)
      else
C!!!         dna = 3.0E9*A1*A2*EXP(A3*bla)
            dna = C1*A1*A2*EXP(A3*bla)
         dna=alog10(dna)
C!!!         dnb = 3.0E9*A1*A2*EXP(A3*blb)
         dnb = C1*A1*A2*EXP(A3*blb)
         dnb=alog10(dnb)
         dnx=dna+(dnb-dna)*deglon/delon
         dn=10.**dnx
      endif
      XXE6=DN
      RETURN
      END
C
      SUBROUTINE ablon(xlati,xlongi,alonn,blonn)
C==== SMI GRID POINTS FOR XXE6(H)       
      DIMENSION xlat(37),xlon(37),zlon(6)
      common /ab/ fimeq,zlon,nl
      DATA xlat /10.62,10.43,10.06,9.60,9.10,8.76,8.60,8.74,9.29,9.82,
     &9.83,9.27,8.38,7.48,7.07,7.05,6.25,4.07,1.76,0.36,-0.42,-1.37,
     &-2.84,-4.01,-4.15,-4.37,-6.34,-9.73,-12.75,-14.46,-14.12,-10.41,
     &-3.90,2.66,7.44,9.95,10.62/
      DATA xlon /0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,
     &130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,
     &260.,270.,280.,290.,300.,310.,320.,330.,340.,350.,360./
      do 1 i=1,36
      k=i
      if ((xlongi.ge.xlon(i)).and.(xlongi.lt.xlon(i+1))) goto 2
    1 CONTINUE
    2 fimeq=xlat(k)+(xlat(k+1)-xlat(k))*(xlongi-xlon(k))*0.1
      if ((xlongi.lt.178.).or.(xlongi.gt.330.)) goto 3
      if (xlati.le.(fimeq+1.)) goto 11
      if (xlati.lt.2.6) goto 4
      alonn=178.
      blonn=330.
      goto 11
    4 nl=0
      do 5 i=1,36
c      k=i
      if (((xlati.ge.xlat(i)).and.(xlati.lt.xlat(i+1))).or.
     &((xlati.le.xlat(i)).and.(xlati.gt.xlat(i+1)))) goto 15
      goto 5
   15 nl=nl+1
      zlon(nl)=xlon(i)+10.*(xlati-xlat(i))/(xlat(i+1)-xlat(i))
    5 continue
      alonn=zlon(1)
      blonn=zlon(2)
      goto 11
    3 if (xlati.ge.(fimeq-1.)) goto 11
      if (xlati.gt.7.6) goto 6
      alonn=340.
      blonn=130.
      goto 11
    6 nl=0
      do 7 i=1,36
c      k=i
      if (((xlati.ge.xlat(i)).and.(xlati.lt.xlat(i+1))).or.
     &((xlati.le.xlat(i)).and.(xlati.gt.xlat(i+1)))) goto 17
      goto 7
  17  nl=nl+1
      zlon(nl)=xlon(i)+10.*(xlati-xlat(i))/(xlat(i+1)-xlat(i))
    7 continue
      if ((xlongi.ge.zlon(1)).and.(xlongi.le.130.)) goto 8
      alonn=zlon(nl)
      blonn=zlon(1)
      goto 11
    8 continue
      do 9 i=1,nl
      k=i
      if ((xlongi.ge.zlon(i)).and.(xlongi.le.zlon(i+1))) goto 10
    9 continue
   10 alonn=zlon(k)
      blonn=zlon(k+1)
   11 RETURN
      END
C
C
C**********************************************************
C***************** ELECTRON TEMPERATURE ********************
C**********************************************************
C
      SUBROUTINE TEBA(DIPL,SLT,NS,TE)
C CALCULATES ELECTRON TEMPERATURES TE(1) TO TE(4) AT ALTITUDES
C 300, 400, 1400 AND 3000 KM FOR DIP-LATITUDE DIPL/DEG AND
C LOCAL SOLAR TIME SLT/H USING THE BRACE-THEIS-MODELS (J. ATMOS.
C TERR. PHYS. 43, 1317, 1981); NS IS SEASON IN NORTHERN
C HEMISOHERE: IS=1 SPRING, IS=2 SUMMER ....
C ALSO CALCULATED ARE THE TEMPERATURES AT 400 KM ALTITUDE FOR
C MIDNIGHT (TE(5)) AND NOON (TE(6)).
      DIMENSION C(4,2,81),A(82),TE(6)
      COMMON/CONST/UMR,PI
      DATA (C(1,1,J),J=1,81)/
     &.3100E1,-.3215E-2,.2440E+0,-.4613E-3,-.1711E-1,.2605E-1,
     &-.9546E-1,.1794E-1,.1270E-1,.2791E-1,.1536E-1,-.6629E-2,
     &-.3616E-2,.1229E-1,.4147E-3,.1447E-2,-.4453E-3,-.1853,
     &-.1245E-1,-.3675E-1,.4965E-2,.5460E-2,.8117E-2,-.1002E-1,
     &.5466E-3,-.3087E-1,-.3435E-2,-.1107E-3,.2199E-2,.4115E-3,
     &.6061E-3,.2916E-3,-.6584E-1,.4729E-2,-.1523E-2,.6689E-3,
     &.1031E-2,.5398E-3,-.1924E-2,-.4565E-1,.7244E-2,-.8543E-4,
     &.1052E-2,-.6696E-3,-.7492E-3,.4405E-1,.3047E-2,.2858E-2,
     &-.1465E-3,.1195E-2,-.1024E-3,.4582E-1,.8749E-3,.3011E-3,
     &.4473E-3,-.2782E-3,.4911E-1,-.1016E-1,.27E-2,-.9304E-3,
     &-.1202E-2,.2210E-1,.2566E-2,-.122E-3,.3987E-3,-.5744E-1,
     &.4408E-2,-.3497E-2,.83E-3,-.3536E-1,-.8813E-2,.2423E-2,
     &-.2994E-1,-.1929E-2,-.5268E-3,-.2228E-1,.3385E-2,
     &.413E-1,.4876E-2,.2692E-1,.1684E-2/
      DATA (C(1,2,J),J=1,81)/.313654E1,.6796E-2,.181413,.8564E-1,
     &-.32856E-1,-.3508E-2,-.1438E-1,-.2454E-1,.2745E-2,.5284E-1,
     &.1136E-1,-.1956E-1,-.5805E-2,.2801E-2,-.1211E-2,.4127E-2,
     &.2909E-2,-.25751,-.37915E-2,-.136E-1,-.13225E-1,.1202E-1,
     &.1256E-1,-.12165E-1,.1326E-1,-.7123E-1,.5793E-3,.1537E-2,
     &.6914E-2,-.4173E-2,.1052E-3,-.5765E-3,-.4041E-1,-.1752E-2,
     &-.542E-2,-.684E-2,.8921E-3,-.2228E-2,.1428E-2,.6635E-2,-.48045E-2,
     &-.1659E-2,-.9341E-3,.223E-3,-.9995E-3,.4285E-1,-.5211E-3,
     &-.3293E-2,.179E-2,.6435E-3,-.1891E-3,.3844E-1,.359E-2,-.8139E-3,
     &-.1996E-2,.2398E-3,.2938E-1,.761E-2,.347655E-2,.1707E-2,.2769E-3,
     &-.157E-1,.983E-3,-.6532E-3,.929E-4,-.2506E-1,.4681E-2,.1461E-2,
     &-.3757E-5,-.9728E-2,.2315E-2,.6377E-3,-.1705E-1,.2767E-2,
     &-.6992E-3,-.115E-1,-.1644E-2,.3355E-2,-.4326E-2,.2035E-1,.2985E-1/
      DATA (C(2,1,J),J=1,81)/.3136E1,.6498E-2,.2289,.1859E-1,-.3328E-1,
     &-.4889E-2,-.3054E-1,-.1773E-1,-.1728E-1,.6555E-1,.1775E-1,
     &-.2488E-1,-.9498E-2,.1493E-1,.281E-2,.2406E-2,.5436E-2,-.2115,
     &.7007E-2,-.5129E-1,-.7327E-2,.2402E-1,.4772E-2,-.7374E-2,
     &-.3835E-3,-.5013E-1,.2866E-2,.2216E-2,.2412E-3,.2094E-2,.122E-2
     &,-.1703E-3,-.1082,-.4992E-2,-.4065E-2,.3615E-2,-.2738E-2,
     &-.7177E-3,.2173E-3,-.4373E-1,-.375E-2,.5507E-2,-.1567E-2,
     &-.1458E-2,-.7397E-3,.7903E-1,.4131E-2,.3714E-2,.1073E-2,
     &-.8991E-3,.2976E-3,.2623E-1,.2344E-2,.5608E-3,.4124E-3,.1509E-3,
     &.5103E-1,.345E-2,.1283E-2,.7238E-3,-.3464E-4,.1663E-1,-.1644E-2,
     &-.71E-3,.5281E-3,-.2729E-1,.3556E-2,-.3391E-2,-.1787E-3,.2154E-2,
     &.6476E-2,-.8282E-3,-.2361E-1,.9557E-3,.3205E-3,-.2301E-1,
     &-.854E-3,-.1126E-1,-.2323E-2,-.8582E-2,.2683E-1/
      DATA (C(2,2,J),J=1,81)/.3144E1,.8571E-2,.2539,.6937E-1,-.1667E-1,
     &.2249E-1,-.4162E-1,.1201E-1,.2435E-1,.5232E-1,.2521E-1,-.199E-1,
     &-.7671E-2,.1264E-1,-.1551E-2,-.1928E-2,.3652E-2,-.2019,.5697E-2,
     &-.3159E-1,-.1451E-1,.2868E-1,.1377E-1,-.4383E-2,.1172E-1,
     &-.5683E-1,.3593E-2,.3571E-2,.3282E-2,.1732E-2,-.4921E-3,-.1165E-2
     &,-.1066,-.1892E-1,.357E-2,-.8631E-3,-.1876E-2,-.8414E-4,.2356E-2,
     &-.4259E-1,-.322E-2,.4641E-2,.6223E-3,-.168E-2,-.1243E-3,.7393E-1,
     &-.3143E-2,-.2362E-2,.1235E-2,-.1551E-2,.2099E-3,.2299E-1,.5301E-2
     &,-.4306E-2,-.1303E-2,.7687E-5,.5305E-1,.6642E-2,-.1686E-2,
     &.1048E-2,.5958E-3,.4341E-1,-.8819E-4,-.333E-3,-.2158E-3,-.4106E-1
     &,.4191E-2,.2045E-2,-.1437E-3,-.1803E-1,-.8072E-3,-.424E-3,
     &-.26E-1,-.2329E-2,.5949E-3,-.1371E-1,-.2188E-2,.1788E-1,
     &.6405E-3,.5977E-2,.1333E-1/
      DATA (C(3,1,J),J=1,81)/.3372E1,.1006E-1,.1436,.2023E-2,-.5166E-1,
     &.9606E-2,-.5596E-1,.4914E-3,-.3124E-2,-.4713E-1,-.7371E-2,
     &-.4823E-2,-.2213E-2,.6569E-2,-.1962E-3,.3309E-3,-.3908E-3,
     &-.2836,.7829E-2,.1175E-1,.9919E-3,.6589E-2,.2045E-2,-.7346E-2
     &,-.89E-3,-.347E-1,-.4977E-2,.147E-2,-.2823E-5,.6465E-3,
     &-.1448E-3,.1401E-2,-.8988E-1,-.3293E-4,-.1848E-2,.4439E-3,
     &-.1263E-2,.317E-3,-.6227E-3,.1721E-1,-.199E-2,-.4627E-3,
     &.2897E-5,-.5454E-3,.3385E-3,.8432E-1,-.1951E-2,.1487E-2,
     &.1042E-2,-.4788E-3,-.1276E-3,.2373E-1,.2409E-2,.5263E-3,
     &.1301E-2,-.4177E-3,.3974E-1,.1418E-3,-.1048E-2,-.2982E-3,
     &-.3396E-4,.131E-1,.1413E-2,-.1373E-3,.2638E-3,-.4171E-1,
     &-.5932E-3,-.7523E-3,-.6883E-3,-.2355E-1,.5695E-3,-.2219E-4,
     &-.2301E-1,-.9962E-4,-.6761E-3,.204E-2,-.5479E-3,.2591E-1,
     &-.2425E-2,.1583E-1,.9577E-2/
      DATA (C(3,2,J),J=1,81)/.3367E1,.1038E-1,.1407,.3622E-1,-.3144E-1,
     &.112E-1,-.5674E-1,.3219E-1,.1288E-2,-.5799E-1,-.4609E-2,
     &.3252E-2,-.2859E-3,.1226E-1,-.4539E-2,.1310E-2,-.5603E-3,
     &-.311,-.1268E-2,.1539E-1,.3146E-2,.7787E-2,-.143E-2,-.482E-2
     &,.2924E-2,-.9981E-1,-.7838E-2,-.1663E-3,.4769E-3,.4148E-2,
     &-.1008E-2,-.979E-3,-.9049E-1,-.2994E-2,-.6748E-2,-.9889E-3,
     &.1488E-2,-.1154E-2,-.8412E-4,-.1302E-1,-.4859E-2,-.7172E-3,
     &-.9401E-3,.9101E-3,-.1735E-3,.7055E-1,.6398E-2,-.3103E-2,
     &-.938E-3,-.4E-3,-.1165E-2,.2713E-1,-.1654E-2,.2781E-2,
     &-.5215E-5,.2258E-3,.5022E-1,.95E-2,.4147E-3,.3499E-3,
     &-.6097E-3,.4118E-1,.6556E-2,.3793E-2,-.1226E-3,-.2517E-1,
     &.1491E-3,.1075E-2,.4531E-3,-.9012E-2,.3343E-2,.3431E-2,
     &-.2519E-1,.3793E-4,.5973E-3,-.1423E-1,-.132E-2,-.6048E-2,
     &-.5005E-2,-.115E-1,.2574E-1/
      DATA (C(4,1,J),J=1,81)/.3574E1,.0,.7537E-1,.0,-.8459E-1,
     &0.,-.294E-1,0.,.4547E-1,-.5321E-1,0.,.4328E-2,0.,.6022E-2,
     &.0,-.9168E-3,.0,-.1768,.0,.294E-1,.0,.5902E-3,.0,-.9047E-2,
     &.0,-.6555E-1,.0,-.1033E-2,.0,.1674E-2,.0,.2802E-3,-.6786E-1
     &,.0,.4193E-2,.0,-.6448E-3,.0,.9277E-3,-.1634E-1,.0,-.2531E-2
     &,.0,.193E-4,.0,.528E-1,.0,.2438E-2,.0,-.5292E-3,.0,.1555E-1
     &,.0,-.3259E-2,.0,-.5998E-3,.3168E-1,.0,.2382E-2,.0,-.4078E-3
     &,.2312E-1,.0,.1481E-3,.0,-.1885E-1,.0,.1144E-2,.0,-.9952E-2
     &,.0,-.551E-3,-.202E-1,.0,-.7283E-4,-.1272E-1,.0,.2224E-2,
     &.0,-.251E-2,.2434E-1/
      DATA (C(4,2,J),J=1,81)/.3574E1,-.5639E-2,.7094E-1,
     &-.3347E-1,-.861E-1,-.2877E-1,-.3154E-1,-.2847E-2,.1235E-1,
     &-.5966E-1,-.3236E-2,.3795E-3,-.8634E-3,.3377E-2,-.1071E-3,
     &-.2151E-2,-.4057E-3,-.1783,.126E-1,.2835E-1,-.242E-2,
     &.3002E-2,-.4684E-2,-.6756E-2,-.7493E-3,-.6147E-1,-.5636E-2
     &,-.1234E-2,-.1613E-2,-.6353E-4,-.2503E-3,-.1729E-3,-.7148E-1
     &,.5326E-2,.4006E-2,.6484E-3,-.1046E-3,-.6034E-3,-.9435E-3,
     &-.2385E-2,.6853E-2,.151E-2,.1319E-2,.9049E-4,-.1999E-3,
     &.3976E-1,.2802E-2,-.103E-2,.5599E-3,-.4791E-3,-.846E-4,
     &.2683E-1,.427E-2,.5911E-3,.2987E-3,-.208E-3,.1396E-1,
     &-.1922E-2,-.1063E-2,.3803E-3,.1343E-3,.1771E-1,-.1038E-2,
     &-.4645E-3,-.2481E-3,-.2251E-1,-.29E-2,-.3977E-3,-.516E-3,
     &-.8079E-2,-.1528E-2,.306E-3,-.1582E-1,-.8536E-3,.1565E-3,
     &-.1252E-1,.2319E-3,.4311E-2,.1024E-2,.1296E-5,.179E-1/
      IF(NS.LT.3) THEN
      IS=NS
      ELSE IF(NS.GT.3) THEN
      IS=2
      DIPL=-DIPL
      ELSE
      IS=1
      ENDIF
      COLAT=UMR*(90.-DIPL)
      AZ=.2618*SLT
      CALL SPHARM(A,8,8,COLAT,AZ)
      IF(IS.EQ.2) THEN
      KEND=3
      ELSE
      KEND=4
      ENDIF
      DO 2 K=1,KEND
      STE=0.
      DO 1 I=1,81
1       STE=STE+A(I)*C(K,IS,I)
2     TE(K)=10.**STE
      IF(IS.EQ.2) THEN
      DIPL=-DIPL
      COLAT=UMR*(90.-DIPL)
      CALL SPHARM(A,8,8,COLAT,AZ)
      STE=0.
      DO 11 I=1,81
11            STE=STE+A(I)*C(4,2,I)
      TE(4)=10.**STE
      ENDIF

C---------- TEMPERATURE AT 400KM AT MIDNIGHT AND NOON
      DO 4 J=1,2
      STE=0.
      AZ=.2618*(J-1)*12.
      CALL SPHARM(A,8,8,COLAT,AZ)
      DO 3 I=1,81
3         STE=STE+A(I)*C(2,IS,I)
4       TE(J+4)=10.**STE
      RETURN
      END
C
      SUBROUTINE SPHARM(C,L,M,COLAT,AZ)
C CALCULATES THE COEFFICIENTS OF THE SPHERICAL HARMONIC
C EXPANSION THAT WAS USED FOR THE BRACE-THEIS-MODELS.
      DIMENSION C(82)
      C(1)=1.
      K=2
      X=COS(COLAT)
      C(K)=X
      K=K+1
      DO 10 I=2,L
      C(K)=((2*I-1)*X*C(K-1)-(I-1)*C(K-2))/I
10    K=K+1
      Y=SIN(COLAT)
      DO 20 MT=1,M
      CAZ=COS(MT*AZ)
      SAZ=SIN(MT*AZ)
      C(K)=Y**MT
      K=K+1
      IF(MT.EQ.L) GOTO 16
      C(K)=C(K-1)*X*(2*MT+1)
      K=K+1
      IF((MT+1).EQ.L) GOTO 16
      DO 15 I=2+MT,L
      C(K)=((2*I-1)*X*C(K-1)-(I+MT-1)*C(K-2))/(I-MT)
15    K=K+1
16    N=L-MT+1
      DO 18 I=1,N
      C(K)=C(K-N)*CAZ
      C(K-N)=C(K-N)*SAZ
18    K=K+1
20    CONTINUE
      RETURN
      END
C
C
      REAL FUNCTION ELTE(H)
c----------------------------------------------------------------
C ELECTRON TEMPERATURE PROFILE BASED ON THE TEMPERATURES AT 120
C HMAX,300,400,600,1400,3000 KM ALTITUDE. INBETWEEN CONSTANT
C GRADIENT IS ASSUMED. ARGMAX IS MAXIMUM ARGUMENT ALLOWED FOR
C EXP-FUNCTION.
c----------------------------------------------------------------
      COMMON /BLOTE/AH(7),ATE1,ST(6),D(5)
C
      SUM=ATE1+ST(1)*(H-AH(1))
      DO 1 I=1,5
      aa = eptr(h    ,d(i),ah(i+1))
      bb = eptr(ah(1),d(i),ah(i+1))
1     SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*D(I)
      ELTE=SUM
      RETURN
      END
C
C
      FUNCTION TEDE(H,DEN,COV)
C ELECTRON TEMEPERATURE MODEL AFTER BRACE,THEIS .
C FOR NEG. COV THE MEAN COV-INDEX (3 SOLAR ROT.) IS EXPECTED.
C DEN IS THE ELECTRON DENSITY IN M-3.
      Y=1051.+(17.01*H-2746.)*
     &EXP(-5.122E-4*H+(6.094E-12-3.353E-14*H)*DEN)
      ACOV=ABS(COV)
      YC=1.+(.117+2.02E-3*ACOV)/(1.+EXP(-(ACOV-102.5)/5.))
      IF(COV.LT.0.)
     &YC=1.+(.123+1.69E-3*ACOV)/(1.+EXP(-(ACOV-115.)/10.))
      TEDE=Y*YC
      RETURN
      END
C
C
C*************************************************************
C**************** ION TEMPERATURE ****************************
C*************************************************************
C
C
      REAL FUNCTION TI(H)
c----------------------------------------------------------------
C ION TEMPERATURE FOR HEIGHTS NOT GREATER 1000 KM AND NOT LESS HS
C EXPLANATION SEE FUNCTION RPID.
c----------------------------------------------------------------
      REAL              MM
      COMMON  /BLOCK8/  HS,TNHS,XSM(4),MM(5),G(4),M

      SUM=MM(1)*(H-HS)+TNHS
      DO 100 I=1,M-1
      aa = eptr(h ,g(i),xsm(i))
      bb = eptr(hs,g(i),xsm(i))
100     SUM=SUM+(MM(I+1)-MM(I))*(AA-BB)*G(I)
      TI=SUM
      RETURN
      END
C
C
      REAL FUNCTION TEDER(H)
C THIS FUNCTION ALONG WITH PROCEDURE REGFA1 ALLOWS TO FIND
C THE  HEIGHT ABOVE WHICH TN BEGINS TO BE DIFFERENT FROM TI
      COMMON    /BLOTN/XSM1,TEX,TLBD,SIG
      TNH = TN(H,TEX,TLBD,SIG)
C#      DTDX = DTNDH(H,TEX,TLBD,SIG)
      DTDX = DTNDH(H,TLBD,SIG)
      TEDER = DTDX * ( XSM1 - H ) + TNH
      RETURN
      END
C
C
      FUNCTION TN(H,TINF,TLBD,S)
C--------------------------------------------------------------------
C       Calculate Temperature for MSIS/CIRA-86 model
C--------------------------------------------------------------------
      ZG2 = ( H - 120. ) * 6476.77 / ( 6356.77 + H )
      TN = TINF - TLBD * EXP ( - S * ZG2 )
      RETURN
      END
C
C
C#      FUNCTION DTNDH(H,TINF,TLBD,S)
      FUNCTION DTNDH(H,TLBD,S)
C---------------------------------------------------------------------
      ZG1 = 6356.77 + H
      ZG2 = 6476.77 / ZG1
      ZG3 = ( H - 120. ) * ZG2
      DTNDH = - TLBD * EXP ( - S * ZG3 ) * ( S / ZG1 * ( ZG3 - ZG2 ) )
      RETURN
      END
C
C
C*************************************************************
C************* ION RELATIVE PRECENTAGE DENSITY *****************
C*************************************************************
C
C
      REAL FUNCTION RPID (H, H0, N0, M, ST, ID, XS)
c------------------------------------------------------------------
C D.BILITZA,1977,THIS ANALYTIC FUNCTION IS USED TO REPRESENT THE
C RELATIVE PRECENTAGE DENSITY OF ATOMAR AND MOLECULAR OXYGEN IONS.
C THE M+1 HEIGHT GRADIENTS ST(M+1) ARE CONNECTED WITH EPSTEIN-
C STEP-FUNCTIONS AT THE STEP HEIGHTS XS(M) WITH TRANSITION
C THICKNESSES ID(M). RPID(H0,H0,N0,....)=N0.
C ARGMAX is the highest allowed argument for EXP in your system.
c------------------------------------------------------------------
      REAL              N0
      DIMENSION         ID(4), ST(5), XS(4)
      COMMON  /ARGEXP/  ARGMAX

      SUM=(H-H0)*ST(1)
      DO 100  I=1,M
         XI=ID(I)
      aa = eptr(h ,xi,xs(i))
      bb = eptr(h0,xi,xs(i))
100           SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*XI
      IF(ABS(SUM).LT.ARGMAX) then
      SM=EXP(SUM)
      else IF(SUM.Gt.0.0) then
      SM=EXP(ARGMAX)
      else
      SM=0.0
      endif
      RPID= n0 * SM
      RETURN
      END
C
c
      SUBROUTINE RDHHE (H,HB,RDOH,RDO2H,RNO,PEHE,RDH,RDHE)
C BILITZA,FEB.82,H+ AND HE+ RELATIVE PERECENTAGE DENSITY BELOW
C 1000 KM. THE O+ AND O2+ REL. PER. DENSITIES SHOULD BE GIVEN
C (RDOH,RDO2H). HB IS THE ALTITUDE OF MAXIMAL O+ DENSITY. PEHE
C IS THE PRECENTAGE OF HE+ IONS COMPARED TO ALL LIGHT IONS.
C RNO IS THE RATIO OF NO+ TO O2+DENSITY AT H=HB.
      RDHE=0.0
      RDH=0.0
      IF(H.LE.HB) GOTO 100
      REST=100.0-RDOH-RDO2H-RNO*RDO2H
      RDH=REST*(1.-PEHE/100.)
      RDHE=REST*PEHE/100.
100   RETURN
      END
C
C
      REAL FUNCTION RDNO(H,HB,RDO2H,RDOH,RNO)
C D.BILITZA, 1978. NO+ RELATIVE PERCENTAGE DENSITY ABOVE 100KM.
C FOR MORE INFORMATION SEE SUBROUTINE RDHHE.
      IF (H.GT.HB) GOTO 200
      RDNO=100.0-RDO2H-RDOH
      RETURN
200   RDNO=RNO*RDO2H
      RETURN
      END
C
C
      SUBROUTINE  KOEFP1(PG1O)
C THIEMANN,1979,COEFFICIENTS PG1O FOR CALCULATING  O+ PROFILES
C BELOW THE F2-MAXIMUM. CHOSEN TO APPROACH DANILOV-
C SEMENOV'S COMPILATION.
      DIMENSION PG1O(80)
      REAL FELD (80)
      DATA FELD/-11.0,-11.0,4.0,-11.0,0.08018,
     &0.13027,0.04216,0.25  ,-0.00686,0.00999,
     &5.113,0.1 ,170.0,180.0,0.1175,0.15,-11.0,
     &1.0 ,2.0,-11.0,0.069,0.161,0.254,0.18,0.0161,
     &0.0216,0.03014,0.1,152.0,167.0,0.04916,
     &0.17,-11.0,2.0,2.0,-11.0,0.072,0.092,0.014,0.21,
     &0.01389,0.03863,0.05762,0.12,165.0,168.0,0.008,
     &0.258,-11.0,1.0,3.0,-11.0,0.091,0.088,
     &0.008,0.34,0.0067,0.0195,0.04,0.1,158.0,172.0,
     &0.01,0.24,-11.0,2.0,3.0, -11.0,0.083,0.102,
     &0.045,0.03,0.00127,0.01,0.05,0.09,167.0,185.0,
     &0.015,0.18/
      K=0
      DO 10 I=1,80
      K=K+1
10    PG1O(K)=FELD(I)
      RETURN
      END
C
C
      SUBROUTINE KOEFP2(PG2O)
C THIEMANN,1979,COEFFICIENTS FOR CALCULATION OF O+ PROFILES
C ABOVE THE F2-MAXIMUM (DUMBS,SPENNER:AEROS-COMPILATION)
      DIMENSION PG2O(32)
      REAL FELD(32)
      DATA FELD/1.0,-11.0,-11.0,1.0,695.0,-.000781,
     &-.00264,2177.0,1.0,-11.0,-11.0,2.0,570.0,
     &-.002,-.0052,1040.0,2.0,-11.0,-11.0,1.0,695.0,
     &-.000786,-.00165,3367.0,2.0,-11.0,-11.0,2.0,
     &575.0,-.00126,-.00524,1380.0/
      K=0
      DO 10 I=1,32
      K=K+1
10    PG2O(K)=FELD(I)
      RETURN
      END
C
C
      SUBROUTINE  KOEFP3(PG3O)
C THIEMANN,1979,COEFFICIENTS FOR CALCULATING O2+ PROFILES.
C CHOSEN AS TO APPROACH DANILOV-SEMENOV'S COMPILATION.
      DIMENSION PG3O(80)
      REAL FELD(80)
      DATA FELD/-11.0,1.0,2.0,-11.0,160.0,31.0,130.0,
     &-10.0,198.0,0.0,0.05922,-0.07983,
     &-0.00397,0.00085,-0.00313,0.0,-11.0,2.0,2.0,-11.0,
     &140.0,30.0,130.0,-10.0,
     &190.0,0.0,0.05107,-0.07964,0.00097,-0.01118,-0.02614,
     &-0.09537,
     &-11.0,1.0,3.0,-11.0,140.0,37.0,125.0,0.0,182.0,
     &0.0,0.0307,-0.04968,-0.00248,
     &-0.02451,-0.00313,0.0,-11.0,2.0,3.0,-11.0,
     &140.0,37.0,125.0,0.0,170.0,0.0,
     &0.02806,-0.04716,0.00066,-0.02763,-0.02247,-0.01919,
     &-11.0,-11.0,4.0,-11.0,140.0,45.0,136.0,-9.0,
     &181.0,-26.0,0.02994,-0.04879,
     &-0.01396,0.00089,-0.09929,0.05589/
      K=0
      DO 10 I=1,80
      K=K+1
10    PG3O(K)=FELD(I)
      RETURN
      END
C
C
      SUBROUTINE SUFE (FIELD,RFE,M,FE)
C SELECTS THE REQUIRED ION DENSITY PARAMETER SET.
C THE INPUT FIELD INCLUDES DIFFERENT SETS OF DIMENSION M EACH
C CARACTERISED BY 4 HEADER NUMBERS. RFE(4) SHOULD CONTAIN THE
C CHOSEN HEADER NUMBERS.FE(M) IS THE CORRESPONDING SET.
      DIMENSION RFE(4),FE(12),FIELD(80),EFE(4)
      K=0
100   DO 101 I=1,4
      K=K+1
101   EFE(I)=FIELD(K)
      DO 111 I=1,M
      K=K+1
111   FE(I)=FIELD(K)
      DO 120 I=1,4
      IF((EFE(I).GT.-10.0).AND.(RFE(I).NE.EFE(I))) GOTO 100
120   CONTINUE
      RETURN
      END
C
C
C*************************************************************
C************* PEAK VALUES ELECTRON DENSITY ******************
C*************************************************************
C
C
      real function FOUT(XMODIP,XLATI,XLONGI,UT,FF0)
C CALCULATES CRITICAL FREQUENCY FOF2/MHZ USING SUBROUTINE GAMMA1.
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME
C (DEC. HOURS), FF0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
      DIMENSION FF0(988)
      INTEGER QFC(9)
      DATA QFC/11,11,8,4,1,0,0,0,0/
      FOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,6,QFC,9,76,13,988,FF0)
      RETURN
      END
C
C
      real function XMOUT(XMODIP,XLATI,XLONGI,UT,XM0)
C CALCULATES PROPAGATION FACTOR M3000 USING THE SUBROUTINE GAMMA1.
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME
C (DEC. HOURS), XM0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
      DIMENSION XM0(441)
      INTEGER QM(7)
      DATA QM/6,7,5,2,1,0,0/
      XMOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,4,QM,7,49,9,441,XM0)
      RETURN
      END
C
C
C
      REAL FUNCTION HMF2EDS(XR,X,XM3,HOUR,DAY,ABML)
C     D. BILITZA, F2 layer peak height (HMF2) for corrected magnetic 
C     latitude (HMLAT) and sunspot number (XR) using coefficient
C     M3000(XM3) and ratio of FOF2/FOE.
C     [Ref. D. BILITZA ET. AL. TELECOMM.J.,46,549-553,1979]
C     hmF2EDS = SMI option includes HOUR (LT), DAY, ABS(MLAT) input and
C     parameters of COMMON Blocks: Solar zenith angle XHI and 
C     corrected magnetic latitute HMLAT
C
      COMMON /BLO11/XHI/B6/HMLAT,H0NE,GLT,FG1,FG2
      F1=(2.32E-3)*XR+0.222
      F2=1.2-(1.16E-2)*EXP((2.39E-2)*XR)
      F3=0.096*(XR-25.0)/150.0
      XMIN=1.7
      IF (ABS(HMLAT).LT.64.) GO TO 1
      XZ=(XHI-65.)/6.
      XZ1=XZ*XZ
      XZ1=XZ1*XZ1
      RBA=(XZ1+16.)**0.25
      RB=0.25*XZ/RBA+1.95
      XZ=(XHI-74.)/6.
      XZ1=XZ*XZ
      XZ1=XZ1*XZ1
      RHA=(XZ1+16.)**0.25
      RH=0.2*XZ/RHA+1.70
      XMIN=(RB-RH)*(XR-10.)/140.+RH
    1 IF (X.LT.XMIN) X=XMIN
      DELM=F1*(1.0-XR/150.0*EXP(-HMLAT*HMLAT/1600.0))/(X-F2)+F3
      HMF2EDS=1490./(XM3+DELM)-176.
      IF (ABML.LT.30.) HMF2EDS=HMF2EDS*DELH(HOUR,DAY,ABML)
      RETURN
      END
C
      FUNCTION DELH(ALST,TD,AFI)
      TD1=TD-14.
      IF(TD1.LE.0.) TD1=TD1+365.
      YD=(TD1+29.5)/30.5
      Y=YD
      IF(Y.GT.6.5)Y=13-Y
      A1=0.02+ZFUN(1.,-1.,(Y-2.9)/1.1)/8.+0.18*EFUN(Y,2.,0.6,2)
      A2=ZFUN(0.07,0.18,(YD-1.)/1.2)-ZFUN(1.,1.,(YD-13.)/0.6)*0.18
      A2=A2-0.08*EFUN(YD,9.,4.-YD/6.,2)
      A3=ZFUN(0.09,0.28,YD-1)-0.28*ZFUN(1.,1.,YD-13)
      A3=A3-0.26*EFUN(YD,7.7,2.7,2)-0.05*EFUN(YD,3.,0.5,2)
      A4=ZFUN(.14,-.12,(YD-3.7)/0.6)+ZFUN(1.,1.,(YD-10.1)/.4)*0.12
      A4=A4-ZFUN(1.,-1.,(Y-1.3)*5)*0.11+EFUN(YD,6.8,1.,2)/10.
      DELH=0.98+A1*EFUN(ALST,0.,1.8,2)+A2*EFUN(ALST,5.2,1.4,2)
      DELH=DELH-(0.6*A3+0.4*A1-A4)*EFUN(ALST,21.,1.3,2)
      DELH=DELH+(5.*A3+(A1-A3)*(ALST-19.))*
     1     ZFUN(1.,1.,(ALST-17.7)/0.6)/10.
      DELH=DELH+(1-DELH)*AFI**2*(90.-2.*AFI)/27000.
      RETURN
      END
C
      REAL FUNCTION EFUN(A,B,C,N)
      D=(A-B)/C
      D2=D*D
      IF (N.GT.2) D2=D2*D2
      EFUN=EXP(-(80.*D2)/(80.+D2))
      RETURN
      END
C
      REAL FUNCTION ZFUN(A,B,C)
      CC=C*C
      C4=CC*CC
      CC=C4*C4
      CC=(CC+256.)**0.125
      ZFUN=A+B*C/CC
      RETURN
      END

C
C
      REAL FUNCTION FOF1ED(YLATI,R,CHI)
c--------------------------------------------------------------
C CALCULATES THE F1 PEAK PLASMA FREQUENCY (FOF1/MHZ)
C FOR   DIP-LATITUDE (YLATI/DEGREE)
c       SMOOTHED ZURICH SUNSPOT NUMBER (R)
c       SOLAR ZENITH ANGLE (CHI/DEGREE)
C REFERENCE:
c       E.D.DUCHARME ET AL., RADIO SCIENCE 6, 369-378, 1971
C                                      AND 8, 837-839, 1973
c       HOWEVER WITH MAGNETIC DIP LATITUDE INSTEAD OF GEOMAGNETIC
c       DIPOLE LATITUDE, EYFRIG, 1979
C--------------------------------------------- D. BILITZA, 1988.
      COMMON/CONST/UMR,PI
      FOF1 = 0.0
      DLA =  YLATI
      CHI0 = 49.84733 + 0.349504 * DLA
      CHI100 = 38.96113 + 0.509932 * DLA
      CHIM = ( CHI0 + ( CHI100 - CHI0 ) * R / 100. )
      IF(CHI.GT.CHIM) GOTO 1
      F0 = 4.35 + DLA * ( 0.0058 - 1.2E-4 * DLA )
      F100 = 5.348 + DLA * ( 0.011 - 2.3E-4 * DLA )
      FS = F0 + ( F100 - F0 ) * R / 100.0
      XMUE = 0.093 + DLA * ( 0.0046 - 5.4E-5 * DLA ) + 3.0E-4 * R
      FOF1 = FS * COS( CHI * UMR ) ** XMUE
1       FOF1ED = FOF1
      RETURN
      END
C
C
      REAL FUNCTION FOEEDI(COV,XHI,XHIM,XLATI)
C-------------------------------------------------------
C CALCULATES FOE/MHZ BY THE EDINBURGH-METHOD.
C INPUT: MEAN 10.7CM SOLAR RADIO FLUX (COV), GEOGRAPHIC
C LATITUDE (XLATI/DEG), SOLAR ZENITH ANGLE (XHI/DEG AND
C XHIM/DEG AT NOON).
C REFERENCE:
C       KOURIS-MUGGELETON, CCIR DOC. 6/3/07, 1973
C       TROST, J. GEOPHYS. RES. 84, 2736, 1979 (was used
C               to improve the nighttime varition)
C D.BILITZA--------------------------------- AUGUST 1986.
      COMMON/CONST/UMR,PI
C variation with solar activity (factor A) ...............
      A=1.0+0.0094*(COV-66.0)
C variation with noon solar zenith angle (B) and with latitude (C)
      SL=COS(XLATI*UMR)
      IF(XLATI.LT.32.0) THEN
      SM=-1.93+1.92*SL
      C=23.0+116.0*SL
      ELSE
      SM=0.11-0.49*SL
      C=92.0+35.0*SL
      ENDIF
      if(XHIM.ge.90.) XHIM=89.999
      B = COS(XHIM*UMR) ** SM
C variation with solar zenith angle (D) ..........................
      IF(XLATI.GT.12.0) THEN
      SP=1.2
      ELSE
      SP=1.31
      ENDIF
C adjusted solar zenith angle during nighttime (XHIC) .............
      XHIC=XHI-3.*ALOG(1.+EXP((XHI-89.98)/3.))
      D=COS(XHIC*UMR)**SP
C determine foE**4 ................................................
      R4FOE=A*B*C*D
C minimum allowable foE (sqrt[SMIN])...............................
      SMIN=0.121+0.0015*(COV-60.)
      SMIN=SMIN*SMIN
      IF(R4FOE.LT.SMIN) R4FOE=SMIN
      FOEEDI=R4FOE**0.25
      RETURN
      END
C
C
      REAL FUNCTION XMDED(XHI,R,YW)
C D. BILITZA, 1978, CALCULATES ELECTRON DENSITY OF D MAXIMUM.
C XHI/DEG. IS SOLAR ZENITH ANGLE, R SMOOTHED ZURICH SUNSPOT NUMBER
C AND YW/M-3 THE ASSUMED CONSTANT NIGHT VALUE.
C [REF.: D.BILITZA, WORLD DATA CENTER A REPORT UAG-82,7,
C       BOULDER,1981]
      COMMON/CONST/UMR,PI
      Y=6.05E8+0.088E8*R
      Z=(-0.1/(ALOG(YW/Y)))**0.3704
      if(abs(z).gt.1.) z=sign(1.,z)
      SUXHI=ACOS(Z)
      IF (SUXHI.LT.1.0472) SUXHI=1.0472
      XXHI=XHI*UMR
      IF (XXHI.GT.SUXHI) GOTO 100
      X=COS(XXHI)
      XMDED=Y*EXP(-0.1/X**2.7)
      RETURN
100   XMDED=YW
      RETURN
      END
C
C
      REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,IHARM,NQ,
     &                          K1,M,MM,M3,SFE)
C CALCULATES GAMMA1=FOF2 OR M3000 USING CCIR NUMERICAL MAP
C COEFFICIENTS SFE(M3) FOR MODIFIED DIP LATITUDE (SMODIP/DEG)
C GEOGRAPHIC LATITUDE (SLAT/DEG) AND LONGITUDE (SLONG/DEG)
C AND UNIVERSIAL TIME (HOUR/DECIMAL HOURS).
C NQ(K1) IS AN INTEGER ARRAY GIVING THE HIGHEST DEGREES IN
C LATITUDE FOR EACH LONGITUDE HARMONIC.
C M=1+NQ1+2(NQ2+1)+2(NQ3+1)+... .
C SHEIKH,4.3.77.
      REAL*8 C(12),S(12),COEF(100),SUM
      DIMENSION NQ(K1),XSINX(13),SFE(M3)
      COMMON/CONST/UMR,PI
      HOU=(15.0*HOUR-180.0)*UMR
      S(1)=SIN(HOU)
      C(1)=COS(HOU)
      DO 250 I=2,IHARM
      C(I)=C(1)*C(I-1)-S(1)*S(I-1)
      S(I)=C(1)*S(I-1)+S(1)*C(I-1)
250   CONTINUE
      DO 300 I=1,M
      MI=(I-1)*MM
      COEF(I)=SFE(MI+1)
      DO 300 J=1,IHARM
      COEF(I)=COEF(I)+SFE(MI+2*J)*S(J)+SFE(MI+2*J+1)*C(J)
300   CONTINUE
      SUM=COEF(1)
      SS=SIN(SMODIP*UMR)
      S3=SS
      XSINX(1)=1.0
      INDEX=NQ(1)
      DO 350 J=1,INDEX
      SUM=SUM+COEF(1+J)*SS
      XSINX(J+1)=SS
      SS=SS*S3
350   CONTINUE
      XSINX(NQ(1)+2)=SS
      NP=NQ(1)+1
      SS=COS(SLAT*UMR)
      S3=SS
      DO 400 J=2,K1
      S0=SLONG*(J-1.)*UMR
      S1=COS(S0)
      S2=SIN(S0)
      INDEX=NQ(J)+1
      DO 450 L=1,INDEX
      NP=NP+1
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S1
      NP=NP+1
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S2
450   CONTINUE
      SS=SS*S3
400   CONTINUE
      GAMMA1=SUM
      RETURN
      END
C
C
C************************************************************
C*************** EARTH MAGNETIC FIELD ***********************
C**************************************************************
C
C
      SUBROUTINE GGM(ART,LONG,LATI,MLONG,MLAT)
C CALCULATES GEOMAGNETIC LONGITUDE (MLONG) AND LATITUDE (MLAT)
C FROM GEOGRAFIC LONGITUDE (LONG) AND LATITUDE (LATI) FOR ART=0
C AND REVERSE FOR ART=1. ALL ANGLES IN DEGREE.
C LATITUDE:-90 TO 90. LONGITUDE:0 TO 360 EAST.
      INTEGER ART,DAYNR
      REAL MLONG,MLAT,LONG,LATI,longi
      COMMON/CONST/FAKTOR,PI
     &/dem/w,glongi,HOUR,daynr,cn1000,alon,blon,enre,HSC
      ZPI=FAKTOR*360.
      CBG=11.4*FAKTOR
      CI=COS(CBG)
      SI=SIN(CBG)
      IF(ART.EQ.0) GOTO 10
      CBM=COS(MLAT*FAKTOR)
      SBM=SIN(MLAT*FAKTOR)
      CLM=COS(MLONG*FAKTOR)
      SLM=SIN(MLONG*FAKTOR)
      SBG=SBM*CI-CBM*CLM*SI
      IF(ABS(SBG).GT.1.) SBG=SIGN(1.,SBG)
      LATI=ASIN(SBG)
      CBG=COS(LATI)
      SLG=(CBM*SLM)/CBG
      CLG=(SBM*SI+CBM*CLM*CI)/CBG
      IF(ABS(CLG).GT.1.) CLG=SIGN(1.,CLG)
      LONG=ACOS(CLG)
      IF(SLG.LT.0.0) LONG=ZPI-LONG
      LATI=LATI/FAKTOR
      LONG=LONG/FAKTOR
      LONG=LONG-69.8
      IF(LONG.LT.0.0) LONG=LONG+360.0
      longi=long
      RETURN
10    YLG=LONG+69.8
      CBG=COS(LATI*FAKTOR)
      SBG=SIN(LATI*FAKTOR)
      CLG=COS(YLG*FAKTOR)
      SLG=SIN(YLG*FAKTOR)
      SBM=SBG*CI+CBG*CLG*SI
      IF(ABS(SBM).GT.1.) SBM=SIGN(1.,SBM)
      MLAT=ASIN(SBM)
      CBM=COS(MLAT)
      SLM=(CBG*SLG)/CBM
      CLM=(-SBG*SI+CBG*CLG*CI)/CBM
      IF(ABS(CLM).GT.1.) CLM=SIGN(1.,CLM)
      MLONG=ACOS(CLM)
      IF(SLM.LT..0) MLONG=ZPI-MLONG
      MLAT=MLAT/FAKTOR
      MLONG=MLONG/FAKTOR
      RETURN
      END
C
C
      SUBROUTINE FIELDG(DLAT,DLONG,ALT,X,Y,Z,F,DIP,DEC,SMODIP)
C THIS IS A SPECIAL VERSION OF THE POGO 68/10 MAGNETIC FIELD
C LEGENDRE MODEL. TRANSFORMATION COEFF. G(144) VALID FOR 1973.
C INPUT: DLAT, DLONG=GEOGRAPHIC COORDINATES/DEG.(-90/90,0/360),
C        ALT=ALTITUDE/KM.
C OUTPUT: F TOTAL FIELD (GAUSS), Z DOWNWARD VERTICAL COMPONENT
C        X,Y COMPONENTS IN THE EQUATORIAL PLANE (X TO ZERO LONGITUDE).
C        DIP INCLINATION ANGLE(DEGREE). SMODIP RAWER'S MODFIED DIP.
C SHEIK,1977.
      DIMENSION H(144),XI(3),G(144),FEL1(72),FEL2(72)
      COMMON/CONST/UMR,PI
      DATA FEL1/0.0, 0.1506723,0.0101742, -0.0286519, 0.0092606,
     & -0.0130846, 0.0089594, -0.0136808,-0.0001508, -0.0093977,
     & 0.0130650, 0.0020520, -0.0121956, -0.0023451, -0.0208555,
     & 0.0068416,-0.0142659, -0.0093322, -0.0021364, -0.0078910,
     & 0.0045586,  0.0128904, -0.0002951, -0.0237245,0.0289493,
     & 0.0074605, -0.0105741, -0.0005116, -0.0105732, -0.0058542,
     &0.0033268, 0.0078164,0.0211234, 0.0099309, 0.0362792,
     &-0.0201070,-0.0046350,-0.0058722,0.0011147,-0.0013949,
     & -0.0108838,  0.0322263, -0.0147390,  0.0031247, 0.0111986,
     & -0.0109394,0.0058112,  0.2739046, -0.0155682, -0.0253272,
     &  0.0163782, 0.0205730,  0.0022081, 0.0112749,-0.0098427,
     & 0.0072705, 0.0195189, -0.0081132, -0.0071889, -0.0579970,
     & -0.0856642, 0.1884260,-0.7391512, 0.1210288, -0.0241888,
     & -0.0052464, -0.0096312, -0.0044834, 0.0201764,  0.0258343,
     &0.0083033,  0.0077187/
      DATA FEL2/0.0586055,0.0102236,-0.0396107,
     & -0.0167860, -0.2019911, -0.5810815,0.0379916,  3.7508268,
     & 1.8133030, -0.0564250, -0.0557352, 0.1335347, -0.0142641,
     & -0.1024618,0.0970994, -0.0751830,-0.1274948, 0.0402073,
     &  0.0386290, 0.1883088,  0.1838960, -0.7848989,0.7591817,
     & -0.9302389,-0.8560960, 0.6633250, -4.6363869, -13.2599277,
     & 0.1002136,  0.0855714,-0.0991981, -0.0765378,-0.0455264,
     &  0.1169326, -0.2604067, 0.1800076, -0.2223685, -0.6347679,
     &0.5334222, -0.3459502,-0.1573697,  0.8589464, 1.7815990,
     &-6.3347645, -3.1513653, -9.9927750,13.3327637, -35.4897308,
     &37.3466339, -0.5257398,  0.0571474, -0.5421217,  0.2404770,
     & -0.1747774,-0.3433644, 0.4829708,0.3935944, 0.4885033,
     &  0.8488121, -0.7640999, -1.8884945, 3.2930784,-7.3497229,
     & 0.1672821,-0.2306652, 10.5782146, 12.6031065, 8.6579742,
     & 215.5209961, -27.1419220,22.3405762,1108.6394043/
      K=0
      DO 10 I=1,72
      K=K+1
      G(K)=FEL1(I)
10    G(72+K)=FEL2(I)
      RLAT=DLAT*UMR
      CT=SIN(RLAT)
      ST=COS(RLAT)
      NMAX=11
      D=SQRT(40680925.0-272336.0*CT*CT)
      RLONG=DLONG*UMR
      CP=COS(RLONG)
      SP=SIN(RLONG)
      ZZZ=(ALT+40408589.0/D)*CT/6371.2
      RHO=(ALT+40680925.0/D)*ST/6371.2
      XXX=RHO*CP
      YYY=RHO*SP
      RQ=1.0/(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      XI(1)=XXX*RQ
      XI(2)=YYY*RQ
      XI(3)=ZZZ*RQ
      IHMAX=NMAX*NMAX+1
      LAST=IHMAX+NMAX+NMAX
      IMAX=NMAX+NMAX-1
      DO 100 I=IHMAX,LAST
100   H(I)=G(I)
      DO 200 K=1,3,2
      I=IMAX
      IH=IHMAX
300   IL=IH-I
      F1=2./(I-K+2.)
      X1=XI(1)*F1
      Y1=XI(2)*F1
      Z1=XI(3)*(F1+F1)
      I=I-2
      IF((I-1).LT.0) GOTO 400
      IF((I-1).EQ.0) GOTO 500
      DO 600 M=3,I,2
      H(IL+M+1)=G(IL+M+1)+Z1*H(IH+M+1)+X1*(H(IH+M+3)-H(IH+M-1))-
     &Y1*(H(IH+M+2)+H(IH+M-2))
      H(IL+M)=G(IL+M)+Z1*H(IH+M)+X1*(H(IH+M+2)-H(IH+M-2))+
     &Y1*(H(IH+M+3)+H(IH+M-1))
600   CONTINUE
500   H(IL+2)=G(IL+2)+Z1*H(IH+2)+X1*H(IH+4)-Y1*(H(IH+3)+H(IH))
      H(IL+1)=G(IL+1)+Z1*H(IH+1)+Y1*H(IH+4)+X1*(H(IH+3)-H(IH))
400   H(IL)=G(IL)+Z1*H(IH)+2.0*(X1*H(IH+1)+Y1*H(IH+2))
700   IH=IL
      IF(I.GE.K) GOTO 300
200   CONTINUE
      S=0.5*H(1)+2.0*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      XT=(RQ+RQ)*SQRT(RQ)
      X=XT*(H(3)-S*XXX)
      Y=XT*(H(4)-S*YYY)
      Z=XT*(H(2)-S*ZZZ)
      F=SQRT(X*X+Y*Y+Z*Z)
      BRH0=Y*SP+X*CP
      Y=Y*CP-X*SP
      X=Z*ST-BRH0*CT
      Z=-Z*CT-BRH0*ST
      zdivf=z/f
      IF(ABS(zdivf).GT.1.) zdivf=SIGN(1.,zdivf)
      DIP=ASIN(zdivf)
      ydivs=y/sqrt(x*x+y*y)
      IF(ABS(ydivs).GT.1.) ydivs=SIGN(1.,ydivs)
      DEC=ASIN(ydivs)
      dipdiv=DIP/SQRT(DIP*DIP+ST)
      IF(ABS(dipdiv).GT.1.) dipdiv=SIGN(1.,dipdiv)
      SMODIP=ASIN(dipdiv)
      DIP=DIP/UMR
      DEC=DEC/UMR
      SMODIP=SMODIP/UMR
      RETURN
      END
C
C
C************************************************************
C*********** INTERPOLATION AND REST ***************************
C**************************************************************
C
C
      SUBROUTINE REGFA1(X11,X22,FX11,FX22,EPS,FW,F,SCHALT,X)
C REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE
C STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL
C HAS BECOME LESS THAN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)
C THEN SCHALT=.TRUE.
      LOGICAL L1,LINKS,K,SCHALT
      SCHALT=.FALSE.
      EP=EPS
      X1=X11
      X2=X22
      F1=FX11-FW
      F2=FX22-FW
      K=.FALSE.
      NG=2
      LFD=0
      IF(F1*F2.LE.0.0) GOTO 200
      X=0.0
      SCHALT=.TRUE.
      RETURN
200   X=(X1*F2-X2*F1)/(F2-F1)
      GOTO 400
300     L1=LINKS
      DX=(X2-X1)/NG
      IF(.NOT.LINKS) DX=DX*(NG-1)
      X=X1+DX
400   FX=F(X)-FW
      LFD=LFD+1
      IF(LFD.GT.20) THEN
      EP=EP*10.
      LFD=0
      ENDIF
      LINKS=(F1*FX.GT.0.0)
      K=.NOT.K
      IF(LINKS) THEN
      X1=X
      F1=FX
      ELSE
      X2=X
      F2=FX
      ENDIF
      IF(ABS(X2-X1).LE.EP) GOTO 800
      IF(K) GOTO 300
      IF((LINKS.AND.(.NOT.L1)).OR.(.NOT.LINKS.AND.L1)) NG=2*NG
      GOTO 200
800   RETURN
      END
C
C
      SUBROUTINE TAL(SHABR,SDELTA,SHBR,SDTDH0,AUS6,SPT)
C CALCULATES THE COEFFICIENTS SPT FOR THE POLYNOMIAL
C Y(X)=1+SPT(1)*X**2+SPT(2)*X**3+SPT(3)*X**4+SPT(4)*X**5
C TO FIT THE VALLEY IN Y, REPRESENTED BY:
C Y(X=0)=1, THE X VALUE OF THE DEEPEST VALLEY POINT (SHABR),
C THE PRECENTAGE DEPTH (SDELTA), THE WIDTH (SHBR) AND THE
C DERIVATIVE DY/DX AT THE UPPER VALLEY BOUNDRY (SDTDH0).
C IF THERE IS AN UNWANTED ADDITIONAL EXTREMUM IN THE VALLEY
C REGION, THEN AUS6=.TRUE., ELSE AUS6=.FALSE..
C FOR -SDELTA THE COEFF. ARE CALCULATED FOR THE FUNCTION
C Y(X)=EXP(SPT(1)*X**2+...+SPT(4)*X**5).
      DIMENSION SPT(4)
      LOGICAL AUS6
      Z1=-SDELTA/(100.0*SHABR*SHABR)
      IF(SDELTA.GT.0.) GOTO 500
      SDELTA=-SDELTA
      Z1=ALOG(1.-SDELTA/100.)/(SHABR*SHABR)
500   Z3=SDTDH0/(2.*SHBR)
      Z4=SHABR-SHBR
      SPT(4)=2.0*(Z1*(SHBR-2.0*SHABR)*SHBR+Z3*Z4*SHABR)/
     &  (SHABR*SHBR*Z4*Z4*Z4)
      SPT(3)=Z1*(2.0*SHBR-3.0*SHABR)/(SHABR*Z4*Z4)-
     &  (2.*SHABR+SHBR)*SPT(4)
      SPT(2)=-2.0*Z1/SHABR-2.0*SHABR*SPT(3)-3.0*SHABR*SHABR*SPT(4)
      SPT(1)=Z1-SHABR*(SPT(2)+SHABR*(SPT(3)+SHABR*SPT(4)))
      AUS6=.FALSE.
      B=4.*SPT(3)/(5.*SPT(4))+SHABR
      C=-2.*SPT(1)/(5*SPT(4)*SHABR)
      Z2=B*B/4.-C
      IF(Z2.LT.0.0) GOTO 300
      Z3=SQRT(Z2)
      Z1=B/2.
      Z2=-Z1+Z3
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.
      IF (ABS(Z3).GT.1.E-15) GOTO 400
      Z2=C/Z2
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.
      RETURN
400   Z2=-Z1-Z3
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.
300   RETURN
      END
C
C

C******************************************************************
C********** ZENITH ANGLE, DAY OF YEAR, TIME ***********************
C******************************************************************
C
C
      subroutine soco (ld,t,flat,Elon,
     &          DECLIN, ZENITH, SUNRSE, SUNSET)
c--------------------------------------------------------------------
c       s/r to calculate the solar declination, zenith angle, and
c       sunrise & sunset times  - based on Newbern Smith's algorithm
c       [leo mcnamara, 1-sep-86, last modified 16-jun-87]
c       {dieter bilitza, 30-oct-89, modified for IRI application}
c
c in:   ld      local day of year
c       t       local hour (decimal)
c       flat    northern latitude in degrees
c       elon    east longitude in degrees
c
c out:  declin      declination of the sun in degrees
c       zenith      zenith angle of the sun in degrees
c       sunrse      local time of sunrise in hours
c       sunset      local time of sunset in hours
c-------------------------------------------------------------------
c
      common/const/   dtr,PI
c amplitudes of Fourier coefficients  --  1955 epoch.................
      data    p1,p2,p3,p4,p6 /
     &  0.017203534,0.034407068,0.051610602,0.068814136,0.103221204 /
c
c s/r is formulated in terms of WEST longitude.......................
      wlon = 360. - Elon
c
c time of equinox for 1980...........................................
      td = ld + (t + Wlon/15.) / 24.
      te = td + 0.9369
c
c declination of the sun..............................................
      dcl = 23.256 * sin(p1*(te-82.242)) + 0.381 * sin(p2*(te-44.855))
     &      + 0.167 * sin(p3*(te-23.355)) - 0.013 * sin(p4*(te+11.97))
     &      + 0.011 * sin(p6*(te-10.41)) + 0.339137
      DECLIN = dcl
      dc = dcl * dtr
c
c the equation of time................................................
      tf = te - 0.5
      eqt = -7.38*sin(p1*(tf-4.)) - 9.87*sin(p2*(tf+9.))
     &      + 0.27*sin(p3*(tf-53.)) - 0.2*cos(p4*(tf-17.))
      et = eqt * dtr / 4.
c
      fa = flat * dtr
      phi = 0.26179939 * ( t - 12.) + et
c
      a = sin(fa) * sin(dc)
      b = cos(fa) * cos(dc)
      cosx = a + b * cos(phi)
      if(abs(cosx).gt.1.) cosx=sign(1.,cosx)
      zenith = acos(cosx) / dtr
c
c calculate sunrise and sunset times --  at the ground...........
c see Explanatory Supplement to the Ephemeris (1961) pg 401......
c sunrise at height h metres is at...............................
c       chi(h) = 90.83 + 0.0347 * sqrt(h)........................
c this includes corrections for horizontal refraction and........
c semi-diameter of the solar disk................................
      ch = cos(90.83 * dtr)
      cosphi = (ch -a ) / b
c if abs(secphi) > 1., sun does not rise/set.....................
c allow for sun never setting - high latitude summer.............
      secphi = 999999.
      if(cosphi.ne.0.) secphi = 1./cosphi
      sunset = 99.
      sunrse = 99.
      if(secphi.gt.-1.0.and.secphi.le.0.) return
c allow for sun never rising - high latitude winter..............
      sunset = -99.
      sunrse = -99.
      if(secphi.gt.0.0.and.secphi.lt.1.) return
c
      if(cosphi.gt.1.) cosphi=sign(1.,cosphi)
      phi = acos(cosphi)
      et = et / 0.26179939
      phi = phi / 0.26179939
      sunrse = 12. - phi - et
      sunset = 12. + phi - et
      if(sunrse.lt.0.) sunrse = sunrse + 24.
      if(sunset.ge.24.) sunset = sunset - 24.
c
      return
      end
c
C
      FUNCTION HPOL(HOUR,TW,XNW,SA,SU,DSA,DSU)
C-------------------------------------------------------
C PROCEDURE FOR SMOOTH TIME-INTERPOLATION USING EPSTEIN
C STEP FUNCTION AT SUNRISE (SA) AND SUNSET (SU). THE
C STEP-WIDTH FOR SUNRISE IS DSA AND FOR SUNSET DSU.
C TW,NW ARE THE DAY AND NIGHT VALUE OF THE PARAMETER TO
C BE INTERPOLATED. SA AND SU ARE TIME OF SUNRIES AND
C SUNSET IN DECIMAL HOURS.
C BILITZA----------------------------------------- 1979.
      IF(ABS(SU).GT.25.) THEN
      IF(SU.GT.0.0) THEN
         HPOL=TW
      ELSE
         HPOL=XNW
      ENDIF
      RETURN
      ENDIF
      HPOL=XNW+(TW-XNW)*EPST(HOUR,DSA,SA)+
     &  (XNW-TW)*EPST(HOUR,DSU,SU)
      RETURN
      END
C
C
      SUBROUTINE MODA(IN,IYEAR,MONTH,IDAY,IDOY,NRDAYMO)
C-------------------------------------------------------------------
C CALCULATES DAY OF YEAR (IDOY, ddd) FROM YEAR (IYEAR, yy or yyyy),
C MONTH (MONTH, mm) AND DAY OF MONTH (IDAY, dd) IF IN=0, OR MONTH
C AND DAY FROM YEAR AND DAY OF YEAR IF IN=1. NRDAYMO is an output
C parameter providing the number of days in the specific month.
C-------------------------------------------------------------------
      DIMENSION       MM(12)
      DATA            MM/31,28,31,30,31,30,31,31,30,31,30,31/

      IMO=0
      MOBE=0

      if((iyear/4*4.eq.iyear).and.(iyear/100*100.ne.iyear)) mm(2)=29

      IF(IN.GT.0) GOTO 5
      mosum=0
      if(month.gt.1) then
         do 1234 i=1,month-1
1234                            mosum=mosum+mm(i)
         endif
      idoy=mosum+iday
      nrdaymo=mm(month)
      RETURN

5       IMO=IMO+1
      IF(IMO.GT.12) GOTO 55
      MOOLD=MOBE
      nrdaymo=mm(imo)
      MOBE=MOBE+nrdaymo
      IF(MOBE.LT.IDOY) GOTO 5
55              MONTH=IMO
      IDAY=IDOY-MOOLD
      RETURN
      END
c
c
      subroutine ut_lt(mode,ut,slt,glong,iyyy,ddd)
c -----------------------------------------------------------------
c Converts Universal Time UT (decimal hours) into Solar Local Time
c SLT (decimal hours) for given date (iyyy is year, e.g. 1995; ddd
c is day of year, e.g. 1 for Jan 1) and geodatic longitude in degrees.
C For mode=0 UT->LT and for mode=1 LT->UT
c Please NOTE that iyyy and ddd are input as well as output parameters
c since the determined LT may be for a day before or after the UT day.
c ------------------------------------------------- bilitza nov 95
      integer     ddd,dddend

        xlong=glong
        if(glong.gt.180) xlong=glong-360
        if(mode.ne.0) goto 1
c
c UT ---> LT
c
        SLT=UT+xlong/15.
        if((SLT.ge.0.).and.(SLT.le.24.)) goto 2
      if(SLT.gt.24.) goto 3
      SLT=SLT+24.
                ddd=ddd-1
                if(ddd.lt.1.) then
                  iyyy=iyyy-1
                        ddd=365
      if((iyyy/4*4.eq.iyyy).and.(iyyy/100*100.ne.iyyy)) ddd=366
                        endif
      goto 2
3     SLT=SLT-24.
      ddd=ddd+1
      dddend=365
      if((iyyy/4*4.eq.iyyy).and.(iyyy/100*100.ne.iyyy)) dddend=366
      if(ddd.gt.dddend) then
         iyyy=iyyy+1
         ddd=1
         endif
      goto 2
c
c LT ---> UT
c
1     UT=SLT-xlong/15.
        if((UT.ge.0.).and.(UT.le.24.)) goto 2
      if(UT.gt.24.) goto 5
      UT=UT+24.
                ddd=ddd-1
                if(ddd.lt.1.) then
                  iyyy=iyyy-1
                        ddd=365
      if((iyyy/4*4.eq.iyyy).and.(iyyy/100*100.ne.iyyy)) ddd=366
                        endif
      goto 2
5     UT=UT-24.
      ddd=ddd+1
      dddend=365
      if((iyyy/4*4.eq.iyyy).and.(iyyy/100*100.ne.iyyy)) dddend=366
      if(ddd.gt.dddend) then
         iyyy=iyyy+1
         ddd=1
         endif
2     return
      end
c
C
C
C
        REAL FUNCTION B0_98 ( HOUR, SAX, SUX, NSEASN, R, ZLO, ZMODIP)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON is northern season with
C ISEASON=1 northern spring), low and high solar activity Rz12=10,
C 100 (IR=1,2), and modified dip latitudes of 0, 18 and 45
C degress (ILATI=1,2,3). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C
C input:
C       hour    LT in decimal hours
C       SAX     time of sunrise in decimal hours
C       SUX     time of sunset in decimal hours
C       nseasn  season in northern hemisphere (1=spring)
C       R       12-month running mean of sunspot number
C       ZLO     longitude
C       ZMODIP  modified dip latitude
C
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
C Updates (B0_new -> B0_98):
C
C 01/98 corrected to include a smooth transition at the modip equator
C       and no discontinuity at the equatorial change in season.
C 09/98 new B0 values incl values at the magnetic equator
C 10/98 longitude as input to determine if magnetic equator in northern 
C         or southern hemisphere
C
      REAL      NITVAL
      DIMENSION B0F(2,4,2,3),bfr(2,2,3),bfd(2,3),zx(5),g(6),dd(5)
      DATA      B0F/201,68,210,61,192,68,199,67,240,80,245,83,
     &              233,71,230,65,108,65,142,81,110,68,77,75,
     &              124,98,164,100,120,94,96,112,78,81,94,84,
     &              81,81,65,70,102,87,127,91,109,88,81,78/
        data    zx/45.,72.,90.,108.,135./,dd/5*3.0/

        num_lat=3

C jseasn is southern hemisphere season
        jseasn=nseasn+2
        if(jseasn.gt.4) jseasn=jseasn-4

        zz = zmodip + 90.
        zz0 = 0.

C Interpolation in Rz12: linear from 10 to 100
        DO 7035 ISL=1,num_lat
          DO 7034 ISD=1,2
            bfr(isd,1,isl) = b0f(isd,nseasn,1,isl) +
     &      (b0f(isd,nseasn,2,isl) - b0f(isd,nseasn,1,isl))/90.*(R-10.)
            bfr(isd,2,isl) = b0f(isd,jseasn,1,isl) +
     &      (b0f(isd,jseasn,2,isl) - b0f(isd,jseasn,1,isl))/90.*(R-10.)
7034      continue
C Interpolation day/night with transitions at SAX (sunrise)
C and SUX (sunset) for northern/southern hemisphere iss=1/2
          do 7033 iss=1,2
                DAYVAL = BFR(1,ISS,ISL)
                NITVAL = BFR(2,ISS,ISL)
                BFD(iss,ISL) = HPOL(HOUR,DAYVAL,NITVAL,SAX,SUX,1.,1.)
7033      continue
7035    continue

C Interpolation with epstein-transitions in modified dip latitude.
C Transitions at +/-18 and +/-45 degrees; constant above +/-45.
C
C g(1:5) are the latitudinal slopes of B0;
C       g(1) is for the region from -90 to -45 degrees
C       g(2) is for the region from -45 to -18 degrees
C       g(3) is for the region from -18 to   0 degrees
C       g(4) is for the region from   0 to  18 degrees
C       g(5) is for the region from  18 to  45 degrees
C       g(6) is for the region from  45 to  90 degrees
C
C B0 =  bfd(2,3) at modip = -45,
C       bfd(2,2) at modip = -18,
C       bfd(2,1) or bfd(1,1) at modip = 0,
C       bfd(1,2) at modip = 20,
C       bfd(1,3) at modip = 45.
C If the Longitude is between 200 and 320 degrees than the modip 
C equator is in the southern hemisphere and bfd(2,1) is used at the 
C equator, otherwise bfd(1,1) is used.
c
        zx1=bfd(2,3)
        zx2=bfd(2,2)
        zx3=bfd(1,1)
        if(zlo.gt.200.0.and.zlo.lt.320) zx3=bfd(2,1)
        zx4=bfd(1,2)
        zx5=bfd(1,3)
        g(1) = 0.
        g(2) = ( zx2 - zx1 ) / 27.
        g(3) = ( zx3 - zx2 ) / 18.
        g(4) = ( zx4 - zx3 ) / 18.
        g(5) = ( zx5 - zx4 ) / 27.
        g(6) = 0.

c        bb0 = bfd(2,3)
c      SUM = bb0
        sum=zx1
      DO 1 I=1,5
        aa = eptr(zz ,dd(i),zx(i))
        bb = eptr(zz0,dd(i),zx(i))
        DSUM = (G(I+1) - G(I)) * (AA-BB) * dd(i)
        SUM = SUM + DSUM
1       continue
      B0_98 = SUM

        RETURN
        END
c
C
      REAL FUNCTION B0POL ( HOUR, SAX, SUX, ISEASON, R, DELA)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON=1 spring), low and high
C solar activity (IR=1,2), and low and middle modified dip
C latitudes (ILATI=1,2). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
      REAL              NITVAL
      DIMENSION         B0F(2,4,2,2),SIPH(2),SIPL(2)
      DATA      B0F/114.,64.0,134.,77.0,128.,66.0,75.,73.0,
     &              113.,115.,150.,116.,138.,123.,94.,132.,
     &              72.0,84.0,83.0,89.0,75.0,85.0,57.,76.0,
     &              102.,100.,120.,110.,107.,103.,76.,86.0/

      DO 7033 ISR=1,2
      DO 7034 ISL=1,2
         DAYVAL   = B0F(1,ISEASON,ISR,ISL)
         NITVAL = B0F(2,ISEASON,ISR,ISL)

C Interpolation day/night with transitions at SAX (sunrise) and SUX (sunset)
7034                    SIPH(ISL) = HPOL(HOUR,DAYVAL,NITVAL,
     &                          SAX,SUX,1.,1.)

C Interpolation low/middle modip with transition at 30 degrees modip
7033            SIPL(ISR) = SIPH(1) + (SIPH(2) - SIPH(1)) / DELA

C Interpolation low/high Rz12: linear from 10 to 100
      B0POL=SIPL(1)+(SIPL(2)-SIPL(1))/90.*(R-10.)
      RETURN
      END
c
C
C *********************************************************************
C
      SUBROUTINE SUN (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C-----------------------------------------------------------------------------
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C  LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION DJ,FDAY
      DATA RAD/57.295779513/
C
      IF(IYEAR.LT.1901.OR.IYEAR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
      G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END
C      
C
C ************************ EPSTEIN FUNCTIONS **************************
C *********************************************************************
C REF:  H. G. BOOKER, J. ATMOS. TERR. PHYS. 39, 619-623, 1977
C       K. RAWER, ADV. SPACE RES. 4, #1, 11-15, 1984
C *********************************************************************
C
C
      REAL FUNCTION  RLAY ( X, XM, SC, HX )
C -------------------------------------------------------- RAWER  LAYER
      Y1  = EPTR ( X , SC, HX )
      Y1M = EPTR ( XM, SC, HX )
      Y2M = EPST ( XM, SC, HX )
      RLAY = Y1 - Y1M - ( X - XM ) * Y2M / SC
      RETURN
      END
C
C
      REAL FUNCTION D1LAY ( X, XM, SC, HX )
C ------------------------------------------------------------ dLAY/dX
      D1LAY = ( EPST(X,SC,HX) - EPST(XM,SC,HX) ) /  SC
      RETURN
      END
C
C
C#      REAL FUNCTION D2LAY ( X, XM, SC, HX )
         REAL FUNCTION D2LAY ( X,SC, HX )
C ---------------------------------------------------------- d2LAY/dX2
      D2LAY = EPLA(X,SC,HX) /  (SC * SC)
      RETURN
      END
C
C
      REAL FUNCTION EPTR ( X, SC, HX )
C ------------------------------------------------------------ TRANSITION
      COMMON/ARGEXP/ARGMAX
      D1 = ( X - HX ) / SC
      IF (ABS(D1).LT.ARGMAX) GOTO 1
      IF (D1.GT.0.0) THEN
        EPTR = D1
      ELSE
        EPTR = 0.0
      ENDIF
      RETURN
1       EPTR = ALOG ( 1. + EXP( D1 ))
      RETURN
      END
C
C
      REAL FUNCTION EPST ( X, SC, HX )
C -------------------------------------------------------------- STEP
      COMMON/ARGEXP/ARGMAX
      D1 = ( X - HX ) / SC
      IF (ABS(D1).LT.ARGMAX) GOTO 1
      IF (D1.GT.0.0) THEN
        EPST = 1.
      ELSE
        EPST = 0.
      ENDIF
      RETURN
1       EPST = 1. / ( 1. + EXP( -D1 ))
      RETURN
      END
C
C
      REAL FUNCTION EPSTEP ( Y2, Y1, SC, HX, X)
C---------------------------------------------- STEP FROM Y1 TO Y2
      EPSTEP = Y1 + ( Y2 - Y1 ) * EPST ( X, SC, HX)
      RETURN
      END
C
C
      REAL FUNCTION EPLA ( X, SC, HX )
C ------------------------------------------------------------ PEAK
      COMMON/ARGEXP/ARGMAX
      D1 = ( X - HX ) / SC
      IF (ABS(D1).LT.ARGMAX) GOTO 1
      EPLA = 0
      RETURN
1       D0 = EXP ( D1 )
      D2 = 1. + D0
      EPLA = D0 / ( D2 * D2 )
      RETURN
      END
c
c
      FUNCTION XE2TO5(H,HMF2,NL,HX,SC,AMP)
C----------------------------------------------------------------------
C NORMALIZED ELECTRON DENSITY (N/NMF2) FOR THE MIDDLE IONOSPHERE FROM
C HME TO HMF2 USING LAY-FUNCTIONS.
C----------------------------------------------------------------------
      DIMENSION       HX(NL),SC(NL),AMP(NL)
      SUM = 1.0
      DO 1 I=1,NL
      YLAY = AMP(I) * RLAY( H, HMF2, SC(I), HX(I) )
      zlay=10.**ylay
1          sum=sum*zlay
      XE2TO5 = sum
      RETURN
      END
C
C
      REAL FUNCTION XEN(H,HMF2,XNMF2,HME,NL,HX,SC,AMP)
C----------------------------------------------------------------------
C ELECTRON DENSITY WITH NEW MIDDLE IONOSPHERE
C----------------------------------------------------------------------
      DIMENSION       HX(NL),SC(NL),AMP(NL)
C
      IF(H.LT.HMF2) GOTO 100
      XEN = XE1(H)
      RETURN
100     IF(H.LT.HME) GOTO 200
      XEN = XNMF2 * XE2TO5(H,HMF2,NL,HX,SC,AMP)
      RETURN
200     XEN = XE6(H)
      RETURN
      END
C
C
      SUBROUTINE VALGUL(XHI,HVB,VWU,VWA,VDP)
C ---------------------------------------------------------------------
C   CALCULATES E-F VALLEY PARAMETERS; T.L. GULYAEVA, ADVANCES IN
C   SPACE RESEARCH 7, #6, 39-48, 1987.
C
C       INPUT:  XHI     SOLAR ZENITH ANGLE [DEGREE]
C
C       OUTPUT: VDP     VALLEY DEPTH  (NVB/NME)
C               VWU     VALLEY WIDTH  [KM]
C               VWA     VALLEY WIDTH  (SMALLER, CORRECTED BY RAWER)
C               HVB     HEIGHT OF VALLEY BASE [KM]
C -----------------------------------------------------------------------
C
      COMMON  /CONST/UMR,PI
C
      CS = 0.1 + COS(UMR*XHI)
      ABC = ABS(CS)
      VDP = 0.45 * CS / (0.1 + ABC ) + 0.55
      ARL = ( 0.1 + ABC + CS ) / ( 0.1 + ABC - CS)
      ZZZ = ALOG( ARL )
      VWU = 45. - 10. * ZZZ
      VWA = 45. -  5. * ZZZ
      HVB = 1000. / ( 7.024 + 0.224 * CS + 0.966 * ABC )
      RETURN
      END
C
C
      SUBROUTINE ROGUL(IDAY,XHI,SX,GRO)
C ---------------------------------------------------------------------
C   CALCULATES RATIO H0.5/HMF2 FOR HALF-DENSITY POINT (NE(H0.5)=0.5*NMF2)
C   T.L. GULYAEVA, ADVANCES IN SPACE RESEARCH 7, #6, 39-48, 1987.
C
C       INPUT:  IDAY    DAY OF YEAR
C               XHI     SOLAR ZENITH ANGLE [DEGREE]
C
C       OUTPUT: GRO     RATIO OF HALF DENSITY HEIGHT TO F PEAK HEIGHT
C               SX      SMOOTHLY VARYING SEASON PARAMTER (SX=1 FOR
C                       DAY=1; SX=3 FOR DAY=180; SX=2 FOR EQUINOX)
C -----------------------------------------------------------------------
C
      SX = 2. - COS ( IDAY * 0.017214206 )
      XS = ( XHI - 20. * SX) / 15.
      GRO = 0.8 - 0.2 / ( 1. + EXP(XS) )
c same as gro=0.6+0.2/(1+exp(-xs))
      RETURN
      END
C
C
      SUBROUTINE LNGLSN ( N, A, B, AUS)
C --------------------------------------------------------------------
C SOLVES QUADRATIC SYSTEM OF LINEAR EQUATIONS:
C
C       INPUT:  N       NUMBER OF EQUATIONS (= NUMBER OF UNKNOWNS)
C               A(N,N)  MATRIX (LEFT SIDE OF SYSTEM OF EQUATIONS)
C               B(N)    VECTOR (RIGHT SIDE OF SYSTEM)
C
C       OUTPUT: AUS     =.TRUE.   NO SOLUTION FOUND
C                       =.FALSE.  SOLUTION IS IN  A(N,J) FOR J=1,N
C --------------------------------------------------------------------
C
      DIMENSION       A(5,5), B(5), AZV(10)
      LOGICAL         AUS
C
      AUS = .FALSE.
      DO 1 K=1,N-1
      IMAX = K
      L    = K
      IZG  = 0
      AMAX  = ABS( A(K,K) )
110             L = L + 1
      IF (L.GT.N) GOTO 111
      HSP = ABS( A(L,K) )
      IF (HSP.LT.1.E-8) IZG = IZG + 1
      IF (HSP.LE.AMAX) GOTO 110
111             IF (ABS(AMAX).GE.1.E-10) GOTO 133
         AUS = .TRUE.
         RETURN
133             IF (IMAX.EQ.K) GOTO 112
      DO 2 L=K,N
         AZV(L+1)  = A(IMAX,L)
         A(IMAX,L) = A(K,L)
2                       A(K,L)    = AZV(L+1)
      AZV(1)  = B(IMAX)
      B(IMAX) = B(K)
      B(K)    = AZV(1)
112             IF (IZG.EQ.(N-K)) GOTO 1
      AMAX = 1. / A(K,K)
      AZV(1) = B(K) * AMAX
      DO 3 M=K+1,N
3                       AZV(M+1) = A(K,M) * AMAX
      DO 4 L=K+1,N
         AMAX = A(L,K)
         IF (ABS(AMAX).LT.1.E-8) GOTO 4
         A(L,K) = 0.0
         B(L) = B(L) - AZV(1) * AMAX
         DO 5 M=K+1,N
5                               A(L,M) = A(L,M) - AMAX * AZV(M+1)
4               CONTINUE
1       CONTINUE
      DO 6 K=N,1,-1
      AMAX = 0.0
      IF (K.LT.N) THEN
         DO 7 L=K+1,N
7                               AMAX = AMAX + A(K,L) * A(N,L)
         ENDIF
      IF (ABS(A(K,K)).LT.1.E-6) THEN
         A(N,K) = 0.0
      ELSE
         A(N,K) = ( B(K) - AMAX ) / A(K,K)
      ENDIF
6       CONTINUE
      RETURN
      END
C
C
      SUBROUTINE LSKNM ( N, M, M0, M1, HM, SC, HX, W, X, Y, VAR, SING)
C --------------------------------------------------------------------
C   DETERMINES LAY-FUNCTIONS AMPLITUDES FOR A NUMBER OF CONSTRAINTS:
C
C       INPUT:  N       NUMBER OF AMPLITUDES ( LAY-FUNCTIONS)
C               M       NUMBER OF CONSTRAINTS
C               M0      NUMBER OF POINT CONSTRAINTS
C               M1      NUMBER OF FIRST DERIVATIVE CONSTRAINTS
C               HM      F PEAK ALTITUDE  [KM]
C               SC(N)   SCALE PARAMETERS FOR LAY-FUNCTIONS  [KM]
C               HX(N)   HEIGHT PARAMETERS FOR LAY-FUNCTIONS  [KM]
C               W(M)    WEIGHT OF CONSTRAINTS
C               X(M)    ALTITUDES FOR CONSTRAINTS  [KM]
C               Y(M)    LOG(DENSITY/NMF2) FOR CONSTRAINTS
C
C       OUTPUT: VAR(M)  AMPLITUDES
C               SING    =.TRUE.   NO SOLUTION
C ------------------------------------------------------------------------
C
      LOGICAL         SING
      DIMENSION       VAR(N), HX(N), SC(N), W(M), X(M), Y(M),
     &                  BLI(5), ALI(5,5), XLI(5,10)
C
      M01=M0+M1
      DO 1 J=1,5
      BLI(J) = 0.
      DO 1 I=1,5
1                       ALI(J,I) = 0.
      DO 2 I=1,N
      DO 3 K=1,M0
3                       XLI(I,K) = RLAY( X(K), HM, SC(I), HX(I) )
      DO 4 K=M0+1,M01
4                       XLI(I,K) = D1LAY( X(K), HM, SC(I), HX(I) )
      DO 5 K=M01+1,M
C#5                       XLI(I,K) = D2LAY( X(K), HM, SC(I), HX(I) )
5                       XLI(I,K) = D2LAY( X(K), SC(I), HX(I) )
2       CONTINUE
      DO 7 J=1,N
      DO 6 K=1,M
         BLI(J) = BLI(J) + W(K) * Y(K) * XLI(J,K)
         DO 6 I=1,N
6                               ALI(J,I) = ALI(J,I) + W(K) * XLI(I,K)
     &                                  * XLI(J,K)
7       CONTINUE
      CALL LNGLSN( N, ALI, BLI, SING )
      IF (.NOT.SING) THEN
      DO 8 I=1,N
8                       VAR(I) = ALI(N,I)
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE INILAY(NIGHT,XNMF2,XNMF1,XNME,VNE,HMF2,HMF1,
     &                          HME,HV1,HV2,HHALF,HXL,SCL,AMP,IQUAL)
C-------------------------------------------------------------------
C CALCULATES AMPLITUDES FOR LAY FUNCTIONS
C D. BILITZA, DECEMBER 1988
C
C INPUT:        NIGHT   LOGICAL VARIABLE FOR DAY/NIGHT DISTINCTION
C               XNMF2   F2 PEAK ELECTRON DENSITY [M-3]
C               XNMF1   F1 PEAK ELECTRON DENSITY [M-3]
C               XNME    E  PEAK ELECTRON DENSITY [M-3]
C               VNE     ELECTRON DENSITY AT VALLEY BASE [M-3]
C               HMF2    F2 PEAK ALTITUDE [KM]
C               HMF1    F1 PEAK ALTITUDE [KM]
C               HME     E  PEAK ALTITUDE [KM]
C               HV1     ALTITUDE OF VALLEY TOP [KM]
C               HV2     ALTITUDE OF VALLEY BASE [KM]
C               HHALF   ALTITUDE OF HALF-F2-PEAK-DENSITY [KM]
C
C OUTPUT:       HXL(4)  HEIGHT PARAMETERS FOR LAY FUNCTIONS [KM]
C               SCL(4)  SCALE PARAMETERS FOR LAY FUNCTIONS [KM]
C               AMP(4)  AMPLITUDES FOR LAY FUNCTIONS
C               IQUAL   =0 ok, =1 ok using second choice for HXL(1)
C                       =2 NO SOLUTION
C---------------------------------------------------------------
      DIMENSION       XX(8),YY(8),WW(8),AMP(4),HXL(4),SCL(4)
      LOGICAL         SSIN,NIGHT
c
c constants --------------------------------------------------------
      NUMLAY=4
      NC1 = 2
      ALG102=ALOG10(2.)
c
c constraints: xx == height     yy == log(Ne/NmF2)    ww == weights
c -----------------------------------------------------------------

      ALOGF = ALOG10(XNMF2)
      ALOGEF = ALOG10(XNME) - ALOGF
      XHALF=XNMF2/2.
      XX(1) = HHALF
      XX(2) = HV1
      XX(3) = HV2
      XX(4) = HME
      XX(5) = HME - ( HV2 - HME )
      YY(1) = -ALG102
      YY(2) = ALOGEF
      YY(3) = ALOG10(VNE) - ALOGF
      YY(4) = ALOGEF
      YY(5) = YY(3)
      YY(7) = 0.0
      WW(2) = 1.
      WW(3) = 2.
      WW(4) = 5.
c
c geometric paramters for LAY -------------------------------------
c difference to earlier version:  HXL(3) = HV2 + SCL(3)
c
      SCL0 = 0.7 * ( 0.216 * ( HMF2 - HHALF ) + 56.8 )
      SCL(1) = 0.8 * SCL0
      SCL(2) = 10.
      SCL(3) = 9.
      SCL(4) = 6.
      HXL(3) = HV2
c
C DAY CONDITION--------------------------------------------------
c earlier tested:       HXL(2) = HMF1 + SCL(2)
c
       IF(NIGHT) GOTO 7711
      NUMCON = 8
      HXL(1) = 0.9 * HMF2
        HXL1T  = HHALF
      HXL(2) = HMF1
      HXL(4) = HME - SCL(4)
      XX(6) = HMF1
      XX(7) = HV2
      XX(8) = HME
      YY(8) = 0.0
      WW(5) = 1.
      WW(7) = 50.
      WW(8) = 500.
c without F-region ----------------------------------------------
      IF(XNMF1.GT.0) GOTO 100
         HXL(2)=(HMF2+HHALF)/2.
         YY(6) = 0.
         WW(6) = 0.
         WW(1) = 1.
         GOTO 7722
c with F-region --------------------------------------------
100   CONTINUE

                YY(6) = ALOG10(XNMF1) - ALOGF
      WW(6) = 3.
      IF((XNMF1-XHALF)*(HMF1-HHALF).LT.0.0) THEN
        WW(1)=0.5
      ELSE
        ZET = YY(1) - YY(6)
        WW(1) = EPST( ZET, 0.1, 0.15)
      ENDIF
      IF(HHALF.GT.HMF1) THEN
        HFFF=HMF1
        XFFF=XNMF1
      ELSE
        HFFF=HHALF
        XFFF=XHALF
      ENDIF
      GOTO 7722
c
C NIGHT CONDITION---------------------------------------------------
c different HXL,SCL values were tested including:
c       SCL(1) = HMF2 * 0.15 - 27.1     HXL(2) = 200.
c       HXL(2) = HMF1 + SCL(2)          HXL(3) = 140.
c       SCL(3) = 5.                     HXL(4) = HME + SCL(4)
c       HXL(4) = 105.
c
7711            NUMCON = 7
      HXL(1) = HHALF
        HXL1T  = 0.4 * HMF2 + 30.
      HXL(2) = ( HMF2 + HV1 ) / 2.
      HXL(4) = HME
      XX(6) = HV2
      XX(7) = HME
      YY(6) = 0.0
      WW(1) = 1.
      WW(3) = 3.
      WW(5) = 0.5
      WW(6) = 50.
      WW(7) = 500.
      HFFF=HHALF
      XFFF=XHALF
c
C are valley-top and bottomside point compatible ? -------------
C
7722    IF((HV1-HFFF)*(XNME-XFFF).LT.0.0) WW(2)=0.5
      IF(HV1.LE.HV2+5.0) WW(2)=0.5
c
C DETERMINE AMPLITUDES-----------------------------------------
C
       NC0=NUMCON-NC1
       IQUAL=0
2299        CALL LSKNM(NUMLAY,NUMCON,NC0,NC1,HMF2,SCL,HXL,WW,XX,YY,
     &          AMP,SSIN)
         IF(IQUAL.gt.0) GOTO 1937
       IF((ABS(AMP(1)).GT.10.0).OR.(SSIN)) THEN
         IQUAL=1
         HXL(1)=HXL1T
         GOTO 2299
         ENDIF
1937        IF(SSIN) IQUAL=2
       RETURN
       END
c
c
        subroutine ioncom(h,z,f,fs,t,cn)
c---------------------------------------------------------------
c ion composition model
c A.D. Danilov and A.P. Yaichnikov, A New Model of the Ion
c   Composition at 75 to 1000 km for IRI, Adv. Space Res. 5, #7,
c   75-79, 107-108, 1985
c
c       h       altitude in km
c       z       solar zenith angle in radians
c       f       latitude in radians
c       fs      10.7cm solar radio flux
c       t       season (decimal month)
c       cn(1)   O+  relative density in percent
c       cn(2)   H+  relative density in percent
c       cn(3)   N+  relative density in percent
c       cn(4)   He+ relative density in percent
c       cn(5)   NO+ relative density in percent
c       cn(6)   O2+ relative density in percent
c       cn(7)   cluster ions  relative density in percent
c---------------------------------------------------------------
c
        dimension       cn(7),cm(7),hm(7),alh(7),all(7),beth(7),
     &                  betl(7),p(5,6,7),var(6),po(5,6),ph(5,6),
     &                  pn(5,6),phe(5,6),pno(5,6),po2(5,6),pcl(5,6)

        common  /argexp/argmax
        data po/4*0.,98.5,4*0.,320.,4*0.,-2.59E-4,2.79E-4,-3.33E-3,
     &          -3.52E-3,-5.16E-3,-2.47E-2,4*0.,-2.5E-6,1.04E-3,
     &          -1.79E-4,-4.29E-5,1.01E-5,-1.27E-3/
        data ph/-4.97E-7,-1.21E-1,-1.31E-1,0.,98.1,355.,-191.,
     &          -127.,0.,2040.,4*0.,-4.79E-6,-2.E-4,5.67E-4,
     &          2.6E-4,0.,-5.08E-3,10*0./
        data pn/7.6E-1,-5.62,-4.99,0.,5.79,83.,-369.,-324.,0.,593.,
     &          4*0.,-6.3E-5,-6.74E-3,-7.93E-3,-4.65E-3,0.,-3.26E-3,
     &          4*0.,-1.17E-5,4.88E-3,-1.31E-3,-7.03E-4,0.,-2.38E-3/
        data phe/-8.95E-1,6.1,5.39,0.,8.01,4*0.,1200.,4*0.,-1.04E-5,
     &          1.9E-3,9.53E-4,1.06E-3,0.,-3.44E-3,10*0./
        data pno/-22.4,17.7,-13.4,-4.88,62.3,32.7,0.,19.8,2.07,115.,
     &          5*0.,3.94E-3,0.,2.48E-3,2.15E-4,6.67E-3,5*0.,
     &          -8.4E-3,0.,-3.64E-3,2.E-3,-2.59E-2/
        data po2/8.,-12.2,9.9,5.8,53.4,-25.2,0.,-28.5,-6.72,120.,
     &          5*0.,-1.4E-2,0.,-9.3E-3,3.3E-3,2.8E-2,5*0.,4.25E-3,
     &          0.,-6.04E-3,3.85E-3,-3.64E-2/
        data pcl/4*0.,100.,4*0.,75.,10*0.,4*0.,-9.04E-3,-7.28E-3,
     &          2*0.,3.46E-3,-2.11E-2/

        DO 8 I=1,5
        DO 8 J=1,6
                p(i,j,1)=po(i,j)
                p(i,j,2)=ph(i,j)
                p(i,j,3)=pn(i,j)
                p(i,j,4)=phe(i,j)
                p(i,j,5)=pno(i,j)
                p(i,j,6)=po2(i,j)
                p(i,j,7)=pcl(i,j)
8       continue

        s=0.
        do 5 i=1,7
          do 7 j=1,6
                var(j) = p(1,j,i)*cos(z) + p(2,j,i)*cos(f) +
     &                   p(3,j,i)*cos(0.013*(300.-fs)) +
     &                   p(4,j,i)*cos(0.52*(t-6.)) + p(5,j,i)
7         continue
          cm(i)  = var(1)
          hm(i)  = var(2)
          all(i) = var(3)
          betl(i)= var(4)
          alh(i) = var(5)
          beth(i)= var(6)
          hx=h-hm(i)
          if(hx) 1,2,3
1               arg = hx * (hx * all(i) + betl(i))
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
                goto 4
2               cn(i) = cm(i)
                goto 4
3               arg = hx * (hx * alh(i) + beth(i))
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
4         continue
          if(cn(i).LT.0.005*cm(i)) cn(i)=0.
          if(cn(i).GT.cm(i)) cn(i)=cm(i)
          s=s+cn(i)
5       continue
        do 6 i=1,7
6               cn(i)=cn(i)/s*100.
        return
        end
c
c

      Subroutine ionco2(h,z,it,F,R1,R2,R3,R4)
*----------------------------------------------------------------
*     INPUT DATA :
*      h -  altitude in km
*      z -  solar zenith angle in degree
*      F -  10.7cm solar radio flux
*      it-  season (month)
*     OUTPUT DATA :
*     R1 -  NO+ concentration (in percent)
*     R2 -  O2+ concentration (in percent)
*     R3 -  Cb+ concentration (in percent)
*     R4 -  O+  concentration (in percent)
*-----------------------------------------------------------------
      dimension j1ms70(7),j2ms70(7),h1s70(13,7),h2s70(13,7),
     *       R1ms70(13,7),R2ms70(13,7),rk1ms70(13,7),rk2ms70(13,7),
     *       j1ms140(7),j2ms140(7),h1s140(13,7),h2s140(13,7),
     *       R1ms140(13,7),R2ms140(13,7),rk1ms140(13,7),rk2ms140(13,7),
     *       j1mw70(7),j2mw70(7),h1w70(13,7),h2w70(13,7),
     *       R1mw70(13,7),R2mw70(13,7),rk1mw70(13,7),rk2mw70(13,7),
     *       j1mw140(7),j2mw140(7),h1w140(13,7),h2w140(13,7),
     *       R1mw140(13,7),R2mw140(13,7),rk1mw140(13,7),rk2mw140(13,7),
     *       j1mr70(7),j2mr70(7),h1r70(13,7),h2r70(13,7),moind(12),
     *       R1mr70(13,7),R2mr70(13,7),rk1mr70(13,7),rk2mr70(13,7),
     *       j1mr140(7),j2mr140(7),h1r140(13,7),h2r140(13,7),
     *       R1mr140(13,7),R2mr140(13,7),rk1mr140(13,7),rk2mr140(13,7)
      data moind/1,1,2,2,3,3,3,3,2,2,1,1/
      data j1ms70/11,11,10,10,11,9,11/
      data j2ms70/13,11,10,11,11,9,11/
      data h1s70/75,85,90,95,100,120,130,200,220,250,270,0,0,
     *        75,85,90,95,100,120,130,200,220,250,270,0,0,
     *        75,85,90,95,100,115,200,220,250,270,0,0,0,
     *        75,80,95,100,120,140,200,220,250,270,0,0,0,
     *        75,80,95,100,120,150,170,200,220,250,270,0,0,
     *        75,80,95,100,140,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,110,145,200,220,250,270,0,0/
      data h2s70/75,80,90,95,100,120,130,140,150,200,220,250,270,
     *        75,80,90,95,100,120,130,200,220,250,270,0,0,
     *        75,80,90,95,100,115,200,220,250,270,0,0,0,
     *        75,80,95,100,120,140,150,200,220,250,270,0,0,
     *        75,80,95,100,120,150,170,200,220,250,270,0,0,
     *        75,80,95,100,140,200,220,250,270,0,0,0,0,
     *        75,80,90,95,100,110,145,200,220,250,270,0,0/
      data R1ms70/6,30,60,63,59,59,66,52,20,4,2,0,0,
     *         6,30,60,63,69,62,66,52,20,4,2,0,0,
     *         6,30,60,63,80,68,53,20,4,2,0,0,0,
     *         4,10,60,85,65,65,52,25,12,4,0,0,0,
     *         4,10,60,89,72,60,60,52,30,20,10,0,0,
     *         4,10,60,92,68,54,40,25,13,0,0,0,0,
     *         1,8,20,60,95,93,69,65,45,30,20,0,0/
      data R2ms70/4,10,30,32,41,41,32,29,34,28,15,3,1,
     *         4,10,30,32,31,38,32,28,15,3,1,0,0,
     *         4,10,30,32,20,32,28,15,3,1,0,0,0,
     *         2,6,30,15,35,30,34,26,19,8,3,0,0,
     *         2,6,30,11,28,38,29,29,25,12,5,0,0,
     *         2,6,30,8,32,30,20,14,8,0,0,0,0,
     *         1,2,10,20,5,7,31,23,18,15,10,0,0/
      data rk1ms70/2.4,6.,.6,-.8,0,.7,-.2,-1.6,-.533,-.1,-.067,0,0,
     *         2.4,6.,.6,1.2,-.35,.4,-.2,-1.6,-.533,-.1,-.067,0,0,
     *         2.4,6.,.6,3.4,-.8,-.176,-1.65,-.533,-.1,-.067,0,0,0,
     *         1.2,3.333,5.,-1.,0,-.216,-1.35,-.433,-.4,-.1,0,0,0,
     *         1.2,3.333,5.8,-.85,-.4,0,-.267,-1.1,-.333,-.4,-.2,0,0,
     *         1.2,3.333,6.4,-.6,-.233,-.7,-.5,-.6,-.267,0,0,0,0,
     *         1.4,2.4,4.,7.,-.2,-.686,-.072,-1.,-.5,-.5,-.5,0,0/
      data rk2ms70/1.2,2.,.4,1.8,0,-.9,-.3,.5,-.12,-.65,-.4,-.1,-.033,
     *         1.2,2.,.4,-.2,.35,-.6,-.057,-.65,-.4,-.1,-.033,0,0,
     *         1.2,2.,.4,-2.4,.8,-.047,-.65,-.4,-.1,-.033,0,0,0,
     *         .8,1.6,-3.,1.,-.25,.4,-.16,-.35,-.367,-.25,-.1,0,0,
     *         .8,1.6,-3.8,.85,.333,-.45,0,-.2,-.433,-.35,-.1,0,0,
     *         .8,1.6,-4.4,.6,-.033,-.5,-.2,-.3,-.2,0,0,0,0,
     *         .2,.8,2.,-3.,.2,.686,-.145,-.25,-.1,-.25,-.2,0,0/
      data j1ms140/11,11,10,10,9,9,12/
      data j2ms140/11,11,10,9,10,10,12/
      data h1s140/75,85,90,95,100,120,130,140,200,220,250,0,0,
     *        75,85,90,95,100,120,130,140,200,220,250,0,0,
     *        75,85,90,95,100,120,140,200,220,250,0,0,0,
     *        75,80,95,100,120,140,200,220,250,270,0,0,0,
     *        75,80,95,100,120,200,220,250,270,0,0,0,0,
     *        75,80,95,100,130,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,110,140,180,200,220,250,270,0/
      data h2s140/75,80,90,95,100,120,130,155,200,220,250,0,0,
     *        75,80,90,95,100,120,130,160,200,220,250,0,0,
     *        75,80,90,95,100,120,165,200,220,250,0,0,0,
     *        75,80,95,100,120,180,200,250,270,0,0,0,0,
     *        75,80,95,100,120,160,200,220,250,270,0,0,0,
     *        75,80,95,100,130,160,200,220,250,270,0,0,0,
     *        75,80,90,95,100,110,140,180,200,220,250,270,0/
      data R1ms140/6,30,60,63,59,59,66,66,38,14,1,0,0,
     *         6,30,60,63,69,62,66,66,38,14,1,0,0,
     *         6,30,60,63,80,65,65,38,14,1,0,0,0,
     *         4,10,60,85,66,66,38,22,9,1,0,0,0,
     *         4,10,60,89,71,42,26,17,10,0,0,0,0,
     *         4,10,60,93,71,48,35,22,10,0,0,0,0,
     *         1,8,20,60,95,93,72,60,58,40,26,13,0/
      data R2ms140/4,10,30,32,41,41,30,30,10,6,1,0,0,
     *         4,10,30,32,31,38,31,29,9,6,1,0,0,
     *         4,10,30,32,20,35,26,9,6,1,0,0,0,
     *         2,6,30,15,34,24,10,5,1,0,0,0,0,
     *         2,6,30,11,28,37,21,14,8,5,0,0,0,
     *         2,6,30,7,29,36,29,20,13,5,0,0,0,
     *         1,2,10,20,5,7,28,32,28,20,14,7,0/
      data rk1ms140/2.4,6.,.6,-.8,0,.7,0,-.467,-1.2,-.433,0,0,0,
     *         2.4,6.,.6,1.2,-.35,.4,0,-.467,-1.2,-.433,0,0,0,
     *         2.4,6.,.6,3.4,-.75,0,-.45,-1.2,-.433,0,0,0,0,
     *         1.2,3.333,5.,-.95,0,-.467,-.8,-.433,-.4,0,0,0,0,
     *         1.2,3.333,5.8,-.9,-.363,-.8,-.3,-.35,-.3,0,0,0,0,
     *         1.2,3.333,6.6,-.733,-.329,-.65,-.433,-.6,-.267,0,0,0,0,
     *         1.4,2.4,4.,7.,-.2,-.7,-.3,-.1,-.9,-.467,-.65,-.333,0/
      data rk2ms140/1.2,2.,.4,1.8,0,-1.1,0,-.444,-.2,-.166,0,0,0,
     *         1.2,2.,.4,-.2,.35,-.7,-.067,-.5,-.15,-.166,0,0,0,
     *         1.2,2.,.4,-2.4,.75,-.2,-.486,-.15,-.166,0,0,0,0,
     *         .8,1.6,-3.,.95,-.167,-.7,-.1,-.2,0,0,0,0,0,
     *         .8,1.6,-3.8,.85,.225,-.4,-.35,-.2,-.15,-.133,0,0,0,
     *         .8,1.6,-4.6,.733,.233,-.175,-.45,-.233,-.4,-.1,0,0,0,
     *         .2,.8,2.,-3.,.2,.7,.1,-.2,-.4,-.2,-.35,-.167,0/
      data j1mr70/12,12,12,9,10,11,13/
      data j2mr70/9,9,10,13,12,11,11/
      data h1r70/75,80,90,95,100,120,140,180,200,220,250,270,0,
     *        75,80,90,95,100,120,145,180,200,220,250,270,0,
     *        75,80,90,95,100,120,145,180,200,220,250,270,0,
     *        75,95,100,110,140,180,200,250,270,0,0,0,0,
     *        75,95,125,150,185,195,200,220,250,270,0,0,0,
     *        75,95,100,150,160,170,190,200,220,250,270,0,0,
     *        75,80,85,95,100,140,160,170,190,200,220,250,270/
      data h2r70/75,95,100,120,180,200,220,250,270,0,0,0,0,
     *        75,95,100,120,180,200,220,250,270,0,0,0,0,
     *        75,95,100,120,130,190,200,220,250,270,0,0,0,
     *        75,80,85,95,100,110,130,180,190,200,220,250,270,
     *        75,80,85,95,100,125,150,190,200,220,250,270,0,
     *        75,80,85,95,100,150,190,200,220,250,270,0,0,
     *        75,85,95,100,140,180,190,200,220,250,270,0,0/
      data R1mr70/13,17,57,57,30,53,58,38,33,14,6,2,0,
     *         13,17,57,57,37,56,56,38,33,14,6,2,0,
     *         13,17,57,57,47,58,55,37,33,14,6,2,0,
     *         5,65,54,58,58,38,33,9,1,0,0,0,0,
     *         5,65,65,54,40,40,45,26,17,10,0,0,0,
     *         5,65,76,56,57,48,44,51,35,22,10,0,0,
     *         3,11,35,75,90,65,63,54,54,50,40,26,13/
      data R2mr70/7,43,70,47,15,17,10,4,0,0,0,0,0,
     *         7,43,63,44,17,17,10,4,0,0,0,0,0,
     *         7,43,53,42,42,13,17,10,4,0,0,0,0,
     *         3,5,26,34,46,42,41,23,16,16,10,1,0,
     *         3,5,26,34,35,35,42,25,22,14,8,5,0,
     *         3,5,26,34,24,41,31,26,20,13,5,0,0,
     *         3,15,15,10,35,35,30,34,20,14,7,0,0/
      data rk1mr70/.8,4.,0,-5.4,1.15,.25,-.5,-.25,-.95,-.267,-.2,
     *             -.067,0,
     *         .8,4.,0,-4.,.95,0,-.514,-.25,-.95,-.267,-.2,-.067,0,
     *         .8,4.,0,-2.,.55,-.12,-.514,-.2,-.95,-.267,-.2,-.067,0,
     *         3.,-2.2,.4,0,-.5,-.25,-.48,-.4,-.033,0,0,0,0,
     *         3.,0,-.44,-.466,0,1.0,-.95,-.3,-.35,-.3,0,0,0,
     *         3.,2.2,-.4,0.1,-.9,-.2,.7,-.8,-.433,-.6,-.267,0,0,
     *         1.6,4.8,4.,3.,-.625,-.1,-.9,0,-.4,-.5,-.467,-.65,-.3/
      data rk2mr70/1.8,5.4,-1.15,-.533,.1,-.35,-.2,-.2,0,0,0,0,0,
     *         1.8,4.,-.95,-.45,0,-.35,-.2,-.2,0,0,0,0,0,
     *         1.8,2.,-.55,0,-.483,.4,-.35,-.2,-.2,0,0,0,0,
     *         .4,4.2,.8,2.4,-.4,-.05,-.36,-.7,0,-.3,-.3,-.05,0,
     *         .4,4.2,.8,.2,0,.28,-.425,-.3,-.4,-.2,-.15,-.133,0,
     *         .4,4.2,.8,-2.,.34,-.25,-.5,-.3,-.233,-.4,-.1,0,0,
     *         1.2,0,-1.,.625,0,-.5,.4,-.7,-.2,-.35,-.167,0,0/
      data j1mr140/12,12,11,12,9,9,13/
      data j2mr140/10,9,10,12,13,13,12/
      data h1r140/75,80,90,95,100,115,130,145,200,220,250,270,0,
     *        75,80,90,95,100,110,120,145,200,220,250,270,0,
     *        75,80,90,95,100,115,150,200,220,250,270,0,0,
     *        75,95,100,120,130,140,150,190,200,220,250,270,0,
     *        75,95,120,150,190,200,220,250,270,0,0,0,0,
     *        75,95,100,145,190,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,120,160,170,190,200,220,250,270/
      data h2r140/75,95,100,115,130,175,200,220,250,270,0,0,0,
     *        75,95,100,110,175,200,220,250,270,0,0,0,0,
     *        75,95,100,115,130,180,200,220,250,270,0,0,0,
     *        75,80,85,95,100,120,130,190,200,220,250,270,0,
     *        75,80,85,95,100,120,140,160,190,200,220,250,270,
     *        75,80,85,95,100,145,165,180,190,200,220,250,270,
     *        75,85,95,100,120,145,170,190,200,220,250,270,0/
      data R1mr140/13,17,57,57,28,51,56,56,12,8,1,0,0,
     *         13,17,57,57,36,46,55,56,10,8,1,0,0,
     *         13,17,57,57,46,56,55,12,8,1,0,0,0,
     *         5,65,54,59,56,56,53,23,16,13,3,1,0,
     *         5,65,65,54,29,16,16,10,2,0,0,0,0,
     *         5,65,76,58,36,25,20,12,7,0,0,0,0,
     *         3,11,35,75,91,76,58,49,45,32,28,20,12/
      data R2mr140/7,43,72,49,44,14,7,4,1,0,0,0,0,
     *         7,43,64,51,14,7,4,1,0,0,0,0,0,
     *         7,43,54,44,44,13,7,4,1,0,0,0,0,
     *         3,5,26,34,46,41,44,9,11,7,2,1,0,
     *         3,5,26,34,35,35,40,40,16,14,9,5,2,
     *         3,5,26,34,24,40,40,32,19,20,10,7,3,
     *         3,15,15,9,24,35,40,28,28,20,10,8,0/
      data rk1mr140/.8,4.,0,-5.8,1.533,.333,0,-.8,-.2,-.233,-.05,0,0,
     *         .8,4.,0,-4.2,1.3,.6,.04,-.836,-.1,-.233,-.05,0,0,
     *         .8,4.,0,-2.2,.667,-.029,-.86,-.2,-.233,-.05,0,0,0,
     *         3.,-2.2,.25,-.3,0,-.3,-.75,-.7,-.15,-.333,-.1,-.033,0,
     *         3.,0,-.367,-.625,-1.3,0,-.2,-.4,-.067,0,0,0,0,
     *         3.,2.2,-.4,-.489,-1.1,-.25,-.267,-.25,-.2,0,0,0,0,
     *         1.6,4.8,4.,3.2,-.75,-.45,-.9,-.2,-1.3,-.2,-.267,-.4,-.3/
      data rk2mr140/1.8,5.8,-1.533,-.333,-.667,-.28,-.15,-.1,-.05,
     *              0,0,0,0,
     *         1.8,4.2,-1.3,-.569,-.28,-.15,-.1,-.05,0,0,0,0,0,
     *         1.8,2.2,-.667,0,-.62,-.3,-.15,-.1,-.05,0,0,0,0,
     *         .4,4.2,.8,2.4,-.25,.3,-.583,.2,-.2,-.167,-.05,-.033,0,
     *         .4,4.2,.8,.02,0,.25,0,-.6,-.2,-.25,-.133,-.15,-.067,
     *         .4,4.2,.8,-2.,.356,0,-.533,-1.3,.1,-.5,-.1,-.2,-.1,
     *         1.2,0,-1.2,.75,.44,.2,-.6,0,-.4,-.333,-.1,-.2,0/
      data j1mw70/13,13,13,13,9,8,9/
      data j2mw70/10,10,11,11,9,8,11/
      data h1w70/75,80,85,95,100,110,125,145,180,200,220,250,270,
     *        75,80,85,95,100,110,120,150,180,200,220,250,270,
     *        75,80,85,95,100,110,120,155,180,200,220,250,270,
     *        75,80,90,100,110,120,140,160,190,200,220,250,270,
     *        75,80,90,110,150,200,220,250,270,0,0,0,0,
     *        75,80,90,100,150,200,250,270,0,0,0,0,0,
     *        75,80,90,100,120,130,140,200,270,0,0,0,0/
      data h2w70/75,90,95,100,110,125,190,200,250,270,0,0,0,
     *        75,90,95,100,110,125,190,200,250,270,0,0,0,
     *        75,90,95,100,110,120,145,190,200,250,270,0,0,
     *        75,80,95,100,110,120,150,200,220,250,270,0,0,
     *        75,80,90,95,110,145,200,250,270,0,0,0,0,
     *        75,80,90,100,140,150,200,250,0,0,0,0,0,
     *        75,80,85,90,100,120,130,140,160,200,270,0,0/
      data R1mw70/28,35,65,65,28,44,46,50,25,25,10,5,0,
     *         28,35,65,65,36,49,47,47,25,25,10,5,0,
     *         28,35,65,65,48,54,51,43,25,25,10,5,0,
     *         16,24,66,54,58,50,50,38,25,25,10,5,0,
     *         16,24,66,66,46,30,20,6,3,0,0,0,0,
     *         16,24,66,76,49,32,12,7,0,0,0,0,0,
     *         6,19,67,91,64,68,60,40,12,0,0,0,0/
      data R2mw70/5,35,35,72,56,54,12,12,2,0,0,0,0,
     *         5,35,35,64,51,53,12,12,2,0,0,0,0,
     *         5,35,35,52,46,49,41,12,12,2,0,0,0,
     *         4,10,40,46,42,50,41,12,7,2,0,0,0,
     *         4,10,30,34,34,51,14,4,2,0,0,0,0,
     *         4,10,30,24,45,48,20,5,0,0,0,0,0,
     *         2,6,17,23,9,36,32,40,40,20,6,0,0/
      data rk1mw70/1.4,6.,0,-7.4,1.6,.133,.2,-.714,0,-.75,-.167,-.25,0,
     *         1.4,6.,0,-5.8,1.3,-.2,0,-.733,0,-.75,-.167,-.25,0,
     *         1.4,6.,0,-3.4,.6,-.3,-.229,-.72,0,-.75,-.167,-.25,0,
     *         1.6,4.2,-1.2,.4,-.8,0,-.6,-.433,0,-.75,-.167,-.25,0,
     *         1.6,4.2,0,-.5,-.32,-.5,-.467,-.15,-.1,0,0,0,0,
     *         1.6,4.2,1.,-.54,-.34,-.4,-.25,-.2,0,0,0,0,0,
     *         2.6,4.8,2.4,-1.35,.4,-.8,-.333,-.4,-.3,0,0,0,0/
      data rk2mw70/2.,0,7.4,-1.6,-.133,-.646,0,-.2,-.1,0,0,0,0,
     *         2.,0,5.8,-1.3,.133,-.631,0,-.2,-.1,0,0,0,0,
     *         2.,0,3.4,-.6,.3,-.32,-.644,0,-.2,-.1,0,0,0,
     *         1.2,2.,1.2,-.4,.8,-.3,-.58,-.25,-.167,-.1,0,0,0,
     *         1.2,2.,.8,0,.486,-.673,-.2,-.1,-.066,0,0,0,0,
     *         1.2,2.,-.6,.525,.3,-.56,-.3,-.1,0,0,0,0,0,
     *         .8,2.2,1.2,-1.4,1.35,-.4,.8,0,-.5,-.2,-.167,0,0/
      data j1mw140/12,11,11,11,11,10,12/
      data j2mw140/10,11,11,11,11,10,12/
      data h1w140/75,80,85,95,100,110,125,145,190,200,220,250,0,
     *        75,80,85,95,100,110,120,150,190,220,250,0,0,
     *        75,80,85,95,100,110,120,155,190,220,250,0,0,
     *        75,80,90,100,110,120,140,160,190,220,250,0,0,
     *        75,80,90,110,150,160,190,200,220,250,270,0,0,
     *        75,80,90,100,150,160,190,200,250,270,0,0,0,
     *        75,80,90,100,120,130,140,160,190,200,250,270,0/
      data h2w140/75,90,95,100,110,125,190,200,220,250,0,0,0,
     *        75,90,95,100,110,120,125,190,200,220,250,0,0,
     *        75,90,95,100,110,120,145,190,200,220,250,0,0,
     *        75,80,95,100,110,120,150,190,200,220,250,0,0,
     *        75,80,90,95,110,145,190,200,220,250,270,0,0,
     *        75,80,90,100,140,150,200,220,250,270,0,0,0,
     *        75,80,85,90,100,120,130,140,160,180,200,220,0/
      data R1mw140/28,35,65,65,28,44,46,50,9,6,2,0,0,
     *         28,35,65,65,36,49,47,47,8,2,0,0,0,
     *         28,35,65,65,48,54,51,43,8,2,0,0,0,
     *         16,24,66,54,58,50,50,42,8,2,0,0,0,
     *         16,24,66,66,46,49,9,10,7,2,0,0,0,
     *         16,24,66,76,49,54,10,14,4,1,0,0,0,
     *         6,19,67,91,64,68,60,58,11,20,5,2,0/
      data R2mw140/5,35,35,72,56,54,5,5,1,0,0,0,0,
     *         5,35,35,64,51,53,53,5,5,1,0,0,0,
     *         5,35,35,52,46,49,41,5,5,1,0,0,0,
     *         4,10,40,46,42,50,41,5,5,1,0,0,0,
     *         4,10,30,34,34,51,10,5,3,1,0,0,0,
     *         4,10,30,24,45,48,4,2,1,0,0,0,0,
     *         2,6,17,23,9,36,32,40,39,29,1,0,0/
      data rk1mw140/1.4,6.,0,-7.4,1.6,.133,.2,-.911,-.3,-.2,-.066,0,0,
     *         1.4,6.,0,-5.8,1.3,-.2,0,-.975,-.2,-.066,0,0,0,
     *         1.4,6.,0,-3.4,.6,-.3,-.229,-1.,-.2,-.066,0,0,0,
     *         1.6,4.2,-1.2,.4,-.8,0,-.4,-1.133,-.2,-.066,0,0,0,
     *         1.6,4.2,0,-.5,.3,-1.133,.1,-.15,-.166,-.1,0,0,0,
     *         1.6,4.2,1.,-.54,.5,-1.466,.4,-.2,-.15,-.0333,0,0,0,
     *         2.6,4.8,2.4,-1.35,.4,-.8,-.1,-1.566,.9,-.3,-.15,-.05,0/
      data rk2mw140/2.,0,7.4,-1.6,-.133,-.754,0,-.2,-.033,0,0,0,0,
     *         2.,0,5.8,-1.3,.2,0,-.738,0,-.2,-.033,0,0,0,
     *         2.,0,3.4,-.6,.3,-.32,-.8,0,-.2,-.033,0,0,0,
     *         1.2,2.,1.2,-.4,.8,-.3,-.9,0,-.2,-.033,0,0,0,
     *         1.2,2.,.8,0,.486,-.911,-.5,-.1,-.066,-.05,0,0,0,
     *         1.2,2.,-.6,.525,.3,-.88,-.1,-.033,-.05,0,0,0,0,
     *         .8,2.2,1.2,-1.4,1.35,-.4,.8,-.05,-.5,-1.4,-.05,0,0/

       if(z.lt.20)z=20
       if(z.gt.90)z=90
       itzz=moind(it)
      if(itzz.eq.1)then
       if(f.lt.140)then
      Call aprok(j1mw70,j2mw70,h1w70,h2w70,R1mw70,R2mw70,
     *                rk1mw70,rk2mw70,h,z,R1,R2)
      R170=R1
      R270=R2
       endif
       if(f.gt.70)then
      Call aprok(j1mw140,j2mw140,h1w140,h2w140,R1mw140,R2mw140,
     *                rk1mw140,rk2mw140,h,z,R1,R2)
      R1140=R1
      R2140=R2
       endif
       if((f.gt.70).and.(f.lt.140))then
         R1=R170+(R1140-R170)*(f-70)/70
         R2=R270+(R2140-R270)*(f-70)/70
       endif
      endif
      if(itzz.eq.3)then
       if(f.lt.140)then
      Call aprok(j1ms70,j2ms70,h1s70,h2s70,R1ms70,R2ms70,
     *                rk1ms70,rk2ms70,h,z,R1,R2)
      R170=R1
      R270=R2
       endif
       if(f.gt.70)then
      Call aprok(j1ms140,j2ms140,h1s140,h2s140,R1ms140,R2ms140,
     *                rk1ms140,rk2ms140,h,z,R1,R2)
      R1140=R1
      R2140=R2
       endif
       if((f.gt.70).and.(f.lt.140))then
      R1=R170+(R1140-R170)*(f-70)/70
      R2=R270+(R2140-R270)*(f-70)/70
       endif
      endif
      if(itzz.eq.2)then
       if(f.lt.140)then
      Call aprok(j1mr70,j2mr70,h1r70,h2r70,R1mr70,R2mr70,
     *                rk1mr70,rk2mr70,h,z,R1,R2)
      R170=R1
      R270=R2
       endif
       if(f.gt.70)then
      Call aprok(j1mr140,j2mr140,h1r140,h2r140,R1mr140,R2mr140,
     *                rk1mr140,rk2mr140,h,z,R1,R2)
      R1140=R1
      R2140=R2
       endif
       if((f.gt.70).and.(f.lt.140))then
      R1=R170+(R1140-R170)*(f-70)/70
      R2=R270+(R2140-R270)*(f-70)/70
       endif
      endif
      R3=0
      R4=0
      if (h.lt.100) R3=100-(R1+R2)
      if (h.ge.100) R4=100-(R1+R2)
       if(R3.lt.0) R3=0
       if(R4.lt.0) R4=0
      R1=ANINT(R1)
      R2=ANINT(R2)
      R3=ANINT(R3)
      R4=ANINT(R4)
 300   continue
      return
      end
c
c
      Subroutine aprok(j1m,j2m,h1,h2,R1m,R2m,
     *                rk1m,rk2m,h,z,R1,R2)
      dimension zm(7),j1m(7),j2m(7),h1(13,7),h2(13,7),R1m(13,7),
     *          R2m(13,7),rk1m(13,7),rk2m(13,7)
      data zm/20,40,60,70,80,85,90/

       j1=1
       j2=1
       i1=1
       do 1 i=1,7
       i1=i
      if(z.eq.zm(i)) j1=0
c        if(z.le.zm(i)) exit
        if(z.le.zm(i)) goto 1234
 1     continue
 11    continue
      i2=1
       do 2 i=2,j1m(i1)
        i2=i-1
c          if(h.lt.h1(i,i1)) exit
          if(h.lt.h1(i,i1)) goto 1234
        i2=j1m(i1)
 2       continue
        i3=1
       do 3 i=2,j2m(i1)
        i3=i-1
          if(h.lt.h2(i,i1)) goto 1234
c          if(h.lt.h2(i,i1)) exit
        i3=j2m(i1)
 3       continue
      R01=R1m(i2,i1)
      R02=R2m(i3,i1)
      rk1=rk1m(i2,i1)
      rk2=rk2m(i3,i1)
      h01=h1(i2,i1)
      h02=h2(i3,i1)
      R1=R01+rk1*(h-h01)
      R2=R02+rk2*(h-h02)
      if(j1.eq.1)then
      j1=0
      j2=0
      i1=i1-1
      R11=R1
      R12=R2
      goto 11
      endif
      if(j2.eq.0)then
      rk=(z-zm(i1))/(zm(i1+1)-zm(i1))
      R1=R1+(R11-R1)*rk
      R2=R2+(R12-R2)*rk
      endif
1234    continue
       return
       end
c
c
      subroutine ioncom_new(h,xhi,xlati,cov,zmosea,dion)
c-------------------------------------------------------
c       see IONCO1 for explanation of i/o parameters
c       NOTE: xhi,xlati are in DEGREES !!!
c       NOTE: zmosea is the seasonal northern month, so
c             for the southern hemisphere zmosea=month+6
c-------------------------------------------------------
        dimension       dion(7),dup(4),diont(7)
      common  /const/ umr,PI

        do 1122 i=1,7
1122                diont(i)=0.

      if (h.gt.300.) then
      xhirad=xhi*umr
      xlarad=xlati*umr
                call ionco1(h,xhirad,xlarad,cov,zmosea,dup)
                diont(1)=dup(1)
                diont(2)=dup(2)
                diont(3)=dup(3)
                diont(4)=dup(4)
      else
      monsea=int(zmosea)
                call ionco2(h,xhi,monsea,cov,rno,ro2,rcl,ro)
                diont(5)=rno
                diont(6)=ro2
                diont(7)=rcl
                diont(1)=ro
      endif
        do 1 i=1,7
                dion(i)=diont(i)
1               continue
      return
      end
c
c
        subroutine ionco1(h,z,f,fs,t,cn)
c---------------------------------------------------------------
c ion composition model
c A.D. Danilov and A.P. Yaichnikov, A New Model of the Ion
c   Composition at 75 to 1000 km for IRI, Adv. Space Res. 5, #7,
c   75-79, 107-108, 1985
c
c       h       altitude in km
c       z       solar zenith angle in radians
c       f       latitude in radians
c       fs      10.7cm solar radio flux
c       t       season (decimal month)
c       cn(1)   O+  relative density in percent
c       cn(2)   H+  relative density in percent
c       cn(3)   N+  relative density in percent
c       cn(4)   He+ relative density in percent
c Please note: molecular ions are now computed in IONCO2
c       [cn(5)   NO+ relative density in percent
c       [cn(6)   O2+ relative density in percent
c       [cn(7)   cluster ions  relative density in percent
c---------------------------------------------------------------
c
c        dimension       cn(7),cm(7),hm(7),alh(7),all(7),beth(7),
c     &                  betl(7),p(5,6,7),var(6),po(5,6),ph(5,6),
c     &                  pn(5,6),phe(5,6),pno(5,6),po2(5,6),pcl(5,6)
      dimension       cn(4),cm(4),hm(4),alh(4),all(4),beth(4),
     &                  betl(4),p(5,6,4),var(6),po(5,6),ph(5,6),
     &                  pn(5,6),phe(5,6)

      common  /argexp/argmax
      data po/4*0.,98.5,4*0.,320.,4*0.,-2.59E-4,2.79E-4,-3.33E-3,
     &          -3.52E-3,-5.16E-3,-2.47E-2,4*0.,-2.5E-6,1.04E-3,
     &          -1.79E-4,-4.29E-5,1.01E-5,-1.27E-3/
      data ph/-4.97E-7,-1.21E-1,-1.31E-1,0.,98.1,355.,-191.,
     &          -127.,0.,2040.,4*0.,-4.79E-6,-2.E-4,5.67E-4,
     &          2.6E-4,0.,-5.08E-3,10*0./
      data pn/7.6E-1,-5.62,-4.99,0.,5.79,83.,-369.,-324.,0.,593.,
     &          4*0.,-6.3E-5,-6.74E-3,-7.93E-3,-4.65E-3,0.,-3.26E-3,
     &          4*0.,-1.17E-5,4.88E-3,-1.31E-3,-7.03E-4,0.,-2.38E-3/
      data phe/-8.95E-1,6.1,5.39,0.,8.01,4*0.,1200.,4*0.,-1.04E-5,
     &          1.9E-3,9.53E-4,1.06E-3,0.,-3.44E-3,10*0./
c       data pno/-22.4,17.7,-13.4,-4.88,62.3,32.7,0.,19.8,2.07,115.,
c    &          5*0.,3.94E-3,0.,2.48E-3,2.15E-4,6.67E-3,5*0.,
c    &          -8.4E-3,0.,-3.64E-3,2.E-3,-2.59E-2/
c       data po2/8.,-12.2,9.9,5.8,53.4,-25.2,0.,-28.5,-6.72,120.,
c    &          5*0.,-1.4E-2,0.,-9.3E-3,3.3E-3,2.8E-2,5*0.,4.25E-3,
c    &          0.,-6.04E-3,3.85E-3,-3.64E-2/
c       data pcl/4*0.,100.,4*0.,75.,10*0.,4*0.,-9.04E-3,-7.28E-3,
c    &          2*0.,3.46E-3,-2.11E-2/

      DO 8 I=1,5
      DO 8 J=1,6
      p(i,j,1)=po(i,j)
      p(i,j,2)=ph(i,j)
      p(i,j,3)=pn(i,j)
      p(i,j,4)=phe(i,j)
c               p(i,j,5)=pno(i,j)
c               p(i,j,6)=po2(i,j)
c               p(i,j,7)=pcl(i,j)
8       continue

      s=0.
c       do 5 i=1,7
      do 5 i=1,4
      do 7 j=1,6
      var(j) = p(1,j,i)*cos(z) + p(2,j,i)*cos(f) +
     &                   p(3,j,i)*cos(0.013*(300.-fs)) +
     &                   p(4,j,i)*cos(0.52*(t-6.)) + p(5,j,i)
7         continue
        cm(i)  = var(1)
        hm(i)  = var(2)
        all(i) = var(3)
        betl(i)= var(4)
        alh(i) = var(5)
        beth(i)= var(6)
        hx=h-hm(i)
        if(hx) 1,2,3
1               arg = hx * (hx * all(i) + betl(i))
         cn(i) = 0.
      if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
      goto 4
2               cn(i) = cm(i)
      goto 4
3               arg = hx * (hx * alh(i) + beth(i))
      cn(i) = 0.
      if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
4         continue
        if(cn(i).LT.0.005*cm(i)) cn(i)=0.
        if(cn(i).GT.cm(i)) cn(i)=cm(i)
        s=s+cn(i)
5       continue
c       do 6 i=1,7
      do 6 i=1,4
6               cn(i)=cn(i)/s*100.
      return
      end
c
           subroutine tcon(yr,mm,day,idn,rz,ig,rsn,nmonth)
c----------------------------------------------------------------
c input:        yr,mm,day       year(yyyy),month(mm),day(dd)
c               idn             day of year(ddd)
c output:       rz(3)           12-month-smoothed solar sunspot number
c               ig(3)           12-month-smoothed IG index
c               rsn             interpolation parameter
c               nmonth          previous or following month depending
c                               on day
c
c Requires I/O UNIT=12 to read the Rz12 and IG12 indices file IG_RZ.DAT 
c 
c rz(1) & ig(1) contain the indices for the month mm and rz(2) & ig(2)
c for the previous month (if day less than 15) or for the following
c month (otherwise). These indices are for the mid of the month. The
c indices for the given day are obtained by linear interpolation and
c are stored in rz(3) and ig(3).
c
c The indices file IG_RZ.DAT is structured as follows (values are 
c separated by comma): 
c   day, month, year of the last update of this file,
c   a blank line
c   start month, start year, end month, end year,
c   a blank line
c   the IG index for December of start year - 1 (needed for interpolation)
c   the 12 IG indices (13-months running mean) for start year, 
c   the 12 IG indices for the second year 
c       .. and so on until the end year,
c   the IG index for January of end year + 1 (needed for interpolation)
c   a blank line
c   the Rz index for December of start year - 1 (needed for interpolation)
c   the 12 Rz indices (13-months running mean) for the start year,
c   the 12 Rz indices for the second year 
c       .. and so on until the end year.
c   the Rz index for January of end year + 1 (needed for interpolation)
c 
c A negative Rz index means that the given index is the 13-months-
C running mean of the solar radio flux (F10.7). The close correlation 
C between (Rz)12 and (F10.7)12 is used to compute the (Rz)12 indices.
c
c An IG index of -111 indicates that no IG values are available for the
c time period. In this case a correlation function between (IG)12 and 
C (Rz)12 is used to obtain (IG)12.
c
c The computation of the 13-month-running mean for month M requires the
c indices for the six months preceeding M and the six months following 
C M (month: M-6, ..., M+6). To calculate the current running mean one 
C therefore requires predictions of the indix for the next six months. 
C Starting from six months before the UPDATE DATE (listed at the top of 
c the file) and onward the indices are therefore based on indices 
c predictions.
c----------------------------------------------------------------

           integer      yr, mm, day, iflag, iyst, iyend,iymst
           integer      imst,iymend
           real         ionoindx(746),indrz(746)
           real         ig(3),rz(3)
           logical      mess
           
           common /iounit/konsol,mess
      COMMON /path/datapath
      character datapath*200

           save         ionoindx,indrz,iflag,iyst,iymst,iymend,imst

        if(iflag.eq.0) then      
           open(unit=12,FILE=TRIM(ADJUSTL(datapath))//'ig_rz.dat',
     &      status='old')

c-web- special for web version
c            open(unit=12,FILE=TRIM(ADJUSTL(datapath))//
c     *         '/usr/local/etc/httpd/cgi-bin/models/IRI/ig_rz.dat',
c     *         status='old')

c Read the update date, the start date and the end date (mm,yyyy), and
c get number of data points to read.

            read(12,*) iupd,iupm,iupy
            read(12,*) imst,iyst,imend,iyend
            iymst=iyst*100+imst
            iymend=iyend*100+imend

c inum_vals= 12-imst+1+(iyend-iyst-1)*12 +imend + 2
c 1st year \ full years       \last y\ before & after

            inum_vals= 3-imst+(iyend-iyst)*12 +imend

c read all the IG12 (ionoindx) and Rz12 (indrz) values

            read(12,*) (ionoindx(i),i=1,inum_vals)
            read(12,*) (indrz(i),i=1,inum_vals)
            do 1 jj=1,inum_vals
                rrr=indrz(jj)
                if(rrr.lt.0.0) then
                    covr=abs(rrr)
                    rrr=33.52*sqrt(covr+85.12)-408.99
                    if(rrr.lt.0.0) rrr=0.0
                    indrz(jj)=rrr
                    endif
                if(ionoindx(jj).gt.-90.) goto 1
                    zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
                    if(zi.gt.274.0) zi=274.0
                    ionoindx(jj)=zi
1               continue
            close(unit=12)
            iflag = 1
        endif

        iytmp=yr*100+mm
        if (iytmp .lt. iymst .or. iytmp .gt. iymend) then
               if(mess) write(konsol,8000) iytmp,iymst,
     &                                            iymend
 8000          format(1x,I10,'** OUT OF RANGE **'/,5x,
     &  'The file IG_RZ.DAT which contains the indices Rz12',
     &  ' and IG12'/5x,'currently only covers the time period',
     &  ' (yymm) : ',I6,'-',I6)
               nmonth=-1
               return
               endif

c       num=12-imst+1+(yr-iyst-1)*12+mm+1
        num=2-imst+(yr-iyst)*12+mm

        rz(1)=indrz(num)
        ig(1)=ionoindx(num)
        midm=15
        if(mm.eq.2) midm=14
        call MODA(0,yr,mm,midm,idd1,nrdaym)
        if(day.lt.midm) goto 1926
c
c day is at or after mid of month
c
                imm2=mm+1
                if(imm2.gt.12) then
                        imm2=1
                        iyy2=yr+1
                        idd2=380            
c               if((yr/4*4.eq.yr).and.(yr/100*100.ne.yr)) idd2=381
                        if(yr/4*4.eq.yr) idd2=381
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num+1)
                ig(2)=ionoindx(num+1)
                rsn=(idn-idd1)*1./(idd2-idd1)                
                rz(3)=rz(1)+(rz(2)-rz(1))*rsn
                ig(3)=ig(1)+(ig(2)-ig(1))*rsn
                goto 1927
1926            imm2=mm-1
                if(imm2.lt.1) then
                        imm2=12
                        idd2=-16
                        iyy2=yr-1
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num-1)
                ig(2)=ionoindx(num-1)
                rsn=(idn-idd2)*1./(idd1-idd2)
                rz(3)=rz(2)+(rz(1)-rz(2))*rsn
                ig(3)=ig(2)+(ig(1)-ig(2))*rsn

1927    nmonth=imm2
            return
            end
C
C
      subroutine tcongec(yr,mm,day,idn,rz,ig,rsn,nmonth)
c----------------------------------------------------------------
c input:        yr,mm,day       year(yyyy),month(mm),day(dd)
c               idn             day of year(ddd)
c output:       rz(3)           12-month-smoothed solar sunspot number
c               ig(3)           12-month-smoothed GEC index
c               rsn             interpolation parameter
c               nmonth          previous or following month depending
c                               on day
c
c rz(1), ig(1) contain the indices for the month mm and rz(2), ig(2)
c for the previous month (if day less than 15) or for the following
c month (otherwise). These indices are for the mid of the month. The
c indices for the given day are obtained by linear interpolation and
c are stored in rz(3) and ig(3).
c
c the indices are obtained from the indices file ig_rz.dat that is 
c read in subroutine initialize and stored in COMMON/indices/
c----------------------------------------------------------------
      integer      yr, mm, day, iflag, iyst, iyend,iymst
      integer      imst,iymend
      real         ionoindx(890),indrz(890)
      real         ig(3),rz(3)
      save         ionoindx,indrz,iflag,iyst,iymst,
     &         iymend,imst
      COMMON /path/datapath
      character datapath*200
c
C NEW TLG Aug. 2013 This Procedure is using Global Electron Content, 
C                   GEC12c index instead of ionospheric IG12 index
C
c Rz12 and GEC are determined from the file GEC_RZ.DAT which has the
c following structure: 
c day, month, year of the last update of this file,
c start month, start year, end month, end year,
c the 12 GEC indices (13-months running mean) for the first year, 
c the 12 GEC indices for the second year and so on until the end year,
c the 12 Rz indices (13-months running mean) for the first year,
c the 12 Rz indices for the second year and so on until the end year.
c The inteporlation procedure also requires the IG and Rz values for
c the month preceeding the start month and the IG and Rz values for the
c month following the end month. These values are also included in 
c GEC_RZ.
c 
c A negative Rz index means that the given index is the 13-months-
C running mean of the solar radio flux (F10.7). The close correlation 
C between (Rz)12 and (F10.7)12 is used to derive the (Rz)12 indices.
c
c An GEC index of -111 indicates that no GEC values are available for the
c time period. In this case a correlation function between (GEC)12 and 
C (Rz)12 is used to obtain (GEC)12.
c
c The computation of the 13-month-running mean for month M requires the
c indices for the six months preceeding M and the six months following 
C M (month: M-6, ..., M+6). To calculate the current running mean one 
C therefore requires predictions of the indix for the next six months. 
C Starting from six months before the UPDATE DATE (listed at the top of 
c the file) and onward the indices are therefore based on indices 
c predictions.
      if(iflag.eq.0) then
         open(unit=12,FILE=TRIM(ADJUSTL(datapath))//'gec_rz.dat',
     &    ERR=10,status='old')
c-web- special for web version
c          open(unit=12,FILE=TRIM(ADJUSTL(datapath))//
c     *'/usr/local/etc/httpd/cgi-bin/models/IRI/gec_rz.dat',
c     *status='old')
      GOTO 20
   10 PAUSE ' CANNOT FIND THE FILE "GEC_RZ.DAT". EXECUTION TEMINATED.'
      STOP
   20 CONTINUE  
c Read the update date, the start date and the end date (mm,yyyy), and
c get number of data points to read.
          read(12,*) iupd,iupm,iupy
          read(12,*) imst,iyst, imend, iyend
          iymst=iyst*100+imst
          iymend=iyend*100+imend
c inum_vals= 12-imst+1+(iyend-iyst-1)*12 +imend + 2
c            1st year \ full years       \last y\ before & after
          inum_vals= 3-imst+(iyend-iyst)*12 +imend
c Read all the ionoindx and indrz values
          read(12,*) (ionoindx(i),i=1,inum_vals)
          read(12,*,err=15,end=15) (indrz(i),i=1,inum_vals)
       goto 16
  15  write(*,*) 'End-of-file i= ',i,indrz(i)
C      pause ' '
      stop
   16       do 1 jj=1,inum_vals
                rrr=indrz(jj)
                if(rrr.lt.0.0) then
                        covr=abs(rrr)
                        rrr=33.52*sqrt(covr+85.12)-408.99
                        if(rrr.lt.0.0) rrr=0.0
                        indrz(jj)=rrr
                        endif
                if(ionoindx(jj).gt.-90.) goto 1
                  zi=((0.0194*rrr+1.0906)-1.0)*50.0       
                  if(zi.gt.274.0) zi=274.0
                  ionoindx(jj)=zi
1               continue
          close(unit=12)
          iflag = 1
        endif

        iytmp=yr*100+mm
        if (iytmp .lt. iymst .or. iytmp .gt. iymend) then
               if(konsol.gt.1) write(konsol,8000) iytmp,iymst,
     &                                            iymend
 8000          format(1x,I10,'** OUT OF RANGE **'/,5x,
     &  'The file GEC_RZ.DAT which contains the indices Rz12',
     &  ' and GEC12'/5x,'currently only covers the time period',
     &  ' (yymm) : ',I6,'-',I6)
               nmonth=-1
               return
               endif

c       num=12-imst+1+(yr-iyst-1)*12+mm+1
        num=2-imst+(yr-iyst)*12+mm

        rz(1)=indrz(num)
        ig(1)=ionoindx(num)
        midm=15
        if(mm.eq.2) midm=14
        call MODA(0,yr,mm,midm,idd1,nrdaym)
        if(day.lt.midm) goto 1926
                imm2=mm+1
                if(imm2.gt.12) then
                        imm2=1
                        iyy2=yr+1
                        idd2=380
c               if((yr/4*4.eq.yr).and.(yr/100*100.ne.yr)) idd2=381
                        if(yr/4*4.eq.yr) idd2=381

                        if(yr/4*4.eq.yr) idd2=381
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num+1)
                ig(2)=ionoindx(num+1)
                rsn=(idn-idd1)*1./(idd2-idd1)
                rz(3)=rz(1)+(rz(2)-rz(1))*rsn
                ig(3)=ig(1)+(ig(2)-ig(1))*rsn
                goto 1927
1926            imm2=mm-1
                if(imm2.lt.1) then
                        imm2=12
                        idd2=-16
                        iyy2=yr-1
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num-1)
                ig(2)=ionoindx(num-1)
                rsn=(idn-idd2)*1./(idd1-idd2)
                rz(3)=rz(2)+(rz(1)-rz(2))*rsn
                ig(3)=ig(2)+(ig(1)-ig(2))*rsn

1927    nmonth=imm2
            return
            end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine tconind(yr,mm,day,idn,rz,sf,rsn,nmonth,nind)
c----------------------------------------------------------------
C  TL Gulyaeva........................................................ May, 2017
C     nind =1 <DEFAULT> ssn1_12.dat  =>  rz(3); F10.7  =>  sf(3) 
C     nind =2 ssn2_12.dat  SSN2 converted to SSN1  =>  rz(3); F10.7  => sf(3)
C     nind =3 f107_12.dat  F10.7    =>  sf(3) converted to SSN1  =>  rz(3) 
C     nind =4 geci_12.dat  GEC => rz(3); F10.7  => sf(3)
C     nind =5 teci_12.dat  TEC => rz(3); F10.7  => sf(3)
C     nind =6 igin_12.dat  IG  => rz(3); F10.7  => sf(3)
C     nind =7 mgii_12.dat  MgII => rz(3); F10.7  => sf(3)
C     nind =8 lyma_12.dat  Lyman-alpha => rz(3); F10.7  => sf(3)
C Ref. T.L.Gulyaeva et al., 2017. TEC proxy index of solar activity for the
C      International Reference Ionosphere IRI and its extension to 
C      Plasmasphere IRI-PLAS model.Int. J. Sci. Eng. Applied Sci.,3,5,144-150,
C      http://ijseas.com/index.php/issue-archive-2/volume3/issue-5/
C
c input:        yr,mm,day       year(yyyy),month(mm),day(dd)
c               idn             day of year(ddd)
c               nind            type of solar proxy index
c output:       rz(3)           12-month-smoothed solar sunspot number or proxy index specified by user
c               sf(3)           12-month-smoothed solar radio flux F10.7
c               rsn             interpolation parameter
c               nmonth          previous or following month depending
c                               on day
c
c rz(1), sf(1) contain the indices for the month mm and rz(2), sf(2)
c for the previous month (if day less than 15) or for the following
c month (otherwise). These indices are for the mid of the month. The
c indices for the given day are obtained by linear interpolation and
c are stored in rz(3) and sf(3).
c
c the indices are obtained from the indices files according to nind  
c 
c----------------------------------------------------------------
      integer      yr, mm, day
      integer      nind
      real         sf(3),rz(3)
      character*4  xxpr,yypr
      COMMON /BXY/xxpr,yypr
c
C
c rz(3) is read from xxpr_12.dat file (xxpr is defined according to nind) 
c and the sf(3) solar radio flux index is determined from yypr_12.dat file 
c which have the following structure: 
c day, month, year of the last update of this file,
c start month, start year, end month, end year,
c the 12 indices (13-months running mean) for the first year, 
c the 12 indices for the second year and so on until the end year,
c The inteporlation procedure also requires the proxy index values for
c the month preceeding the start month and the proxy index values for the
c month following the end month. These values are also included in 
c xxpr_12.dat and yypr_12.dat file.
c 
c Input of nind=1 (default) means that the given index is the 13-months-running mean of 
c SSN1 (former Rz12) and the solar radio flux (F10.7). The close correlation between Rz12 
c and other solar proxies is used to derive the rz(3) and sf(3) indices.
c
c
c The computation of the 13-month-running mean for month M requires the
c indices for the six months preceeding M and the six months following 
C M (month: M-6, ..., M+6). To calculate the current running mean one 
C therefore requires predictions of the index for the next six months. 
C Starting from six months before the UPDATE DATE (listed at the top of 
c the file) and onward the indices are therefore based on indices 
c predictions.
C
C Check nind = 1,...,8::
         if ((nind.eq.0).or.(nind.gt.8)) then 
      write(*,*) ' WRONG OPTION FOR SOLAR PROXY JIND = ',nind,
     +   ' EXECUTION TEMINATED.'
      PAUSE  ' '
      STOP
      endif
C Start to determine prec. month or next month:
      midm=15
        if(mm.eq.2) midm=14
        call MODA(0,yr,mm,midm,idd1,nrdaym)
        if(day.lt.midm) goto  926
                imm2=mm+1
                if(imm2.gt.12) then
                        imm2=1
                        iyy2=yr+1
                        idd2=380
                        if(yr/4*4.eq.yr) idd2=381
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                goto 927
  926            imm2=mm-1
                if(imm2.lt.1) then
                      imm2=12
                        idd2=-16
                        iyy2=yr-1
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
C
  927 continue
C
C
        select case (nind)
C 
      case (1)   
C        
    1 xxpr='ssn1'      
      call subproxy(xxpr,yr,mm,iyy2,imm2,rz) 
      yypr='f107'
      call subproxy(yypr,yr,mm,iyy2,imm2,sf)  
      goto 33
C 
      case (2)  
C        
    2 xxpr='ssn2'      
      call subproxy(xxpr,yr,mm,iyy2,imm2,rz) 
      rz(1)=0.7*rz(1)
      rz(2)=0.7*rz(2)
      yypr='f107'
      call subproxy(yypr,yr,mm,iyy2,imm2,sf)  
      goto 33
C        
      case (3)   
C 
    3 yypr='f107'      
      call subproxy(yypr,yr,mm,iyy2,imm2,sf) 
      rz(1)=1.076*sf(1)-65.7817 
      rz(2)=1.076*sf(2)-65.7817  
       if (rz(1).lt.0.) rz(1)=0.
       if (rz(2).lt.0.) rz(2)=0.
      goto 33
C 
      case (4)  
C 
    4 xxpr='geci'      
      call subproxy(xxpr,yr,mm,iyy2,imm2,rz) 
      yypr='f107'
      call subproxy(yypr,yr,mm,iyy2,imm2,sf) 
      goto 33
C 
      case (5)   
C 
    5 xxpr='teci'      
      call subproxy(xxpr,yr,mm,iyy2,imm2,rz) 
      yypr='f107'
      call subproxy(yypr,yr,mm,iyy2,imm2,sf) 
      goto 33
C 
      case (6)  
C 
    6 xxpr='igin'      
      call subproxy(xxpr,yr,mm,iyy2,imm2,rz) 
      yypr='f107'
      call subproxy(yypr,yr,mm,iyy2,imm2,sf)  
      goto 33
C 
      case (7)  
C 
    7 xxpr='mgii'      
      call subproxy(xxpr,yr,mm,iyy2,imm2,rz) 
      yypr='f107'
      call subproxy(yypr,yr,mm,iyy2,imm2,sf)  
      goto 33
C 
      case (8)   
C 
    8 xxpr='lyma'      
      call subproxy(xxpr,yr,mm,iyy2,imm2,rz)  
      yypr='f107'
      call subproxy(yypr,yr,mm,iyy2,imm2,sf) 
      goto 33
C+++++++
      end select 
   33 continue
        if(day.lt.midm) goto 34
      rsn=(idn-idd1)*1./(idd2-idd1)
      rz(3)=rz(1)+(rz(2)-rz(1))*rsn
      sf(3)=sf(1)+(sf(2)-sf(1))*rsn
      goto 35
C
   34 rsn=(idn-idd2)*1./(idd1-idd2)
      rz(3)=rz(2)+(rz(1)-rz(2))*rsn
      sf(3)=sf(2)+(sf(1)-sf(2))*rsn
   35 nmonth=imm2
      return
      end
C
C----------------------------------------------------------------------------
C
      subroutine subproxy(xxxx,kyr1,kmn1,kyr2,kmn2,zz)
C
C  TL Gulyaeva....................................................May, 2017
C Subroutine to read 12-monthly smoothed solar proxy data
C for given year=kyr1 and month = kmn1 
C kyr2,kmn2 denoted either preceding month or the following month
C zz(1) results for the given month kyr1,kmn1, and zz(2) for kyr2,kmn2  
C
      CHARACTER*4 xxxx
      CHARACTER*12 infile
      real zz(3),dat(12)
      COMMON /path/datapath
      character datapath*200
C
      infile='xxxx_12.dat'
      infile(1:4)=xxxx
      zz(1)=-1.0
      zz(2)=-1.0
      zz(3)=0.0
C
      open(unit=12,FILE=TRIM(ADJUSTL(datapath))//infile,
     & ERR=10,status='old')
      GOTO 20
   10 write(*,*) ' CANNOT FIND THE FILE ',infile,' EXECUTION TEMINATED.'
      PAUSE  ' '
      STOP
   20 CONTINUE  
c Read the update date, the start date and the end date (mm,yyyy)
c 
          read(12,*) iupd,iupm,iupy
          read(12,*) imst,iyst, imend, iyend
c
            if (kyr2.lt.iyst) then
         kyr2=iyst
         kmn2=1
         endif
         if (kyr2.gt.iyend) then
         kyr2=iyend
         kmn2=12
         endif
c
      if ((kyr1.lt.iyst).or.(kyr1.gt.iyend)) then
      write (*,*) ' REQUIED YEAR OUTSIDE OF DATA PERIOD IN THE FILE '
     +,infile,' EXECUTION TEMINATED.'
      pause ' '
      STOP
      endif
c
c Read year and monthly input values
   11 read(12,*,err=15,end=15) kyyyy, (dat(i),i=1,12)
      IF (kyr1.eq.kyr2) THEN               
            if (kyyyy.eq.kyr1) then
            zz(1)=dat(kmn1)
            zz(2)=dat(kmn2)
            goto 16
                     else
               goto 11
            endif
         ENDIF                        
            IF (kyr1.lt.kyr2) THEN            !+
C kyr1 < kyr2
           IF (kyr1.eq.kyyyy) THEN          !
            zz(1)=dat(kmn1)
             if (kyr2.gt.iyend) then        !!
                zz(2)=zz(1)
               goto 16
                 else                       !!
              read(12,*,err=15,end=15) kyyyy,zz(2)
                goto 16
                endif                       !!
            ELSE                      !
              goto 11
          ENDIF                             !
                           ELSE           !+
C kyr2 < kyr1
             IF (kyr2.ge.iyst) THEN         !
C  kyr2 > iyst             
          if (kyr2.eq.kyyyy) then           !!
            zz(2)=dat(kmn2)
                read(12,*,err=15,end=15) kyyyy,zz(1)
                goto 16
            else                      !!
              goto 11
          endif                             !!
                              ELSE           !
C kyr2 < iyst
              if (kyr1.eq.kyyyy) then        !!!
         zz(1)=dat(kmn1)
         zz(2)=zz(1)
                        else           !!!
              goto 11
                          endif           !!!                     
         ENDIF                         !!
               ENDIF                         !+
C
  15  write(*,*) 'End-of-file ',infile
      pause ' '
      stop
C  
   16     close(unit=12)
            return
            end
C -----------------------------------------------------------
      subroutine blet1(ilet,alet1)
c TL Gulyaeva........................................May, 2017
C
      CHARACTER*1 IN(0:9)
      DATA IN/'0','1','2','3','4','5','6','7','8','9'/
      CHARACTER*1 alet1
      alet1='0'
    5 do i=0,9
      if (ilet.eq.i) then
      alet1=IN(i)
      endif
      enddo
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine LSTID(FI,ICEZ,R,AE,TM,SAX,SUX,TS70,DF0F2,DHF2)
C*****************************************************************

C   COMPUTER PROGRAM FOR UPDATING FOF2 AND HMF2 FOR EFFECTS OF
C   THE LARGE SCALE SUBSTORM.
C
C   P.V.Kishcha, V.M.Shashunkina, E.E.Goncharova, Modelling of the
C   ionospheric effects of isolated and consecutive substorms on
C   the basis of routine magnetic data, Geomagn. and Aeronomy v.32,
C   N.3, 172-175, 1992.
C
C   P.V.Kishcha et al. Updating the IRI ionospheric model for
C   effects of substorms, Adv. Space Res.(in press) 1992.
C
C   Address: Dr. Pavel V. Kishcha,
C            Institute of Terrestrial Magnetism,Ionosphere and Radio
C            Wave Propagation, Russian Academy of Sciences,
C            142092, Troitsk, Moscow Region, Russia
C
C***       INPUT PARAMETERS:
C       FI ------ GEOMAGNETIC LATITUDE,
C       ICEZ ---- INDEX OF SEASON(1-WINTER AND EQUINOX,2-SUMMER),
C       R ------- ZURICH SUNSPOT NUMBER,
C       AE ------ MAXIMUM AE-INDEX REACHED DURING SUBSTORM,
C       TM ------ LOCAL TIME,
C       SAX,SUX - TIME OF SUNSET AND SUNRISE,
C       TS70 ---- ONSET TIME (LT) OF SUBSTORMS ONSET
C                        STARTING ON FI=70 DEGR.
C***      OUTPUT PARAMETERS:
C       DF0F2,DHF2- CORRECTIONS TO foF2 AND hmF2 FROM IRI OR
C                   OBSERVATIONAL MEDIAN  OF THOSE VALUES.
C*****************************************************************


      INTEGER ICEZ
      REAL A(7,2,3,2),B(7,2,3,2),C(7,2,3,2),D(7,2,3,2),A1(7,2,2),
     *B1(7,2,2),Y1(84),Y2(84),Y3(84),Y4(84),Y5(28),Y6(28)
      DATA Y1/
     *150.,250.,207.8,140.7,158.3,87.2,158.,
     *150.,250.,207.8,140.7,158.3,87.2,158.,
     *115.,115.,183.5,144.2,161.4,151.9,272.4,
     *115.,115.,183.5,144.2,161.4,151.9,272.4,
     *64.,320.,170.6,122.3,139.,79.6,180.6,
     *64.,320.,170.6,122.3,139.,79.6,180.6,
     *72.,84.,381.9,20.1,75.1,151.2,349.5,
     *120.,252.,311.2,241.,187.4,230.1,168.7,
     *245.,220.,294.7,181.2,135.5,237.7,322.,
     *170.,110.,150.2,136.3,137.4,177.,114.,
     *170.,314.,337.8,155.5,157.4,196.7,161.8,
     *100.,177.,159.8,165.6,137.5,132.2,94.3/
      DATA Y2/
     *2.5,2.0,1.57,2.02,2.12,1.46,2.46,
     *2.5,2.0,1.57,2.02,2.12,1.46,2.46,
     *2.3,1.6,1.68,1.65,2.09,2.25,2.82,
     *2.3,1.6,1.68,1.65,2.09,2.25,2.82,
     *0.8,2.0,1.41,1.57,1.51,1.46,2.2,
     *0.8,2.0,1.41,1.57,1.51,1.46,2.2,
     *3.7,1.8,3.21,3.31,2.61,2.82,2.34,
     *2.8,3.2,3.32,3.33,2.96,3.43,2.44,
     *3.5,2.8,2.37,2.79,2.26,3.4,2.28,
     *3.9,2.,2.22,1.98,2.33,3.07,1.56,
     *3.7,3.,3.3,2.99,3.57,2.98,3.02,
     *2.6,2.8,1.66,2.04,1.91,1.49,0.43/
      DATA Y3/
     *-1.8,-1.9,-1.42,-1.51,-1.53,-1.05,-1.66,
     *-1.8,-1.9,-1.42,-1.51,-1.53,-1.05,-1.66,
     *-1.5,-1.3,-1.46,-1.39,-1.53,-1.59,-1.9,
     *-1.5,-1.3,-1.46,-1.39,-1.53,-1.59,-1.9,
     *-0.7,-2.,-1.41,-1.09,-1.22,-0.84,-1.32,
     *-0.7,-2.,-1.41,-1.09,-1.22,-0.84,-1.32,
     *-1.7,-1.0,-2.08,-1.80,-1.35,-1.55,-1.79,
     *-1.5,-2.,-2.08,-2.16,-1.86,-2.19,-1.70,
     *-2.2,-1.7,-1.57,-1.62,-1.19,-1.89,-1.47,
     *-1.9,-1.5,-1.26,-1.23,-1.52,-1.89,-1.02,
     *-1.7,-1.7,-1.76,-1.43,-1.66,-1.54,-1.24,
     *-1.1,-1.5,-1.09,-1.23,-1.11,-1.14,-0.4/
      DATA Y4/
     *-2.,-5.,-5.,0.,0.,0.,2.,
     *-2.,-5.,-5.,0.,0.,0.,2.,
     *-5.,-5.,6.,0.,1.,5.,2.,
     *-5.,-5.,6.,0.,1.,5.,2.,
     *0.,-7.,-3.,-6.,2.,2.,3.,
     *0.,-7.,-3.,-6.,2.,2.,3.,
     *-5.,-1.,-11.,-6.,0.,-5.,-6.,
     *-5.,-10.,1.,4.,-6.,-2.,1.,
     *2.,-13.,-10.,0.,-8.,10.,-16.,
     *0.,-3.,-7.,-2.,-2.,4.,2.,
     *-11.,-12.,-13.,0.,0.,7.,0.,
     *-8.,6.,-1.,-5.,-7.,4.,-4./
      DATA Y5/
     *0.,0.,-0.1,-0.19,-0.19,-0.25,-0.06,
     *0.,0.,-0.31,-0.28,-0.27,-0.06,0.02,
     *0.,0.,0.18,-0.07,-0.2,-0.1,0.3,
     *0.,0.,-0.24,-0.5,-0.4,-0.27,-0.48/
      DATA Y6/
     *0.,0.,-0.00035,-0.00028,-0.00033,-0.00023,-0.0007,
     *0.,0.,-0.0003,-0.00025,-0.0003,-0.0006,-0.00073,
     *0.,0.,-0.0011,-0.0006,-0.0003,-0.0005,-0.0015,
     *0.,0.,-0.0008,-0.003,-0.0002,-0.0005,-0.0003/

      INN=0
      IF(TS70.GT.12. .AND. TM.LT.SAX)INN=1
      IF(FI.LT.0.)FI=ABS(FI)

      N=0
      DO 001 M=1,2
      DO 001 K=1,3
      DO 001 J=1,2
      DO 001 I=1,7
      N=N+1
      A(I,J,K,M)=Y1(N)
      B(I,J,K,M)=Y2(N)
      C(I,J,K,M)=Y3(N)
 001    D(I,J,K,M)=Y4(N)
      N1=0
      DO 002 M=1,2
      DO 002 J=1,2
      DO 002 I=1,7
      N1=N1+1
      A1(I,J,M)=Y5(N1)
 002    B1(I,J,M)=Y6(N1)
      IF(FI.GT.65..OR.AE.LT.500.)THEN
      WRITE(*,*)'LSTID are for AE>500. and ABS(FI)<65.'
      GOTO 004
      ENDIF
      TS=TS70+(-1.5571*FI+109.)/60.
      IF(TS.LT.SUX.AND.TS.GT.SAX)THEN
      WRITE(*,*)' LSTID are only at night'
      GOTO 004
      ENDIF
      IF(INN.EQ.1)TM=TM+24.

      IF(TS.GE.TM.OR.TS.LT.TM-5.)THEN
C        WRITE(*,*)'LSTID are onli if  TM-5.<TS<TM ;Here TS=',TS,'TM=',TM
      GOTO 004
      ENDIF
      DO 007 I=1,7
      IF(FI.GE.-5.+10.*(I-1) .AND. FI.LT.5.+10.*(I-1))GOTO 008
 007    CONTINUE
 008    J=ICEZ
      IF(AE.GE.500. .AND. AE.LE.755.)K=1
      IF(AE.GT.755. .AND. AE.LT.1000.)K=2
      IF(AE.GE.1000.)K=3
      M=-1
      IF(R.LE.20.)M=1
      IF(R.GE.120.)M=2
      T=TM-TS
      IF(M.LT.0)GOTO 003
C        WRITE(*,*)'A1=',A1(I,J,M),' B1=',B1(I,J,M)
C        WRITE(*,*)'A=',A(I,J,K,M),' B=',B(I,J,K,M),' C=',C(I,J,K,M),
C     *'D=',D(I,J,K,M)
      DF0F2=A1(I,J,M)+B1(I,J,M)*AE
      DHF2=A(I,J,K,M)*(T**B(I,J,K,M))*EXP(C(I,J,K,M)*T)+D(I,J,K,M)
      GOTO 005
 003    DF1=A1(I,J,1)+B1(I,J,1)*AE
      DF2=A1(I,J,2)+B1(I,J,2)*AE
      DF0F2=DF1+(DF2-DF1)*(R-20.)/100.
      DH1=A(I,J,K,1)*(T**B(I,J,K,1))*EXP(C(I,J,K,1)*T)+D(I,J,K,1)
      DH2=A(I,J,K,2)*(T**B(I,J,K,2))*EXP(C(I,J,K,2)*T)+D(I,J,K,2)
      DHF2=DH1+(DH2-DH1)*(R-20.)/100.
      GOTO 005
 004    DHF2=0.
      DF0F2=0.
 005    CONTINUE
      IF(INN.EQ.1)TM=TM-24.
      RETURN
      END
c
c  to get ionospheric correction for frequency f/Hz use the following
C  formula   ioncorr/m = 40.3 * TEC / (f*f)
c
c
c
      real function ioncorr(tec,f)
c-------------------------------------------------------------------
c computes ionospheric correction IONCORR (in m) for given vertical
c ionospheric electron content TEC (in m-2) and frequency f (in Hz)
c-------------------------------------------------------------------
      ioncorr = 40.3 * tec / (f*f)
      return
      end
c
c
c
        subroutine iri_tec (hstart,hend,istep,tectot,tectop,tecbot)
c-----------------------------------------------------------------------        
C subroutine to compute the total ionospheric content
C INPUT:      
C   hstart  altitude (in km) where integration should start
C   hend    altitude (in km) where integration should end
C   istep   =0 [fast, but higher uncertainty <5%]
C           =1 [standard, recommended]
C           =2 [stepsize of 1 km; best TEC, longest CPU time]
C OUTPUT:
C   tectot  total ionospheric content in tec-units (10^16 m^-2)
C   tectop  topside content (in %)
C   tecbot  bottomside content (in %)
C
C The different stepsizes for the numerical integration are 
c defined as follows (h1=100km, h2=hmF2-10km, h3=hmF2+10km, 
c h4=hmF2+150km, h5=hmF2+250km):
C       istep   h1-h2   h2-h3   h3-h4   h4-h5   h5-hend
C       0       2.0km   1.0km   2.5km   exponential approximation
C       1       2.0km   1.0km   2.5km   10.0km  30.0km
C       2       1.0km   0.5km   1.0km   1.0km   1.0km   
C
c-----------------------------------------------------------------------        

        logical         expo
        dimension       step(5),hr(6)
        common  /block1/hmf2,xnmf2,hmf1,f1reg
        logical     f1reg

ctest   
        save

        expo = .false.
        numstep = 5
        xnorm = xnmf2/1000.

        hr(1) = 100.
        hr(2) = hmf2-10.
        hr(3) = hmf2+10.
        hr(4) = hmf2+150.
        hr(5) = hmf2+250.
        hr(6) = hend
        do 2918 i=2,6 
2918            if (hr(i).gt.hend) hr(i)=hend

        if (istep.eq.0) then 
                step(1)=2.0
                step(2)=1.0
                step(3)=2.5
                step(4)=5.
                if (hend.gt.hr(5)) expo=.true.
                endif

        if (istep.eq.1) then
                step(1)=2.0
                step(2)=1.0
                step(3)=2.5
                step(4)=10.0
                step(5)=30.0
                endif

        if (istep.eq.2) then
                step(1)=1.0
                step(2)=0.5
                step(3)=1.0
                step(4)=1.0
                step(5)=1.0
                endif

        sumtop = 0.0
        sumbot = 0.0
C
C find the starting point for the integration
C

        i=0
        ia=1
3       i=i+1
        h=hr(i)
        if(hstart.gt.h) then
                hr(i)=hstart
                ia=i
                goto 3
                endif
C
C start the numerical integration
C
        i=ia
        h=hr(i)
        hu=hr(i+1)
        delx = step(i)
1       h = h + delx
        hh = h
        if (h.ge.hu) then
                delx = hu - h + delx
                hx = hu - delx/2.
cc                YNE = XE_1(hx)
                YNE = XE(hx)
                if((hx.gt.hmf2).and.(yne.gt.xnmf2)) yne=xnmf2
                yyy = yne * delx / xnorm
                i=i+1
            if(i.lt.6) then
                      h = hr(i)
                      hu = hr(i+1)
                            delx = step(i)
                  endif
        else
              hx = h - delx/2.
cc                YNE = XE_1(hx)
                YNE = XE(hx)
                if((hx.gt.hmf2).and.(yne.gt.xnmf2)) yne=xnmf2
                yyy = yne * delx / xnorm
        endif
        if (hx.le.hmf2) then
                sumbot = sumbot + yyy
        else
                sumtop = sumtop + yyy
        endif
        if (expo.and.(hh.ge.hr(4))) goto 5
        if (hh.lt.hend.and.i.lt.6) goto 1

        zzz = sumtop + sumbot
        tectop = sumtop / zzz * 100.
        tecbot = sumbot / zzz * 100.
        tectot = zzz * xnmf2    
        return

5       num_step = 3
        hei_top = hr(4)
        hei_end = hend
        top_end = hei_end - hei_top
        del_hei = top_end / num_step
cc        xntop = xe_1(hei_end)/xnmf2
        xntop = xe(hei_end)/xnmf2

        if(xntop.gt.0.9999) then
                ss_t = top_end  
                goto 2345
                endif

        hei_2 = hei_top
        hei_3 = hei_2 + del_hei
        hei_4 = hei_3 + del_hei
        hei_5 = hei_end

        hss = top_end / 4.
C       hss = 360.
        xkk = exp ( - top_end / hss ) - 1.
        x_2 = hei_2
        x_3 =hei_top-hss*alog(xkk*(hei_3 - hei_top)/top_end + 1.) 
        x_4 =hei_top-hss*alog(xkk*(hei_4 - hei_top)/top_end + 1.)
        x_5 = hei_end

cc        ed_2 = xe_1(x_2)/xnmf2
        ed_2 = xe(x_2)/xnmf2
          if(ed_2.gt.1.) ed_2=1.
cc        ed_3 = xe_1(x_3)/xnmf2
        ed_3 = xe(x_3)/xnmf2
          if(ed_3.gt.1.) ed_3=1.
cc        ed_4 = xe_1(x_4)/xnmf2
        ed_4 = xe(x_4)/xnmf2
          if(ed_4.gt.1.) ed_4=1.
        ed_5 = xntop
        if(ed_3.eq.ed_2) then
         ss_2 = ed_3 * (x_3 - x_2)
        else
         ss_2=( ed_3 - ed_2 ) * ( x_3 - x_2 ) / alog ( ed_3 / ed_2 )
        endif
        if(ed_4.eq.ed_3) then
         ss_3 = ed_4 * (x_4 - x_3)
        else
         ss_3=( ed_4 - ed_3 ) * ( x_4 - x_3 ) / alog ( ed_4 / ed_3 )
        endif
        if(ed_5.eq.ed_4) then
         ss_4 = ed_5 * (x_5 - x_4)
        else
         ss_4=( ed_5 - ed_4 ) * ( x_5 - x_4 ) / alog ( ed_5 / ed_4 )
        endif

        ss_t = ss_2 + ss_3 + ss_4 

2345    sumtop = sumtop + ss_t * 1000.
        
        zzz = sumtop + sumbot
        tectop = sumtop / zzz * 100.
        tecbot = sumbot / zzz * 100.
        tectot = zzz * xnmf2

      RETURN
      END

C           SMI TECPL INTEGRATION
C=====================================
      SUBROUTINE dinteg(ah,hmin,dnemin,tiec)
C==== SMI TECpl INTEGRATION PROCEDURE
C      dimension outf(0:11,0:120)
      integer daynr
      logical F1REG
      common /block1/hmf2,xnmf2,hmf1,F1REG
     &/dem/w,glongi,HOUR,daynr,cn1000,alon,blon,enre,HSC
      dnemin=xe(hmin)
c      k=ifix((ah-hmin)/20.+0.5)
      K=INT((AH-HMIN)/20.+0.5)
      x1=hmin
      y1=dnemin
      st=20.0
      SUM=0.0
      do 1 l=1,k
      x2=x1+st
      y2=xe(x2)
      SUM=SUM+(y2+y1)*st
      x1=x2
      y1=y2
    1 continue
      tiec=sum*500.
      return
      end
c
C========================================================
      SUBROUTINE dintegr(ah,ah0,tiec)
C calculate TEC for given range [
C==== SMI TECpl INTEGRATION PROCEDURE
C      dimension outf(0:11,0:120)
      integer daynr
      logical F1REG
      common /block1/hmf2,xnmf2,hmf1,F1REG
     &/dem/w,glongi,HOUR,daynr,cn1000,alon,blon,enre,HSC
      dnemin=xe(hmin)
c      k=ifix((ah-hmin)/20.+0.5)
      K=INT((AH-AH0)/20.+0.5)
      x1=ah0
      y1=xe(x1)
      st=20.0
      SUM=0.0
      do 1 l=1,k
      x2=x1+st*l
      y2=xe(x2)
      SUM=SUM+(y2+y1)*st
      y1=y2
    1 continue
      tiec=sum*500.
      return
      end
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CIRA86(IDAY,SEC,GLAT,GLONG,STL,F107A,TINF,TLB,SIGMA)
C*******************************************************************
C   Calculates neutral temperature parameters for IRI using the
C   MSIS-86/CIRA 1986 Neutral Thermosphere Model. The subroutines
C   GTS5, GLOBE5 and GLOBL5 developed by A.E. Hedin (2/26/87) were
C   modified for use in IRI --------- D. Bilitza -------- March 1991
C
C   CHANGES:
C     11/09/99 always calculated Legendre; 'if glat' and 'if stl' taken out;
C     11/09/99 use UMR, dumr and humr from COMMON 
C
C     INPUT:
C        IDAY - DAY OF YEAR 
C        SEC - UT(SEC)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C
C     OUTPUT: 
C        TINF - EXOSPHERIC TEMPERATURE (K)
C        TLB - TEMPERATURE AT LOWER BOUNDARY (K)
C        SIGMA - SHAPE PARAMETER FOR TEMPERATURE PROFILE
C
C **********************************************************************
      DIMENSION         PLG(9,4)
     *,indap(13),akp(13)
        common  /const/umr,PI   /const1/hr,dr,xfap,indap,akp,cov,rzs
      DATA      XL/1000./,TLL/1000./

c  data umr/1.74E-2/,hr/0.2618/,dr/1.74e-2
cDR,DR2/1.72142E-2,0.0344284/,
cSR/7.2722E-5/,
c,HR/.2618/
c,DGTR/1.74533E-2/
c       dr = hr * 24. / 365.
c
        dr2 = dr * 2.
        sr = hr / 3600.

C
C CALCULATE LEGENDRE POLYNOMIALS
C
      IF(XL.EQ.GLAT)   GO TO 15
      C = SIN(GLAT*umr)
      S = COS(GLAT*umr)
      C2 = C*C
      C4 = C2*C2
      S2 = S*S
      PLG(2,1) = C
      PLG(3,1) = 0.5*(3.*C2 -1.)
      PLG(4,1) = 0.5*(5.*C*C2-3.*C)
      PLG(5,1) = (35.*C4 - 30.*C2 + 3.)/8.
      PLG(6,1) = (63.*C2*C2*C - 70.*C2*C + 15.*C)/8.
      PLG(2,2) = S
      PLG(3,2) = 3.*C*S
      PLG(4,2) = 1.5*(5.*C2-1.)*S
      PLG(5,2) = 2.5*(7.*C2*C-3.*C)*S
      PLG(6,2) = 1.875*(21.*C4 - 14.*C2 +1.)*S
      PLG(7,2) = (11.*C*PLG(6,2)-6.*PLG(5,2))/5.
      PLG(3,3) = 3.*S2
      PLG(4,3) = 15.*S2*C
      PLG(5,3) = 7.5*(7.*C2 -1.)*S2
      PLG(6,3) = 3.*C*PLG(5,3)-2.*PLG(4,3)
      PLG(4,4) = 15.*S2*S
      PLG(5,4) = 105.*S2*S*C 
      PLG(6,4)=(9.*C*PLG(5,4)-7.*PLG(4,4))/2.
      PLG(7,4)=(11.*C*PLG(6,4)-8.*PLG(5,4))/3.
      XL=GLAT
   15 CONTINUE
      IF(TLL.EQ.STL)   GO TO 16
      STLOC = SIN(HR*STL)
      CTLOC = COS(HR*STL)
      S2TLOC = SIN(2.*HR*STL)
      C2TLOC = COS(2.*HR*STL)
      S3TLOC = SIN(3.*HR*STL)
      C3TLOC = COS(3.*HR*STL)
      TLL = STL
   16 CONTINUE
C
      DFA=F107A-150.
C
C EXOSPHERIC TEMPERATURE
C
C         F10.7 EFFECT
      T1 =  ( 3.11701E-3 - 0.64111E-5 * DFA ) * DFA
        F1 = 1. + 0.426385E-2 * DFA 
        F2 = 1. + 0.511819E-2 * DFA
        F3 = 1. + 0.292246E-2 * DFA
C        TIME INDEPENDENT
      T2 = 0.385528E-1 * PLG(3,1) + 0.303445E-2 * PLG(5,1)
C        SYMMETRICAL ANNUAL AND SEMIANNUAL
        CD14 = COS( DR  * (IDAY+8.45398) )
        CD18 = COS( DR2 * (IDAY-125.818) )
        CD32 = COS( DR  * (IDAY-30.0150) )
        CD39 = COS( DR2 * (IDAY-2.75905) )
      T3 = 0.805486E-2 * CD32 + 0.14237E-1 * CD18
C        ASYMMETRICAL ANNUAL AND SEMIANNUAL
      T5 = F1 * (-0.127371 * PLG(2,1) - 0.302449E-1 * PLG(4,1) ) * CD14
     &       - 0.192645E-1 * PLG(2,1) * CD39
C        DIURNAL
      T71 =  0.123512E-1 * PLG(3,2) * CD14
      T72 = -0.526277E-2 * PLG(3,2) * CD14
      T7 = ( -0.105531 *PLG(2,2) - 0.607134E-2 *PLG(4,2) + T71 ) *CTLOC
     4   + ( -0.115622 *PLG(2,2) + 0.202240E-2 *PLG(4,2) + T72 ) *STLOC
C        SEMIDIURNAL
      T81 = 0.386578E-2 * PLG(4,3) * CD14
      T82 = 0.389146E-2 * PLG(4,3) * CD14
      T8= (-0.516278E-3 *PLG(3,3) - 0.117388E-2 *PLG(5,3) +T81)*C2TLOC
     3   +( 0.990156E-2 *PLG(3,3) - 0.354589E-3 *PLG(5,3) +T82)*S2TLOC
C        TERDIURNAL
      Z1 =  PLG(5,4) * CD14
      Z2 =  PLG(7,4) * CD14
      T14=(0.147284E-2*PLG(4,4)-0.173933E-3*Z1+0.365016E-4*Z2)*S3TLOC
     2   +(0.341345E-3*PLG(4,4)-0.153218E-3*Z1+0.115102E-3*Z2)*C3TLOC 
      T7814 = F2 * ( T7 + T8 + T14 )
C        LONGITUDINAL
      T11= F3 * (( 0.562606E-2 * PLG(3,2) + 0.594053E-2 * PLG(5,2) + 
     $       0.109358E-2 * PLG(7,2) - 0.301801E-2 * PLG(2,2) - 
     $       0.423564E-2 * PLG(4,2) - 0.248289E-2 * PLG(6,2) + 
     $      (0.189689E-2 * PLG(2,2) + 0.415654E-2 * PLG(4,2)) * CD14
     $     ) * COS(umr*GLONG) +
     $     ( -0.11654E-1 * PLG(3,2) - 0.449173E-2 * PLG(5,2) - 
     $       0.353189E-3 * PLG(7,2) + 0.919286E-3 * PLG(2,2) + 
     $       0.216372E-2 * PLG(4,2) + 0.863968E-3 * PLG(6,2) +
     $      (0.118068E-1 * PLG(2,2) + 0.331190E-2 * PLG(4,2)) * CD14
     $     ) * SIN(umr*GLONG) )
C        UT AND MIXED UT,LONGITUDE
      T12 = ( 1. - 0.565411 * PLG(2,1) ) * COS( SR*(SEC-31137.0) ) *
     $ (-0.13341E-1*PLG(2,1)-0.243409E-1*PLG(4,1)-0.135688E-1*PLG(6,1))
     $ +    ( 0.845583E-3 * PLG(4,3) + 0.538706E-3  * PLG(6,3) ) *
     $      COS( SR * (SEC-247.956) + 2.*umr*GLONG )
C  Exospheric temperature TINF/K  [Eq. A7]
      TINF = 1041.3 * ( 1. + T1+T2+T3+T5+T7814+T11+T12 ) * 0.99604
C
C TEMPERATURE DERIVATIVE AT LOWER BOUNDARY
C
C         F10.7 EFFECT
      T1 =  0.252317E-2 * DFA 
C        TIME INDEPENDENT
      T2 = -0.467542E-1 * PLG(3,1) + 0.12026 * PLG(5,1) 
C        ASYMMETRICAL ANNUAL
        CD14 = COS( DR  * (IDAY+8.45398) )
      T5 = -0.13324 * PLG(2,1)  * CD14
C        SEMIDIURNAL
      ZZ = PLG(4,3) * CD14
      T81 = -0.973404E-2 * ZZ
      T82 = -0.718482E-3 * ZZ
      T8 =(0.191357E-1 *PLG(3,3) + 0.787683E-2 *PLG(5,3) + T81) *C2TLOC
     3  + (0.125429E-2 *PLG(3,3) - 0.233698E-2 *PLG(5,3) + T82) *S2TLOC
C  dTn/dh at lower boundary  [Eq. A6]
      G0 = 0.166728E2 * ( 1. + T1+T2+T5+T8 ) * 0.951363
C
C NEUTRAL TEMPERATURE AT LOWER BOUNDARY 120KM
C
        CD9  = COS( DR2 * (IDAY-89.3820) )
        CD11 = COS( DR  * (IDAY+8.45398) )
      T1 = 0.568478E-3 * DFA
      T4 = 0.107674E-1 * CD9
      T5 =-0.192414E-1 * PLG(2,1) * CD11
      T7 = -0.2002E-1 *PLG(2,2) *CTLOC - 0.195833E-2 *PLG(2,2) *STLOC
      T8 = (-0.938391E-2 * PLG(3,3) - 0.260147E-2 * PLG(5,3)
     $       + 0.511651E-4 * PLG(6,3) * CD11 ) * C2TLOC
     $   + ( 0.131480E-1 * PLG(3,3) - 0.808556E-3 * PLG(5,3)
     $       + 0.255717E-2 * PLG(6,3) * CD11 ) * S2TLOC
C  Tn at lower boundary 120km   [Eq. A8]
      TLB = 386.0 * ( 1. + T1+T4+T5+T7+T8 ) * 0.976619
C  Sigma      [Eq. A5]
      SIGMA = G0 / ( TINF - TLB )
      RETURN
      END
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE TEMODEL(Tein,  Ht, TLoc, gmLat,ggLat, aMon,  rFlx, aKp)
! New Te, Ti [,Tn] model, August 2004,..    8/6/06: TeModela.f90, without debug code.
!                          For use with IRI;  see 
! Adv.Space Res Vol.38,#11, 2587-2595, 2006.
!INPUTS:
!    Ht   = km above the surface of the earth.
!    TLoc = Local solar time in hours (noon= 12.0).
!    gmLat= GeoMagnetic Latitude(deg) of the 'observing' point, at height Ht.
!               Thus calculations with varying Ht give a true, vertical Te profile.
!    ggLat= True (GeoGraphic) Latitude, in deg, corresponding to gmLat.  This is
!               to get correct solar changes.     [A zero value assumes ggLat= gmLat]
!    aMon = 1.0 (Jan 1st) to 12.97 (Dec31).    [A zero value uses aMon= 3.7,=equinox]
!    rFlx = 10.7cm solar Flux sFlx.            [A zero uses a mean value, F10.7= 120]
!    aKp  = Mean value of Kp over last 6 hr.             [A zero value uses aKp= 1.5]
!  Optional: Tein(3) = the neutral exospheric temperature  Tn
!                      (e.g. from MSIS; a value for Tn is used to calculate Ti.)
!OUTPUT: Tein(1) = calculated model value of Electron Temperature  Te.
!        Tein(2) = corresponding value for the Ion Temperature  Ti.
!        Tein(3) = Neutral Exospheric Temperature  Tn.  If a 'reasonable' value of Tn
!                  was supplied (i.e. Tein(3)> 100), this value is used and returned;
!                  otherwise Tn is estimated, and returned as  -Tein(3).
!
!CALCULATED TEMPERATURES correspond to a field-aligned heat flow, agreeing closely
!   with current theory (ASR paper ............).  Results depend on conditions near
!   the 'base' of the field line, at a height of 400 km in the upper ionosphere; the
!   geomagnetic latitude at this point is calculated from Ht and gmlat.  Calculation
!   of sunrise/set at this point assumes that the difference ggLat-gmLat is constant.
!   Internally, calculations use FLat,RLat= GeoMag,GeoGr lats of field line at 400km.
!
!ALTERNATIVE USEAGE:
!(A) For Te along a fixed field line: If the input value of Ht is negative, gmLat is
!       taken as the latitude of the base point(at 400km) of the required field line.
!            Results at varying Ht then give the changes in Te along that field line.
!              If |Ht| > height at the top of the field line (heq), then Te = Ti = 0.
!(B) To use the smoothed sunspot number Rz in place of S10.7 (the normal sFlx), call
!       TeModel with rFlx= -Rz.  This is used to calculate sFlx from standard equns.
!====================================================================================
!-USED: Eqns in JGeoRes Vol103,p2261,1998;  
!-------with updates from  Adv.Space Res Vol.38,#11, 2587-2595, 2006.
!   Comments give eqn nos from JGR,1998;  except Nxx is eqn (xx) from ASR,2006 paper.
!-ADDED: Variables polz and dayz, to ensure smooth varns over pole, equator.
!-ALLOWS for conjugate heating, using mixing factor xh= h/heq to match Te at equator.
!-J.E.Titheridge, March 2006.  Send comments/problems to j.titheridge@auckland.ac.nz.
!Copyright: John Titheridge, Physics Department, University of Auckland, New Zealand.
      real:: Tein(3), gmLat, ggLat                         
      real:: TD(5)=(/1250., 3460.,-3140.,  0.27, -0.60 /)  !Rational fit to T vs Lat.
      real:: GD(5)=(/ 2.09, -3.00,  1.13, -1.20,  0.835/)  !Rational fit to G vs Lat.
      real:: AD(5)=(/ 0.22, 0.016, -0.23, -2.08,  1.22 /)  !xD()= Day fits, to T,G,A.
      real:: TN(5)=(/ 990.,-1520., 2200., -1.13,  0.88 /)  !xN()= Night fits,  T,G,A.
      real:: GN(5)=(/ 0.63, -1.20, 0.615, -2.36,  1.49 /)
      real:: AN(5)=(/0.396, -0.81, 0.417, -2.17,  1.22 /)
      real,save:: TcD,TcN,GcD,GcN,AoD,AoN, TNit,TDay, FLine,Flat,dFlat  !These values
      real,save:: heq,polz,dayz, sRise,slong, Bnu, TKp=0.              !are saved for
      real,save:: zHt,zLoc,zLat,zMon, zKp,zRise, cLat,sL2  
      real,save:: TDv(2), TNv(2), TDiur(2), TNa=0.          !conjugate, diurnal varns
      logical,save:: NewFline=.false.                       !Flags changed field line
      data zHt,zLoc,zLat,zsLat, zMon,zKp,zRise / 7*-99. /   !Save to check if changed
      parameter(Pi=3.14159265,d2r=Pi/180., ex=2./7., ee=1.e-9)
      parameter(Re=6371.2, ho=400, Ro=1+ho/Re, Ro2=Ro**2)  
      rmLat= gmLat ;  if (gmLat==0.) rmLat= ggLat          !GeoMagnetic Latitude, deg
             if (abs(rmLat)>89.8) rmLat= sign(89.8, rmlat)       !6'06 limit to 89.8d
      drLat= ggLat-rmLat   ;   if(ggLat==0.) drLat=0.      !GeoGr-GeoMag (field line)
      qMon= aMon - 3.7     ;   if (aMon==0.) qMon= 0.            !change from equinox
      sFlx = rFlx                                                !-ve =sunspot number
      if(rFlx<0.)sFlx= 67. -0.572*rFlx +(.0575*rFlx)**2 +(.0209*rFlx)**3    !-Rz=>Flx
      dFlx= sFlx -120.     ;   if (sFlx==0.) dFlx= 0.            !change from Flx=120
      dKp = aKp - 1.5      ;   if (aKp ==0.) dKp = 0.            !change from aKp=1.5
      aHt = abs(Ht)        ;   Rh  = 1 + aHt/Re                        !Height in km.
!======================== A: Calculate FIELD LINE constants (once only: if L changes)
      IF(abs(rmLat-zLat)>ee.OR.(abs(aHt-zHt)>ee.AND.Ht>0.))THEN  !=== NEW FIELD LINE ====
        if (Ht >0.) then                             !for true VERTICAL CALCN -------
          cLat2= cos(rmlat*d2r)**2 *Ro/Rh                       !corresp Lat at 400km
          cLat = sqrt(min(cLat2,.999))                          !this Lat varies w ht
          FLat = sign(acos(cLat), rmLat) /d2r                   !Field Lat(400km) deg
        else                                         !calcns along FIELD LINE -------
          FLat = rmLat                               !Field-line Lat, at 400 km (deg)
          cLat = abs(cos(FLat*d2r))                  !cos(Lat400), field-line calculn
        end if
       FLine= Ro /cLat**2                                      !L value of field line
       heq  = (FLine - 1.) *Re                                       !ht over equator
       polz = cLat**0.65                                             !=>zero at poles
       Bnu  = 5./(2*FLine - Ro)                              !New eqn N7, replaces Bh
       dFlat= abs(FLat) - 46.                                        !for hilat corrn
       NewFline=.true.                                       !Flag changed parameters

!=====                       Calculate To, Go, A  for DAY, NIGHT  (eqn 19)
       sL2 = 1. - cLat**2  ;  sL4 = sL2*sL2                        !sL2=sin^2(Lat400)
       TcN = (TN(1) +TN(2)*sL2 +TN(3)*sL4) / (1 +TN(4)*sL2 +TN(5)*sL4)
       GcN = (GN(1) +GN(2)*sL2 +GN(3)*sL4) / (1 +GN(4)*sL2 +GN(5)*sL4)
       AoN = (AN(1) +AN(2)*sL2 +AN(3)*sL4) / (1 +AN(4)*sL2 +AN(5)*sL4)
       TcD = (TD(1) +TD(2)*sL2 +TD(3)*sL4) / (1 +TD(4)*sL2 +TD(5)*sL4)
       GcD = (GD(1) +GD(2)*sL2 +GD(3)*sL4) / (1 +GD(4)*sL2 +GD(5)*sL4)
       AoD = (AD(1) +AD(2)*sL2 +AD(3)*sL4) / (1 +AD(4)*sL2 +AD(5)*sL4)
      END IF   !---------------- END of setup for new field line ------------------------

      Te = 0.  ;  Ti = 0.  ;   if (aHt>heq+1.) goto 99          !Ht is above field line
!---------------------------                   Adjust for SOLAR and MAGNETIC activity
      ToN = TcN*(1 + dFlx/330)  ;   GoN = GcN * TcN/ToN         !Solar Varn, as percent
      ToD = TcD*(1 + dFlx/330)  ;   GoD = GcD * TcD/ToD
      IF (dFLat>0. .AND.abs(dKp)>ee) then                        !FLat >46deg, & Kp/=1.5
       if (abs(aKp-zKp)>ee.OR.abs(rmLat-zlat)>ee) then            !Lat or Kp changed:
          TKp= dKp*(37+1.33*dFLat-37*cos(0.135*dFLat)) ;  endif
       ToN = ToN + TKp  ;   ToD = ToD + TKp          !HiLat: Kp increases To, eqn 18a
       GoN = GoN *(1. - TKp/ToN)**2.5                !& large decrease of Go, eqn 18b
       GoD = GoD *(1. - TKp/ToD)**2.5
      END IF
!=========================== END of CONSTANT (Field-Line) Calculations ==============
!
!======================== B: Calculate SOLAR: sRise [for FLine at 400km] ============

      rslat = (FLat + drLat) *d2r                       !geog/solar=FLat+(ggr-gmag),rad
      IF (abs(aMon-zMon)>ee.OR.abs(rslat-zsLat)>ee) THEN !NEW SUNRISE time -------------
         sunDec= 0.41015*sin(qMon*0.5236)                           !solar decln(rad)
         csLat = max(cos(rslat), .0001)                                   !keep +ve,,
         slong = 1. - cos(rslat-sunDec) /(csLat*cos(sunDec))              !to stop /0
         if (slong>= 1.) then         ;  sRise =12.                           !no day
           else if (slong<=-1.) then  ;  sRise = 0.                           !no nit
           else   ;  sRise = 12. - acos(slong)*3.82                   !sunrise (L.T.)
          end if
         dayz= 0.7/ max(0.7, abs(slong))                            !<1 if no day/nit
         if(slong>.7) dayz= (1+ (dayz-1.)*polz  + dayz) /2.  !6'06 smooth winter pole
      END IF
!======================== C: Calculate DAY, NIGHT temps, at height h [Eqns 13] ======
      IF (abs(aHt-zHt)>ee.OR.NewFline) THEN              !NEW HEIGHT or FIELD LINE =====
       dhr = (heq - ho)/Ro2 - (heq - aHt)/Rh**2                         !from eqn 13a
          if (dhr<0.) dhr= dhr/(1. + (ho-aHt)/300)                 !limits to h>100km
       TNit= ToN *abs(1 + Bnu*dhr *GoN/ToN)**ex                         !New eqn. N6a
       TDay= ToD *abs(1 + Bnu*dhr *GoD/ToD)**ex
          z =(aHt-ho)/(1500*AoN +10) ;  Z = z*sqrt(abs(z))              !-ve below ho
        TNit= TNit*(1. + AoN*Z/(1+abs(Z)))                      !gives Te*A atlarge h
          z =(aHt-ho)/(1500*AoD +10) ;  Z = z*sqrt(abs(z))
        TDay= TDay*(1. + AoD*Z/(1+abs(Z)))                             !in New eqn N7
      END IF
!======================== D: Add DIURNAL VARIATION [Eqns 20-22, +Conjugate effects] =
      TLoc = mod(TLoc+24., 24.)                                                 !in 0-24
      IF ((abs(TLoc-zLoc).gt.ee).or.(abs(sRise-zRise).gt.ee)
     &.or.(abs(zLat-rmLat).gt.ee).or.(abs(TNa).lt.ee)) THEN                 !NEW time, lat
      dtr = 1.2 + 0.5*sL2                                                    !eqn 20b
      sr  = sRise
      TNx = (mod(TLoc+11.,24.) -11.) *(0.8*sL2-1.4)*dayz                     !for dtN
      DO k = 1, 2                                    !k= 1,2 for local,conjugate calcn
        ts = min(1., 0.15*sr)                                                !eqn 20a
        if (TLoc<12.) then
           D = (sr - 0.5*ts - TLoc) / dtr                              !a.m., eqn 21b
         else                                                  !D is <0 day, >0 night
           D = (TLoc - ts + sr - 24.)/(dtr + 0.9*polz)                 !p.m., eqn 21c
         end if                                     !D jumps at 12h; TNv jumps at 13h
        TDiur(k)= 1./(1 + exp(D*3.2*polz) )         !Day=1,Nit=0(eq 21a); 0.5 @ poles
          if (abs(slong)>0.7) TDiur(k)= TDiur(k) *dayz                 !=0 if all nit
          if (slong <-0.7) TDiur(k)= TDiur(k) + 1-dayz                 !=1 if all day
        TDa = (TLoc - 12.5)/ max(12.5-sr, 3.)                !stop growth at sris>9.5
        TDv(k)= 0.97 + 0.22*TDa**2 *dayz *polz                               !eqn 22a
        dtN = exp(TNx /max(5., sr-3.) )                      !Changes at 13h. Can =11
        TNa = 0.83 + 2.1*dtN/(14.+dtN)  *polz                                !eqn 22b
        TNv(k)= min(TNa, 2.)                                                 !eqn 22b
       sr = 12. - sRise                               !store conjugate values, at k=2
       slong= -slong                                  ![changes  for conjugate calcn]
      END DO                                           ![dayz,polz= same at conjugate]
      END IF
      xh2 = 0.5*min(aHt,heq)/heq              !weight for conjugate temps (=> nit incr)
      xh1 = 1 - xh2                                            !!weight for local temps
      TDiurn= TDiur(1)*xh1 + TDiur(2)*xh2                            !=0,night to 1,day
      TeN = TNit * (TNv(1)*xh1 + TNv(2)*xh2)                                   !eqn 22b
      TeD = TDay * (TDv(1)*xh1 + TDv(2)*xh2)                                   !eqn 22a
      Te  = TeN + (TeD - TeN) *TDiurn                                          !eqn 21a
!=========================== END  Te.   Calculate [Tn and] Ti =======================
      Tnn = Tein(3)                           !Positive = 'true' value of Tn, from MSIS
      if (Tnn<100.) then                                         !calculate Tn, eqns 26
       T00 = 956 + 3.1*dFlx - 100*cos(0.028*abs(FLat) - 0.5)                     !26a
       T12 =1030 + 3.8*dFlx -  30*cos(0.055*abs(FLat) - 0.5)                     !26b
       Tnn =0.5*(T00+T12) + 0.68*(T00-T12)*(cos(0.2618*TLoc-0.7) -0.1)           !26c
       Tein(3)= -Tnn                                    !Negative = approximate model
      end if
      if(Te<1.5*Tnn) Te= 0.9952*Te*(1 + (Tnn/Te)**8)**0.125      !limit; match at 1.5Tn
      Rion= (0.75*Tnn - 75) /aHt                  !relative Ht for Ti=(Te+Tn)/2, eqn 24
      A  = 0.08 + 0.18*TDiurn                                       !Ti<Te at large hts
      Ti = Tnn + (Te - Tnn)/(1. + (A + (Rion**2)**2 )*Rion )                    !eqn 25
      Tein(1)= Te  ;  Tein(2)= Ti
!--------------------------- Save parameters, to check for changed values -----------
      zsLat= rslat
99    zHt  = aHt   ;  zLoc = TLoc  ;  zLat= rmLat             !exit (Te=0. if Ht>heq)
      zMon = aMon  ;  zRise= sRise ;  zKp = aKp
      RETURN
      END       !! of TeModel =============================================================
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       subroutine topscale1(hsip,scalh,xxx)
C
C Calculation of topside altitude hsip(NmF2/e) and topside scale height Scalh=hsip-hmF2
C
C T.L. Gulyaeva......................................................Feb. 2010
C
      DIMENSION OARR(50)
      common   /block1/hmf2,xnmf2,hmf1,F1REG
     &         /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU,HTOP,QF0
     &/PLAS/HPL,XNEPL,XKP,ICALLS,RUT,RLT,OARR,TECI,TCB,TCT,TCPL,TEC
     &/DEM/W,GLON,HOUR,DAYNR,CN1000,ALON,BLON,ENRE,HSC
      logical    f1reg
      integer daynr
C
       xxx=xnmf2/2.718282   ! NmF2/e
      istart=int(H05TOP)
      ht1=H05TOP
       xne1=XNETOP

      ht2=float(istart)  
      do i=istart,1336
      ht2=ht2+1.0
      xne2=xxe6(ht2) 
      if ((xxx.le.xne1).and.(xxx.gt.xne2)) then
          hsip=ht1+(ht2-ht1)/(xne2-xne1)*(xxx-xne1)
      exit
      else
          ht1=ht2
          xne1=xne2
      endif
      enddo
      scalh=hsip-hmF2
      RETURN
      END
C =================================================================================================
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       subroutine plashtop
C
C Calculation of topside altitude h05pl=h05top(0.5*NmF2) using Plasmasphere model
C
C T.L. Gulyaeva......................................................Feb. 2010
C
      DIMENSION OARR(50)
      common   /block1/hmf2,xnmf2,hmf1,F1REG
     &         /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU,HTOP,QF0
     &/PLAS/HPL,XNEPL,XKP,ICALLS,RUT,RLT,OARR,TECI,TCB,TCT,TCPL,TEC
     &/DEM/W,GLON,HOUR,DAYNR,CN1000,ALON,BLON,ENRE,HSC
      logical    f1reg
      integer daynr
C
c-       xxx=xnmf2/2.718282   ! NmF2/e
      xnetop=xnmf2/2.0         ! 0.5 NmF2
      istart=int(hmf2)
      ht1=hmf2
       xne1=xnmf2
      step=100.
      ht2=float(istart)  
C
    1 ht2=ht2+step
      if (ht2.lt.1336.) then       !>>>>>>>>>>
      xne2=xe1(ht2)   ! 
      if (abs(xne2-xnetop).lt.0.1) then
      H05TOP=ht2
      goto 2
      endif                      
      if ((xnetop.ge.xne2).and.(xnetop.lt.xne1)) then  !,,,,,,,,,,,,,
      ht2=ht2-step         !
      step=step/10.
      xne2=xne1
      if (step.lt.1.) then
      h05top=ht2
      goto 2
           else
      goto 1
                  endif
                                              else   !............
      xne1=xne2
      ht1=ht2
      goto 1  
      endif                                    !..............
                        endif !>>>>>>>>>>>>>>>>>>>>>
    2 RETURN
      END
C =================================================================================================
      SUBROUTINE SUBDLOGH(dlogne,amlat,sg,rz,dh)
C....................................................................................Mar. 2011
C
C Model of dh=dlog_hmF2, in terms of log10(NeF2/NmF2), abs(MLAT), deg., SG=Sday, Rz
C
C  (1) 2 levels of solar activity: Rz= 20,  120
C  (2) 5th order polinomial for geomagnetic latitudes (N=S):0,10,20,30,40,50,60,70,80,90.
C  (4) 4 seasonal grids: 1 equinox(sg=90deg), 2 summer (sg=180), 
C                        3 equinox (sg=270), 4 winter(sg=360)
C Gulyaeva, T. Empirical model of ionospheric storm effects on the F2 layer peak height 
C        associated with changes of peak electron density, J. Geophys. Res., 117, A02302, 
C        doi:10.1029/2011JA017158, 2012. 
C
      DIMENSION YY(3),XX(3)
     &,BR(6,3,2) ! 6  beta coef, 3 seasons, 2 Rz levels
     &,AR(6,3,2) ! 6 alpha coef, 3 seasons, 2 Rz levels
      COMMON /BLA1B1/ALPHA,BETA

      DATA rad/0.01745329/
      DATA AR/
C Polinomial <alpha> coefficients A1, A2, A3, A4, A5, A6 for Rz=20:
C Equinox:
     * 708.327,-1527.953,  1103.833, -304.492,  34.062,  -4.422
C Summer 
     *,-64.5256, 254.8368, -291.2572,  97.7416,  8.4343, -3.9076
C Winter:
     *,99.9013, -61.4048,-120.2377,  92.2898,  -4.6909, -2.9850
C Polinomial <alpha> coefficients A1, A2, A3, A4, A5, A6 for Rz=120:
C Equinox:
     *,271.9513,-261.3368,-200.8858, 249.5253, -47.7269,  -2.0143
C Summer    
     *,-422.818,1012.842, -890.418,  347.730,  -54.262,    1.859
C Winter:      
     *,145.3128,-128.5541,-131.3022, 148.7525, -29.4966,  -1.1053/
      DATA BR/
C Polinomial <beta> coefficients B1, B2, B3, B4, B5, B6 for Rz=20:
C Equinox:
     * 3.0154,  -6.4793,   4.8852,  -1.3451,   0.0490,  -0.1681
C Summer    
     *,-5.8795,  17.7169, -18.0480,   7.5002,  -1.1586,  -0.1527
C Winter:      
     *,24.7256, -57.3807,  45.7255, -13.6397,   0.9003,  -0.1040
C Polinomial <beta> coefficients B1, B2, B3, B4, B5, B6 for Rz=120:        
C Equinox:  
     *,10.3154, -23.8605,  19.0188,  -5.5255,   0.2329,  -0.1574
C Summer    
     *,-3.3154,   9.3716,  -9.1249,   3.2476,  -0.1863,  -0.2255
C Winter:
     *,12.9551, -28.8785,  22.2281,  -6.2295,   0.1901,  -0.0897/
C
C----------------------------------------------------------------------
C
      ABMLAT=ABS(AMLAT)/100.
      dlogne=dlogne*1000.
      IR=1    ! LSA
      xi=abmlat
  25  CONTINUE
C
      DO LS=1,3
      B1=BR(6,LS,IR)
      B2=BR(5,LS,IR)
      B3=BR(4,LS,IR)
      B4=BR(3,LS,IR)
      B5=BR(2,LS,IR)
      B6=BR(1,LS,IR)
      YY(LS)=B1+xi*(B2+xi*(B3+xi*(B4+xi*(B5+xi*(B6+xi)))))
      A1=AR(6,LS,IR)
      A2=AR(5,LS,IR)
      A3=AR(4,LS,IR)
      A4=AR(3,LS,IR)
      A5=AR(2,LS,IR)
      A6=AR(1,LS,IR)
      XX(LS)=A1+xi*(A2+xi*(A3+xi*(A4+xi*(A5+xi*(A6+xi)))))
      ENDDO
C Apply seasonal interpolation
      p0=(2.*YY(1)+YY(2)+YY(3))/4.
      p1=(YY(3)-YY(2))/2.
      p2=(YY(2)+YY(3)-2.*YY(1))/4.
      YB=p0+p1*cos(sg*rad)+p2*cos(2.*sg*rad)
      r0=(2.*XX(1)+XX(2)+XX(3))/4.
      r1=(XX(3)-XX(2))/2.
      r2=(XX(2)+XX(3)-2.*XX(1))/4.
      XA=r0+r1*cos(sg*rad)+r2*cos(2.*sg*rad)
C
      IF (IR.eq.2) THEN
      goto 26
                ELSE  
      IR=2
      YB1=YB
      XA1=XA
      goto 25
                ENDIF
   26 CONTINUE
C Solar activity interpolation
      alpha=XA1+(XA-XA1)*(RZ-20.)/100.
      beta=YB1+(YB-YB1)*(RZ-20.)/100.
C
      DH=(beta*dlogne+alpha)/1000.
      RETURN
      END
C ===============================================================
