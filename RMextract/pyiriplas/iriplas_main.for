CC    T.L. Gulyaeva     Iriplas20                            Oct. 2020
CC
CC    T.L. Gulyaeva     Iriplas7                             May. 2019
CC
CC    T.L. Gulyaeva     Iriplas7                             May. 2017
CC
CC    T.L. Gulyaeva     Iriplas5c                            Feb. 2016
CC
CC    T.L. Gulyaeva     Iriplas5                             Jun. 2015
CC
CC    T.L. Gulyaeva     Iriplas4                             Aug. 2014
CC
CC    T.L. Gulyaeva     Iriplas3                             Apr. 2014
CC
CC    T.L. Gulyaeva     Iriplas1                             Dec. 2011
CC
CC    T.L. Gulyaeva     Isomain7                             June 2009
CC
CC                         ISO_IRI-Plas MAIN PROGRAM 
CC
CC              EARTH'S MODEL OF IONOSPHERE AND PLASMASPHERE
CC
CC                INCLUDING IRI/SMI/ISOMAIN3 SUBROUTINES
CC
CC
CC     PRODUCTS: ELECTRON DENSITY PROFILES 
CC                           AND TOTAL ELECTRON CONTENT
CC     AT ALTITUDES OF 65 TO 35000 KM AT ANY LOCATION OF THE EARTH
CC
CC     INPUT: YEAR, MONTH, DAY, HOUR, LOCATION, SUNSPOT NUMBER, 
CC     MAGNETIC AP/KP-INDEX, SOLAR RADIO FLUX AND SOLAR PROXY INDEX (8 options)
CC     OPTIONAL: foF2 AND/OR hmF2/OR M3000F2 (IF ZERO, CCIR MAPS ARE USED
CC     INCLUDING IRI_STORM MODEL)
CC     NEW OPTION: TEC INPUT (IF TEC=0, IRI-Plas TEC CALCULATION)
CC                
CC     OUTPUT: ELECTRON DENSITY PROFILES AND 
CC             ELECTRON CONTENT AT SELECTED
CC             ALTITUDES FROM 80km TO THE PLASMAPAUSE (hpl < 35000km).
CC             STANDARD SET OF PARAMETERS / Ne(h),fN(h) profiles /
CC
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ISOMAIN3.FOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Changes ____________________________________________________Oct. 2006
C
C Temodel-2006 subroutine (by J. Titheridge) is included for Te, Ti, Tn above 400 km
C
C Changes of ISOMAIN3.FOR___________________________Feb. 2005
C
C Dependence of foF2 and M3000F2 on sunspot number Rz is used 
C replacing interpolation of CCIR maps versus ionospheric IG-index used by IRI-2000/1995
C
C Chebishev polinomial approximation is introduced (SUBROUTINE CHEBISH) for
C topside ratio of half-width to F2 layer peak height depending 
C on geomagnetic latitude and local time
C
C Changes from ISOMAIN2.FOR_________________________July 2003
C
C Topside IRI-Bent model is corrected with QF-factor by fitting
C
C Ne(h) profile to h05top height at Ne05=0.5*NmF2 deduced from
C ISIS and IK19 topside electron density profiles
C
C contains IRIT13, IONCORR, IRI_TEC----D.Bilitza-- Oct 20, 1995
C INTEGR, CIRA-86 , STORM, APF ______________________ IRI-2001
C SMI-source: CORDM, new: GKPM using results of APF 
C ___________________________________T.L.Gulyaeva___Sep. 2001
C
C ----------------------------------------------------------------
C
C
      SUBROUTINE IRI_PLAS_MAIN(path,ALATI,ALONG,JMAG,IYYYY,MMDD,HOURS,
     & OUTF)

c     change code such that it takes input from array of booleans (if necessary) + above values
c     write output to array instead of file. use irisub as example code
c     Remove all !*
c     Parameter list definition are not running well withf2py, separate 1 per line (PARAMETER (a=1,b=2...)
c     becomes PARAMETER a=1, PARAMETER b=2 ....
c     implied-DO list not supported, change DATA (C(1,2,J),J=1,81)/ to DATA C(1,2,:)/
c     Add predefined settings instead of using read!
c     Add datapath also to igfr*.for and Iris*for
      
C-----------------------------------------------------------------
c     Program for numerical integration of IRI profiles from h=100km
C     to h=alth. 
C     h=hpl IS ALLOWED WITH SMI XXE6(H) AMENDMENT
C     INPUT:  ALATI,ALONG  LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C     jmag     =0 geographic  =1 geomagnetic coordinates
C     iyyyy,mmdd date as yyyy and mmdd (or -ddd)
C     hour     decimal hours LT (or UT+25)
C+ NEW_2020_Proxy2(SSN2) using GMF2(SSN2) => foF2(SSN2) and hmF2(SSN2) models
C* NEW_2017 solar proxy indices driving IRI-Plas model:
C*    jinp =1(SSN1); 2(SSN2); 3(F10.7); 4(GEC); 5(TEC); 6(IG); 7(MgII); 8 (Lym-a); 1(SSN1 & F10.7: defalt)
C*    jinp =0  previous version using input of gec_rz.dat file
C     jout = 0 NO OUTPUT Ne(h)& fN(h) PROFILES
C     jout = 1 OUTPUT Ne(h)& fN(h) PROFILES
C     jstr = 0 no foF2 STORM model 
C     Jstr = 1 foF2 STORM model included 
C     juc  = 0 input foF2>0 and/or hmF2>0 (or M3000F2*100.+1000.)
C     juc  = 1 option of CCIR maps of foF2 and M3000F2
C     juc  = 2 option of URSI foF2 map 
C     hbeg,hend   upper and lower integration IRI limits in km
C     RZS   sunspot number (0 =< SSN1 < 1000.)
C     RZS=-1 : sunspot number (SSN1 or SSN2) from proxy.dat file used
C     RZS+1000 means input of F10.7 
C     XKP geomagnetic  Kp-index (option of input)
C     XKP=-1 : accumulated Kpm produced by GKPM subroutine from aprz.dat
C     input: foF2>0 and/or hmF2>0 and/or TEC>0 [TEC*1E-16, m-2]
C     OUTPUT: 
C     XHI solar zenith angle 
C     COV solar radio flux averaged for 81 days preceeding given day
C     Kpm geomagnetic index driving plasmasphere model
C     Api geomagnetic index driving foF2 storm model
C     YMLAT magnetic latitude
C     YMLONG magnetic longitude
C     YMODIP modified dip latitude
C     foF2, hmF2 after input or CCIR (URSI) calculations
C     Nes   Ne(1336 km) in cm-3
C     HPL,XNEPL PLASMAPAUSE HEIGHT AND ELECTRON DENSITY
C     ECbot Ionospheric Electron Content in M-2 FOR [80km, hmF2]
C     ECtop Ionospheric Electron Content in M-2 FOR [hmF2, 1336km]
C     ECpl  Plasmaspheric Electron Content in M-2 FOR [1336km,HPL]
C     TEC   Total Electron Content in M-2 FOR [80km,HPL]
C     TAU slab thickness, km (ratio TEC/NmF2)
C     h05bot height below hmF2 at Ne=0.5*NmF2
C     h05top height above hmF2 at Ne=0.5*NmF2
C     Ne(h) & fN(h) PROFILES FOR [80km, HPL] // optional //
C     OARR(37) = ALATI  OARR(38) = ALONG  Geographic coordinates
C     OARR(40) = AMLONG  OARR(42) = AMLAT   Geomagnetic coordinates
C
C------------------------------------------------------------------
      dimension   outf(0:11,0:120),oarr(50),indap(13),akp(13)
      DIMENSION  icurkp(0:7),iprekp(0:7)
      logical     jf(30),F1REG
      REAL        RZS,COV,RZS1,RZS2
      integer daynr,numhei
      character*4  xxpr,yypr     
      COMMON /BXY/xxpr,yypr      
      common   /block1/hmf2,xnmf2,hmf1,F1REG
     &/PLAS/HPL,XNEPL,XKP,ICALLS,RUT,RLT,OARR,TECI,TCB,TCT,TCPL,TEC
      COMMON   /CONST/UMR,PI/ARGEXP/ARGMAX
     &         /CONST1/HUMR,DUMR,FAP,INDAP,AKP,COV,RZS1,RZS2
     &/dem/W,glong,hour,daynr,cn1000,alon,blon,ENRE,HSC
      COMMON /FIKP/icurkp,iprekp,iyear_cur,imn_cur,idy_cur TEST
C# new common block:
     &         /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU,HTOP
      common /temod/month            
      CHARACTER(10) DD
      CHARACTER(10) TT
      CHARACTER(5) ZZ
C**********TEST
      CHARACTER*128 INFILE,OUTFILE,HEADER(4)
       CHARACTER*1 AINP
       CHARACTER*4 XYEAR
       CHARACTER*2 XMN,XDY,AMNUT,ADYUT,AMNLT,ADYLT
C+2020++++++++++++++++++++++++++++++++++++
Crem	  character*4 F2param
Crem	character*1 fh                       
Crem        real Lati, longi
Crem     	  integer INY_R12,INY_R
Crem        integer ninp,nout
Crem        data ninp,nout/2,3/		   
	integer dayut, monthut, yearut
Crem	+, ndayut
	integer daylt, monthlt, yearlt
	real hour, SLT, UT
C//	real Ri, SSN_R12, RR 
	real foF2,hmF2
C//	real model_GMF2R
C+2020++++++++++++++++++++++++++++++++++++
      EXTERNAL          XE1,XE2,XE3,XE3_1,XE4,XE5,XE6,XXE6,TEDER,ABLON
C      DATA HEADER/' YEAR MMDD UThr LThr  XHI   SSN   COV  Kpm    L  ',
      DATA HEADER
     & /' UTYRMMDD UThr LTYRMMDD LThr  XHI   SSN   COV  Kpm    L  ',
     &            '  Glati Glong Mlati Mlong MoDip  hmF2  foF2',
     &            '    NmF2   Nes       QF     MLT   ECbot  ECtop',
     &            '   ECpl   TEC    TAU  h05b  h05t   Hsc '/
      CHARACTER path*200
      LOGICAL OUTPUT
      OUTPUT=.FALSE.
      call set_datapath(path)
C**********TEST
c
  111    ICALLS=0
C**** TEST
      CALL DATE_AND_TIME(DATE=DD,TIME=TT,ZONE=ZZ)
      TT(5:)='      '
      WRITE(*,*) 'PC Date: Year,Month,Day = ',DD,'Time = ',TT
C
       XYEAR=DD(1:4)
       XMN=DD(5:6)
       XDY=DD(7:8)
       read(xyear,*) ryear
       read(xmn,*) rmn
       read(xdy,*) rdy
       iyear_cur=int(ryear)
       imn_cur=int(rmn)
       idy_cur=int(rdy)
C
c
    1 CONTINUE
c      READ (11,12,END=2,ERR=2) IYR,MMDD,JINP,JOUT,JSTR,JUC,HOURS,RZS  
c     &,UKP,JMAG,ALATI,ALONG,HMF2I,FOF2,TECI
C*****SOME PREDEFINED SETTINGS SEE INPUT_FILE.FMT FOR MEANING
      JINP=1
      JOUT=1
      JSTR=1
      JUC=2
      RZS=-1
      UKP=-1
      HMF2I=0
      FOF2=0
      TECI=0
      DO 11 I=1,30
      JF(I)=.TRUE.
   11 CONTINUE 

      IF(MMDD.EQ.0.) GOTO 2
C
C Check availability of kpYEAR file:
C   Input file kpYEAR starts from 1948:
              if(MMDD.lt.0) then
                lda=-MMDD
                call MODA(1,iyyyy,imn,idy,lda,nrdaym)
        else
                imn=MMDD/100
                idy=MMDD-imn*100
                call MODA(0,iyyyy,imn,idy,lda,nrdaym)
        endif
C++
C++
      if (HOURS.ge.25.0) then
	UT=hours-25.0
         dayut = idy
	   monthut = imn
	   yearut = IYYYY
	    call rconvtime(0,dayut,monthut,yearut,along,SLT,UT,
     &                                 daylt,monthlt,yearlt)
      	               else
	SLT=hours
         daylt = idy
	   monthlt = imn
	   yearlt = IYYYY
	   call rconvtime(1,daylt,monthlt,yearlt,along,SLT,UT,
     &                                dayut,monthut,yearut)
                         endif
      call blet2(monthut,AMNUT)
	call blet2(monthlt,AMNLT)
	call blet2(dayut,ADYUT)
	call blet2(daylt,ADYLT)
C++     
C
      if((iyyyy.lt.1948).or.(iyyyy.gt.iyear_cur)) goto 21    
      IF ((iyyyy.eq.iyear_cur).and.(imn.gt.imn_cur)) goto 21 
      IF ((iyyyy.eq.iyear_cur).and.(imn.eq.imn_cur).and.     
     *   (idy.gt.idy_cur)) goto 21 
      XKP=UKP/10.
      goto 22
C
   21   write(*,100)
  100 format(1X,'Date is outside range of Kp-Ap indices file.',
     &     'put kp=2.0 ')
        XKP=2.0
C      READ(*,*) HEIBEG,HEIEND
   22      HBEG=80.
COR:
C New:
      jf(25)=.true.
      jf(27)=.true.
      gind=0.                     
C
      IF (rzs.lt.0.) GOTO 26       
C
	if (RZS.gt.999.) then
           cov=rzs-1000.        
           RZS1=1.076*COV-65.7817     
      	 RZS2=RZS1/0.7
      oarr(41)=cov
      jf(25)=.false. 
	goto 25
	endif
      	IF (ninpx.eq.2) THEN
      	rzs2=rzs     
      	rzs1=0.7*rzs2
	oarr(33)=rzs2
 				   ELSE
           rzs1=rzs   
           rzs2=rzs1/0.7
	oarr(33)=rzs1
	              ENDIF
C
   24	COV=0.9066*RZS1+62.6645            
   25 if (jinp.eq.0) then               
       GIND=0.97236*RZS1+0.18836
       oarr(39)=gind                  
      endif                             
C
      jf(27)=.false.                    
C
   26 continue
C++   Add subroutine to convert GEOMag -to- GEOGraf coords
      call subcoords(JMAG,IYYYY,LDA,ALATI,ALONG,AMLAT,AMLON,
     &                       AMODIP,DIP,AMAGBR)
C++
C
c  select various options and choices for IRI
c
C     jf(2)=.false.     ! no temperatures
C     jf(3)=.false.     ! no ion composition
Ctest jf(4)=.false.     ! Gulyaeva-B0
C     jf(5)=.false.     ! URSI-88 for foF2
C     jf(8)=.false.     ! input foF2 or NmF2
C     jf(9)=.false.     ! input hmF2 or M(3000)F2
C     jf(12)=.false.    ! konsol output to file (unit=12)
C     jf(17)=.false.    ! input of sunspot number Rzs
C     jf(18)=.false.    ! no output of Ne(h) & fN(h) profiles
C     jf(25)=.false.    ! OARR(41)=user input for daily F10.7 index
C     jf(26)=.false.    ! no STORM output
C      JF(2)=.FALSE.
      JF(3)=.FALSE.
      JF(4)=.FALSE.
      IF (JUC.EQ.2) JF(5)=.FALSE.
      IF (FOF2.GT.0.) JF(8)=.FALSE.
C++
      IF (jinp.eq.2) THEN
C...      call GMF2R12 for fof2
      JF(8)=.TRUE.
	ENDIF
C++
      IF (HMF2I.GT.0.) JF(9)=.FALSE.
C++
      IF (jinp.eq.2) THEN
C...      call GMF2R12 for hmf2
      JF(9)=.TRUE.
	ENDIF
C++
      IF (RZS.GE.0.) JF(17)=.FALSE.
      IF (JOUT.EQ.0) JF(18)=.FALSE.
      IF (JSTR.EQ.0) JF(26)=.FALSE.
      IF (TECI.GT.0.) JF(30)=.FALSE.
C********************************************************************
      tec = -111.
      tect= -111.
      tecb= -111.
      TCPL= -111.
       IF (HMF2I.LT.999.) THEN
       hmf2=hmf2i
       xm3000=0.
       OARR(2)=HMF2
       ELSE
       hmf2=0.
       xm3000=(hmf2i-1000.)/100.
       OARR(2)=xm3000
       ENDIF
      OARR(1)=FOF2
c
c     initialize IRI parameter in COMMON blocks
c
      abeg=hbeg
C      aend=hend
      AEND=1336.
      astp=aend-hbeg
        if (icalls.eq.0) then
      write(*,*) ' ISO_IRI-Plas parameters are being calculated '
      write(*,*) 'Ne, B0: Bottomside thickness is ',
     &      'obtained with Gulyaeva-1987 model.'
      endif
      call IRIS2020(JF,IYYYY,MMDD,HOURS,abeg,aend,numhei,OUTF,JINP
     &                ,ALATI,ALONG,AMLAT,AMLON,AMODIP,DIP,AMAGBR)                          
C     
C Add:
      if (JF(30)) GOTO 123  
C
C Add: 2nd call IRIS2015c: replace TEC input with 'input' of updated foF2, hmF2:
      jf(30)=.TRUE. 
      TECI=0.
      JF(8)=.FALSE.
      JF(9)=.FALSE.
      tec = -111.
      tect= -111.
      tecb= -111.
      TCPL= -111.
      icalls=icalls-1
C++
      call IRIS2020(JF,IYYYY,MMDD,HOURS,abeg,aend,numhei,OUTF,JINP
     &                ,ALATI,ALONG,AMLAT,AMLON,AMODIP,DIP,AMAGBR)                          
C End Add commands
C
  123    IF(JF(8)) FOF2=OARR(1)
      IF(JF(9)) HMF2=OARR(2)
      IF(JF(17))  RZS=OARR(33)
      AL=OARR(32)   
      IF (AL.GT.99.99) AL=99.99
      XHI=OARR(23)
      YMODI=OARR(27)
CCC
CC ENRE=OARR(32)
      GLATI=OARR(37)
      GLONG=OARR(38)                  
      XMLAT=OARR(42)
      XMLON=OARR(40)
      GMLT=OARR(44)
      HTOP=OARR(46) 
      EN1000=CN1000/1.0E+06
      TOT=TCB+TCT+TCPL
      TAU=(TOT/XNMF2)*1.0E+13
      XNMF2=XNMF2/1.0E+06
      FOF2=SQRT(XNMF2/1.24E+04)
	ifap=nint(fap)
c>>      XNEPL=XNEPL/1.0E+06
C*********TEST
C TWO OPTIONS OF OUTPUT: <7> ONLY KEY PARAMETERS; <8> FULL PROFILE
c      IF (JOUT.EQ.0) THEN 
c         GOTO 7
c      ELSE
c         GOTO 8
c      ENDIF
c    7 IF (ICALLS.EQ.1) THEN
c      WRITE(*,'(A)') TRIM(HEADER(1))//TRIM(HEADER(2))//
c     &                     TRIM(HEADER(3))//TRIM(HEADER(4))
c      ENDIF
c      WRITE(*,9) yearut,AMNUT,ADYUT,RUT,yearlt,AMNLT,ADYLT,RLT,XHI,RZS
c     &,COV,XKP,AL,GLATI,GLONG,XMLAT,XMLON,YMODI,HMF2,FOF2,XNMF2,EN1000
c     &,QF,GMLT,TCB,TCT,TCPL,TOT,TAU,HHALF,H05TOP,HTOP
c    9 FORMAT (2(1X,I4,2A2,F5.1),3F6.1,F4.1,1X,F5.2,6(F6.1),F6.2,F9.0
c     &,F8.0,1X,F6.3,F8.3,4F7.2,F7.1,2(F6.1),F7.1)
c    7    GOTO 1
c    8 WRITE(*,'(A)') TRIM(HEADER(1))//TRIM(HEADER(2))
c      WRITE(*,13) yearut,AMNUT,ADYUT,RUT,yearlt,AMNLT,ADYLT,RLT,XHI,RZS
c     &,COV,XKP,AL,GLATI,GLONG,XMLAT,XMLON,YMODI,HMF2,FOF2
c      WRITE(*,'(A)') TRIM(HEADER(3))//TRIM(HEADER(4))
c      WRITE(*,14) XNMF2,EN1000,QF,GMLT,TCB,TCT,TCPL,TOT,TAU,HHALF
c     &,H05TOP,HTOP
C*********TEST
c   13 FORMAT (2(1X,I4,2A2,F5.1),3F5.0,F4.1,1X,F5.2,6(F6.1),F6.2)
c   14 FORMAT (1X,F8.0,F8.0,1X,F6.3,F8.3,4F7.2,F7.1,2(F6.1),F7.1)
c      WRITE(*,15)
c   15   FORMAT(4X,'H',9X,'NE',9X,'FN',6X,'Te',7X,'Ti',7X,'Tn')
C 
   18 CONTINUE
c     MM:this goto starts over again, for infinite questioning, remove 
c     GOTO 1
    2 CONTINUE
      if (icalls.eq.0) then
      write(*,*) 'INPUT FILE IS NOT IN YOUR DIRECTORY '
      endif
c      write(*,*) 'finished iriplas calculations'

c      close(unit=11)
      close(unit=25)
Crem	GOTO 111
      RETURN
      END
           subroutine set_datapath(path)
c----------------------------------------------------------------
c set common data path for data files
c----------------------------------------------------------------
           CHARACTER path*200
           common /path/datapath
           character datapath*200
           datapath=path
           RETURN
           END
