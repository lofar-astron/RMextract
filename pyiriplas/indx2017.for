C INDEX SUBROUTINES:___INDX2017.FOR___T.L.Gulyaeva___May. 2017  
C
C INDEX SUBROUTINES:___INDX2015c.FOR___T.L.Gulyaeva___Feb. 2016  
C
C INDEX SUBROUTINES:___INDX2015.FOR___T.L.Gulyaeva___Dec. 2014  
C
C ==   AMENDMENTS ...................................Dec. 2014
C STORM RELATED SUBROUTINE FROM IRI2001: = = APF == is replaced
C by: == SUBAPF == which enters geomagnetic Kp and Ap 3-h indices
C and solar daily Rz and F107 indices from standard = kpYEAR = file
C for given 'YEAR' and preceding 'YEAR-1' for calculations of 81-days 
C indices R81 and F81 instead of former APRZ.DAT file (cancelled input)
C
C INDEX SUBROUTINES:___INDX2014.FOR___T.L.Gulyaeva___Apr. 2014  
C
C INDEX SUBROUTINES:___INDX2013.FOR___T.L.Gulyaeva___Dec. 2013  
C
C INDEX SUBROUTINES:___INDX2011.FOR___T.L.Gulyaeva___Dec. 2011  
C
C INDEX SUBROUTINES:___INDX2007.FOR___T.L.Gulyaeva___June 2009  
C
C ******************* ISO_IRI_SMI PROJECT *******************
C
C ________________________ISOMAIN3.FOR_______________________
C
C UPDATED FOR ISO_IRI CORRECTIONS____
C
C Changes from ISOMAIN2.FOR_:________________________Sep. 2005
C
C New SUBROUTINE TOPH05 is composed for modeling h05top height
C at Ne05=0.5*NmF2 deduced from ISIS and IK19 topside electron
C density profiles. Topside IRI-Bent model is corrected with
C QF-factor obtained from h05top, NmF2 and hmF2.
C Dependence of h05top(COV) introduced
C
C ==    AMENDMENTS :               T.L.Gulyaeva    Sep. 2001
C
C == STORM RELATED SUBROUTINES FROM IRI2001: APF, STORM ==
C    Ref.[IRIFUN.FOR, version 2001.1, May 7, 2001 D. Bilitza]
C
C 2001.1 05/07/01 Include storm subroutine STORM and Ap access s/w
C The subroutine APF uses the input unit 13 for reading 3-hourly
C magnetic Ap indices. The subroutine TCON uses input unit 12 for
C reading the 12-month-running means of ionospheric index IG and
C of sunspot number Rz.
C
C ==  GKPM-subroutine using APF-read ap_indices produces:
C ==     instantaneous kp-index from ap-index
C ==     time-weighted accumulation Kpm index from kp-indices
C ==     for plasmasphere model
C ==  F81=COV81-subroutine produces daily F10.7 (Covington) solar
C ==  radio flux index averaged for 81 days preceding given day
C ==  R81=Rz81-subroutine produces daily sunspot number Rz
C ==  averaged for 81 days preceding given day
C**************************************************************  
C**************************************************************  
C
C********** INTERNATIONAL REFERENCE IONOSPHERE ****************  
C**********RUSSIAN STANDARD MODEL OF IONOSPHERE**************** 
C****************  FUNCTIONS,SUBROUTINES  *********************
C**************************************************************
C** initialize: INITIALIZE (needs to be called before using 
C**                 subroutines or functions)
C** NE:         XE1,ZERO,DXE1N,XE2,XE3_1,XE4_1,XE5,XE6,XE_1
C** TE/TI:      ELTEIK,SPHARM_IK,TEBA,SPHARM,ELTE,TEDE,TI,TEDER,
C** TN,DTNDH
C** NI:         RPID,RDHHE,RDNO,KOEFP1,KOEFP2,KOEFP3,SUFE
C**               IONCO2, APROK,IONCOM_NEW,IONCO1
C** PEAKS:      FOUT,XMOUT,HMF2ED,FOF1ED,f1_c1,f1_prob,FOEEDI,XMDED,
C**  GAMMA1
C** PROFILE PAR:B0_98,TAL,VALGUL
C** MAG. FIELD: GGM,FIELDG,CONVER(corrected latitude)
C** FUNCTIONS:  REGFA1
C** TIME:       SOCO,HPOL,MODA,UT_LT
C** EPSTEIN:    RLAY,D1LAY,D2LAY,EPTR,EPST,EPSTEP,EPLA
C** LAY:        XE2TO5,XEN,VALGUL,ROGUL,LNGLSN,LSKNM,INILAY
C** INDICES:    TCON,TCONGEC, SUBAPF,GKPM
C** Storm:      LSTID,STORM, SUBAPF,GKPM
C** ion drift:  vdrift
C**************************************************************  
C  
C**************************************************************  
C***  -------------------ADDRESSES------------------------  ***
C***  I  PROF. K. RAWER             DR. D. BILITZA       I  ***
C***  I  HERRENSTR. 43              GSFC CODE 933        I  ***
C***  I  7801 MARCH 1               GREENBELT MD 20771   I  ***
C***  I  F.R.G.                     USA                  I  ***
C***  ----------------------------------------------------  ***
C***  I                 DR. T.L. GULYAEVA                I  ***
C***  I                 IZMIRAN 142190                   I  ***
C***  I                 TROITSK, MOSCOW                  I  ***
C***  I                 RUSSIA                           I  ***
C***  I                 gulyaeva@izmiran.ru              I  ***
C***  ----------------------------------------------------  ***
C**************************************************************  
C
C
      SUBROUTINE GKPM(IAS,DKP,RKP)
C      
C    T.L. Gulyaeva, Dec. 2013:
C New: put real kp instead of Kpm calculation:
C
C    T.L. Gulyaeva, Sep. 2001:
C ===================================================================
C calculates geomagnetic 3h kp-index, dkp(13), and Kpm index, rkp,
C     used by the plasmasphere model 
C for given year, month, day and UT hour using ap-indices  
C collected by APF-subroutine in the array ias(13).
C The ap-indices are transformed to kp-indices with conventional 
C thresholds for ap - kp .
C Kpm presents forecast of kp-index 3h in advance based on accumulation 
C of instantaneous kp-indices for 12 preceding hours ranked by 
C decreasing order.
C
C Ref. T.L. Gulyaeva, G. De Franceschi, and L. Perrone. 'Electron
C temperature variations at the F2 layer peak height durung the space
C weather month of September 1999'. Adv. Space Res.31(4), 965-970,2003.
C
      DIMENSION ias(13),lap(0:27),tkp(0:27),dkp(13),X1(4),X2(4)
     &,icurkp(0:7),iprekp(0:7)
      DATA lap/0,2,3,4,5,6,7,9,12,15,18,22,27,32,39,48,56,67,80,94,111
     &,132,154,179,207,236,300,400/
      DATA tkp/0.0,.3,.7,1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3
     &,4.7,5.0,5.3,5.7,6.0,6.3,6.7,7.0,7.3,7.7,8.0,8.3,8.7,9.0/
      COMMON /FIKP/icurkp,iprekp,iyear_cur,imn_cur,idy_cur
C
C calculate 3h kp-indices : kp(13) presents kp-index for current 3h UT bin
C
      DO J=1,13
      IF (IAS(J).EQ.0) THEN
         DKP(J)=0.0
         GOTO 2
      ELSE
            DO I=1,27
              IF ((IAS(j).GT.LAP(i-1)).AND.(IAS(j).LE.LAP(i))) THEN
                   DKP(j)=TKP(i)
                   GOTO 2
               ENDIF
            ENDDO
      ENDIF
   2     ENDDO
C New: put real kp instead of Kpm calculation:
       RKP=dkp(13)
       goto 870
C calculate Kpm-index
       RKP=0.0
       STAU=1.0
C       K1=13
        K1=8
          DO N=1,4
            NK=K1+N
            X1(N)=DKP(NK) 
        X2(N)=0.0
          ENDDO
            MM=0
C rank 4 preceding 3h kp-indices by decreasing order

      DO M=1,4

C find naximum of 4 preceding 3h kp-indices (XM)
          XM=0.0
      DO N=1,4
              IF (XM.LT.X1(N)) THEN 
           XM=X1(N)
           ENDIF
      ENDDO

         IF (XM.EQ.0.0) GOTO 860

C ranking 4 preceding kp by decreasing order
      DO N=1,4
              IF (X1(N).LT.XM) GOTO 830
            MM=MM+1 
        X2(MM)=XM 
        X1(N)=0.0
830       CONTINUE
        ENDDO
      ENDDO
C
860   DO K=1,4
           RKP=RKP+STAU*X2(K)
           STAU=STAU*.7788
      ENDDO
        RKP=RKP*0.35
 870     RETURN
      END
C**************************************************************  
C
C
        subroutine initialize
      dimension indap(13),akp(13)
      COMMON /CONST/UMR,PI/ARGEXP/ARGMAX  
     &       /const1/humr,dumr,xfap,indap,akp,cov,rzs
        ARGMAX=88.0
        pi=atan(1.0)*4.
        UMR=pi/180.
        humr=pi/12.
        dumr = pi / 182.5
        return 
        end
C        
C*************************************************************   
CC 
      SUBROUTINE SUBAPF(IYYYY,IMN,IDY,HOUR,IAP,F81,R81)
CC-------------------------------------------------------------------
CC T.L. Gulyaeva........................................May. 2017 
C
CC Source: former IRI-Plas subroutine: 
CC        SUBROUTINE APF(IYYYY,IMN,ID,HOUR,IAP,F81,R81)
CC Modified to produce sunspot number R81, and solar radio flux F81 
CC from two successive kpYEAR annual files
CC Input: 
CC       kpYEAR  (preceding year: IYYYY-1)                              
CC       kpYEAR  (given year: IYYYY)                              
c--------------------------------------------------------------------
c D. Bilitza .......................................... May, 2001.
c
c finds Ap indices for IRI-storm for given year IYYYY (yyyy), month
c (IMN) and day (ID) and UT hour (HOUR decimal hours). The indices are
c stored in IAP(13). IAP(13) is ap index for UT=hour. Intervals are UT:
c (0-3),(3-6),(6-9),(9-12),(12-15),(15-18),(18-21),(21-24). 
c 
c T.L. Gulyaeva ......................................  Sep. 2001. 
c Modified for producing solar radio flux F10.7, F81, averaged for
c 81 days preceeding given day
C 
C------------------------------------------------------------------------
      DIMENSION iap(13),iap0(8),iap1(8),iap2(8),iapt(24),lm(12)
     &,icv(82),irz(82),icurkp(0:7),iprekp(0:7)
      CHARACTER*4 YEAR,YEAR_pre
      CHARACTER*80 infilekp
      COMMON /FIKP/icurkp,iprekp,iyear_cur,imn_cur,idy_cur
      DATA LM/31,28,31,30,31,30,31,31,30,31,30,31/
      COMMON /path/datapath
      character datapath*200
C
C START: 
      ipcnt=0 ! keep cnt -82,...,-1

      do i=1,8
                iap0(i)=0
           iap1(i)=0
                iap2(i)=0
                enddo
         do k=1,82
       irz(k)=0
       icv(k)=0
         enddo

C   Input file kpYEAR starts from 1948:
      if((iyyyy.lt.1948).or.(iyyyy.gt.iyear_cur)) goto 21 !*
      IF ((iyyyy.eq.iyear_cur).and.(imn.gt.imn_cur)) goto 21 !*
      IF ((iyyyy.eq.iyear_cur).and.(imn.eq.imn_cur).and.    !*
     *   (idy.gt.idy_cur)) goto 21 !*TEST
C
        iy=iyyyy-1900

 770       z1=iy/4.0
        lm(2)=28
           if(iyyyy/4*4.eq.iyyyy) lm(2)=29
      if (iy.ge.100) then
      iyr=iy-100
      else
      iyr=iy
      endif
      call blet4(iyyyy,YEAR)

      iYYYY_pre=iYYYY-1
      call blet4(iYYYY_pre,YEAR_pre)
C
      idend=idy
C
C Day-of-year:
      ldaend=ndoy(iyyyy,imn,idend)       ! current day-of-year
C
      ipre=0
       if (ldaend.lt.81) then
        ipre=1    ! Include data for preceding year
       endif
C
      infilekp='kpYEAR'  
C
10    FORMAT(3I2,6X,8(I2),3X,8(I3),7X,I3,I3,1X,I1)

C
C-------------------------------------------------------------------------------
C Go to read current year file kpYEAR (avoid preceding year input) if current day IDY>81 day-of-year:
C
      if (ipre.eq.0) GOTO 203 !avoid input of preceding year file kpYEAR>>>>>
C
C-------------------------------------------------------------------------------
C Include Input file 'kpYEAR' for preceding year 
C
      infilekp(3:6)=YEAR_pre  !!! Tamara

      OPEN(13,FILE=TRIM(ADJUSTL(DATAPATH))//infilekp,STATUS='OLD',
     & ERR=21)  
  201 continue
                do jj=1,8
                        iap2(jj)=iap1(jj)
                        iap1(jj)=iap0(jj)
                        enddo

      READ(13,10,END=202) JYR,JMN,JDY,(icurkp(k),k=0,7)
     &,(iap0(i),i=1,8),ir1,if1,if2
C
C Replace missed F10.7 by regression model:
      if (if1.eq.0) then
      cov=63.8255+0.8868*ir1  ! daily F10.7(Rz) model
      icv(82)=nint(cov*10.)
      else
      icv(82)=if1*10+if2  !F107*10
      endif
      irz(82)=ir1
C Move data 1 line up:
      do k=1,81
      irz(k)=irz(k+1)
      icv(k)=icv(k+1)
      enddo
      GOTO 201
  202 close(unit=13)
      ipcnt=81-ldaend   ! 
C
C For current year:
  203 infilekp(3:6)=YEAR   !
C
      OPEN(13,FILE=TRIM(ADJUSTL(DATAPATH))//infilekp,ERR=21)
C
C Read file kpyear for current year:
C
  204 continue  
      do jj=1,8  
              iprekp(jj-1)=icurkp(jj-1)
                        iap2(jj)=iap1(jj)
                        iap1(jj)=iap0(jj)
                  enddo

      READ(13,10,err=29,END=205) JYR,JMN,JDY,(icurkp(k),k=0,7)
     &,(iap0(i),i=1,8),ir1,if1,if2
      irz(82)=ir1
      if (if1.eq.0) then
       cov=63.8255+0.8868*ir1  ! daily F10.7(Rz) model
       icv(82)=nint(cov*10.)
      else
       icv(82)=if1*10+if2  !F107*10
      endif
C Move data 1 line up:
  140 do  k=1,81
      irz(k)=irz(k+1)
      icv(k)=icv(k+1)
      enddo
C
      ipcnt=ipcnt+1

      if (ipcnt.lt.81) goto 204

      IF ((JYR.eq.IYR).and.(JMN.eq.IMN).and.(JDY.eq.idend)) then 
      goto 207
      ELSE
      ipcnt=ipcnt-1
      goto 204
      ENDIF

C data for 81 preceding days are collected                       
  207 F81=0.0
      R81=0.0
      icnv=0
      do k=1,81
      R81=R81+irz(k)
      if (icv(k).gt.0) then
      F81=F81+icv(k)
       icnv=icnv+1
      endif
         enddo
      F81=F81/(icnv*10.)
      R81=R81/81.0
C
C === collection of ap-indices for 13  3h UT bins:
 205         do jj=1,8
                iapt(jj)=iap2(jj)
                iapt(8+jj)=iap1(jj)
                iapt(16+jj)=iap0(jj)
                enddo
      if(hour.ge.24.0) hour=0.0
        ih=int(hour/3.)+1
        ihb=16+ih
        ilauf=ihb
        do jj=13,1,-1
                iap(jj)=iapt(ilauf)
                ilauf=ilauf-1
                enddo
C
      goto 30
21      write(*,100)
100     format(1X,'Date is outside range of Kp-Ap indices file.',
     &     ' STORM model is turned off.')
        IAP(1)=-iap(1)
      goto 30
  29  write(*,*) 'Error in input k-index file',infilekp 
C      pause ' '
      stop
  30  CLOSE(unit=13)
      RETURN
      END

C __________________________________________________________________________
C
C
      integer function ndoy(nyear,nmn,ndy)
C To define day-of-year
C
          DIMENSION IM(0:12)
      DATA IM/31,31,28,31,30,31,30,31,31,30,31,30,31/
      IF(nyear/4*4.EQ.nyear) THEN
               IM(2)=29
        ELSE
                IM(2)=28
          ENDIF
C Day-of-year LDA1
      mosum=0
C NEW++++++++++++++++++++
      if (nmn.eq.0) then
      nmn=12
      nyear=nyear-1
      endif
      if(nmn.gt.1) then
         do 1234 i=1,nmn-1
1234    mosum=mosum+im(i)
         endif
      ndoy=mosum+ndy 
      END FUNCTION
C -----------------------------------------------------------
      subroutine blet2(jlet,alet2)
c
      CHARACTER*1 IN(0:9)
C      integer*2 ilet
      DATA IN/'0','1','2','3','4','5','6','7','8','9'/
      CHARACTER*2 alet2
      ilet=jlet
      alet2='00'
      j10=ilet/10
      j1=ilet-j10*10
    5 do i=0,9
      if (j10.eq.i) then
      alet2(1:1)=IN(i)
      endif
      enddo
   7  do i=0,9 
      if (j1.eq.i) then
      alet2(2:2)=IN(i)
      endif
      enddo 
      return
      end
C--------------------------------------------------------------------------------
C
      subroutine blet3(ilet,alet3)
c     nn=1,2,3,4 

      CHARACTER*1 IN(0:9)
      DATA IN/'0','1','2','3','4','5','6','7','8','9'/
      CHARACTER*3 alet3
C
      j100=ilet/100
         j10=(ilet-j100*100)/10
      j1=ilet-j100*100-j10*10
    3 do i=0,9
      if (j100.eq.i) then
      alet3(1:1)=IN(i)
      endif
      enddo
    5 do i=0,9
      if (j10.eq.i) then
      alet3(2:2)=IN(i)
      endif
      enddo
   7  do i=0,9 
      if (j1.eq.i) then
      alet3(3:3)=IN(i)
      endif
      enddo 
      return
      end
C
C---------------------------------------------------------------
      subroutine blet4(ilet,alet4)
c     nn=1,2,3,4 
      INTEGER ilet
      CHARACTER*4 alet4
      CHARACTER*1 IN(0:9)
      DATA IN/'0','1','2','3','4','5','6','7','8','9'/
      alet4='0000'
C
      j1000=ilet/1000
         j100=(ilet-j1000*1000)/100
      j10=(ilet-j1000*1000-j100*100)/10
      j1=ilet-j1000*1000-j100*100-j10*10
    3 do i=0,9
      if (j1000.eq.i) then
      alet4(1:1)=IN(i)
      endif
      enddo
    4 do i=0,9
      if (j100.eq.i) then
      alet4(2:2)=IN(i)
      endif
      enddo
    5 do i=0,9
      if (j10.eq.i) then
      alet4(3:3)=IN(i)
      endif
      enddo
   7  do i=0,9 
      if (j1.eq.i) then
      alet4(4:4)=IN(i)
      endif
      enddo 
      return
      end
C ===============================================================
c
      SUBROUTINE STORM(ap,rga,rgo,coor,rgma,ut,doy,cf,rap)
C  E.A. Araujo-Pradere and T.J. Fuller-Rowell .............May, 2001
C       
C     Refs. 
C     (i) T.J.Fuller-Rowell, M.V.Codrescu, and E.A.Araujo-Pradere,
C         Capturing the storm-time ionospheric response in an empirical
C         model. AGU Geophys. Monograph, 125, 393-401, 2001.
C     (ii)E.A. Araujo-Pradere, T.J. Fuller-Rowell, and M.V. Codrescu,
C         STORM: An empirical storm-time ionospheric correction model,
C         Radio Sci.,37(5),1070, DOI:10.1029/2001RS002467, 2002.
C
C      Fortran code to obtain the foF2 storm-time correction factor at a
C      given location and time, using the current and the 12 previous
C      ap values as input.
C
C      ap ---> (13 elements integer array). Array with the preceeding
C              13 value of the 3-hourly ap index. The 13th value
C              in the array will contain the ap at the UT of interest,
C              the 12th value will contain the 1st three hourly interval
C              preceeding the time of interest, and so on to the 1st
C              ap value at the earliest time.
C     coor --> (integer). If coor = 2, rga should contain the geomagnetic
C                         latitude.
C                         If coor = 1, rga should contain the geographic
C                         latitude.
C     rga ---> (real, -90 to 90) geographic or geomagnetic latitude.
C     rgo ---> (real, 0 to 360, positive east from Greenwich.)
C                           geographic longitude, only used if coor=1.
C     rgma --> (corrected geomagnetic latitude)
C     ut  ---> (integer, hours 00 to 23) Universal Time of interest.
C     doy ---> (integer, 1 to 366)Day of the year.
C     cf  ---> (real) The output; the storm-time correction factor used
C              to scale foF2, foF2 * cf.
C     rap ---> integrated ap-index used for producing storm factor cf

C     DIMENSIONS AND COEFFICIENTS VALUES
cc-   COMMON   
cc-     &         /CONST1/HUMR,DUMR,XFAP,INDAP,AKP,COV,RZS
      DIMENSION c4(20)
      DATA c4/0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,
     +0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,
     +0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00/

      DIMENSION c3(20)
      DATA c3/0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-9.44E-12,
     +0.00E+00,3.04E-12,0.00E+00,9.32E-12,-1.07E-11,0.00E+00,0.00E+00,
     +0.00E+00,1.09E-11,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-1.01E-11/

      DIMENSION c2(20)
      DATA c2/1.16E-08,0.00E+00,0.00E+00,-1.46E-08,0.00E+00,9.86E-08,
     +2.25E-08,-1.67E-08,-1.62E-08,-9.42E-08,1.17E-07,4.32E-08,3.97E-08,
     +3.13E-08,-8.04E-08,3.91E-08,2.58E-08,3.45E-08,4.76E-08,1.13E-07/

      DIMENSION c1(20)
      DATA c1/-9.17E-05,-1.37E-05,0.00E+00,7.14E-05,0.00E+00,-3.21E-04,
     +-1.66E-04,-4.10E-05,1.36E-04,2.29E-04,-3.89E-04,-3.08E-04,
     +-2.81E-04,-1.90E-04,4.76E-05,-2.80E-04,-2.07E-04,-2.91E-04,
     +-3.30E-04,-4.04E-04/

      DIMENSION c0(20)
      DATA c0/1.0136E+00,1.0478E+00,1.00E+00,1.0258E+00,1.00E+00,
     +1.077E+00,1.0543E+00,1.0103E+00,9.9927E-01,9.6876E-01,1.0971E+00,
     +1.0971E+00,1.0777E+00,1.1134E+00,1.0237E+00,1.0703E+00,1.0248E+00,
     +1.0945E+00,1.1622E+00,1.1393E+00/

      DIMENSION fap(36)
      DATA fap/0.,0.,0.037037037,0.074074074,0.111111111,0.148148148,
     10.185185185,0.222222222,0.259259259,0.296296296,0.333333333,
     20.37037037,0.407407407,0.444444444,0.481481481,0.518518519,
     30.555555556,0.592592593,0.62962963,0.666666667,0.703703704,
     40.740740741,0.777777778,0.814814815,0.851851852,0.888888889,
     50.925925926,0.962962963,1.,0.66666667,0.33333334,0.,0.333333,
     60.666666,1.,0.7/

      integer code(8,6)
      data code/3,4,5,4,3,2,1,2,3,2,1,2,3,4,5,4,8,7,6,7,8,9,10,9,
     *13,12,11,12,13,14,15,14,18,17,16,17,18,19,20,19,18,17,16,17,
     *18,19,20,19/

      INTEGER ape(39)
      INTEGER ap(13)
      INTEGER ut,doy,dayno,coor,s1,s2,l1,l2
      REAL rgma, rap, rga, rgo, rs, rl
C
C -----------------------------------------------------------------------
      rgo1=rgo
C      CALLING THE PROGRAM TO CONVERT TO GEOMAGNETIC COORDINATES
      IF (coor .EQ. 1) THEN
          CALL CONVER (rga,rgo,rgma)
      ELSE IF (coor .EQ. 2) THEN
               rgma = rga
      ELSE
         WRITE (6,*)' '
         WRITE (6,*)' '
         WRITE (6,*)'    Wrong Coordinates Selection -------- >>', coor
         WRITE (6,*)' '
         GOTO 100
      ENDIF

      IF (COOR.EQ.2) THEN RGMA=RGA

C     FROM 3-HOURLY TO HOURLY ap
      i = 1
      DO 10 k = 1,13
         DO j = 1,3

      ape(i) = ap(k)

      i = i + 1
      END DO
10     CONTINUE

C     TO OBTAIN THE INTEGRAL OF ap.
C     INTEGRAL OF ap

       if(ut.eq.24) ut=0
      IF (ut .EQ. 0 .OR. ut .EQ. 3 .OR. ut .EQ. 6 .OR. ut .EQ. 9 .OR.
     1ut .EQ. 12 .OR. ut .EQ. 15 .OR. ut .EQ. 18 .OR. ut .EQ. 21) THEN
          k = 1
      ELSE IF (ut .EQ. 1 .OR. ut .EQ. 4 .OR. ut .EQ. 7 .OR. ut .EQ. 10
     1.OR.ut .EQ. 13 .OR. ut .EQ. 16 .OR. ut .EQ. 19 .OR. ut .EQ. 22)
     2THEN
          k = 2
      ELSE IF (ut .EQ. 2 .OR. ut .EQ. 5 .OR. ut .EQ. 8 .OR. ut .EQ. 11
     1.OR. ut .EQ. 14 .OR. ut .EQ. 17 .OR. ut .EQ. 20 .OR. ut .EQ. 23)
     2THEN
          k = 3

      ELSE

          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)'  Wrong Universal Time value -------- >>', ut
          WRITE (6,*)' '
          GOTO 100

      END IF

      rap = 0

      DO j = 1,36
      rap = rap + fap(j) * ape(k+j)
      END DO

      if(rap.le.200.)then
      cf=1.0
      goto 100
      end if

C Test TLG:
111   continue
      if(doy.gt.366.or.doy.lt.1)then
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)'      Wrong Day of Year value --- >>', doy
          WRITE (6,*)' '
          GOTO 100
      end if

      if(rgma.gt.90.0.or.rgma.lt.-90.0)then
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)' '
          WRITE (6,*)'   Wrong GEOMAGNETIC LATITUDE value --- >>', rgma
          WRITE (6,*)' '
          GOTO 100
      end if

c      write(6,*)rgma

      dayno=doy
      if(rgma.lt.0.0)then
      dayno=doy+172
      if(dayno.gt.365)dayno=dayno-365
      end if

      if (dayno.ge.82) rs=(dayno-82.)/45.6+1.
      if (dayno.lt.82) rs=(dayno+283.)/45.6+1.
      s1=rs
      facs=rs-s1
      s2=s1+1
      if(s2.eq.9) s2=1
c      write(6,*)s1,s2,rs

      rgma = abs(rgma)

      rl=(rgma+10.)/20.+1
      if(rl.eq.6.0)rl=5.9
      l1=rl
      facl=rl-l1
      l2=l1+1
c      write(6,*)l1,l2,rl

C     FACTORS CALCULATIONS

      if(rap.lt.300.)then
      rapf=300.
      n1=code(s1,l1)
      cf1=c4(n1)*(rapf**4)+c3(n1) * (rapf**3) + c2(n1) * (rapf**2) +
     1c1(n1) * rapf + c0(n1)
      n2=code(s1,l2)
      cf2=c4(n2)*(rapf**4)+c3(n2) * (rapf**3) + c2(n2) * (rapf**2) +
     1c1(n2) * rapf + c0(n2)
      n3=code(s2,l1)
      cf3=c4(n3)*(rapf**4)+c3(n3) * (rapf**3) + c2(n3) * (rapf**2) +
     1c1(n3) * rapf + c0(n3)
      n4=code(s2,l2)
      cf4=c4(n4)*(rapf**4)+c3(n4) * (rapf**3) + c2(n4) * (rapf**2) +
     1c1(n4) * rapf + c0(n4)

C     INTERPOLATION

      cf300=cf1*(1 - facs) * (1 - facl) + cf2 * (1 - facs) * (facl) +
     *cf3 * (facs) * (1 - facl) + cf4 * (facs) * (facl)

      cf = (cf300-1.0)*rap/100.-2.*cf300+3.
      goto 100
      end if

      n1=code(s1,l1)
c      write(6,*)n1
      cf1 = c4(n1) * (rap**4) + c3(n1) * (rap**3) + c2(n1) * (rap**2) +
     1c1(n1) * rap + c0(n1)
      n2=code(s1,l2)
      cf2 = c4(n2) * (rap**4) + c3(n2) * (rap**3) + c2(n2) * (rap**2) +
     1c1(n2) * rap + c0(n2)
      n3=code(s2,l1)
      cf3 = c4(n3) * (rap**4) + c3(n3) * (rap**3) + c2(n3) * (rap**2) +
     1c1(n3) * rap + c0(n3)
      n4=code(s2,l2)
      cf4 = c4(n4) * (rap**4) + c3(n4) * (rap**3) + c2(n4) * (rap**2) +
     1c1(n4) * rap + c0(n4)

c     INTERPOLATION

      cf = cf1 * (1 - facs) * (1 - facl) + cf2 * (1 - facs) * (facl) +
     *cf3 * (facs) * (1 - facl) + cf4 * (facs) * (facl)
100   CONTINUE

      RETURN

      END


C T.L. Gulyaeva ............................................ ........ 2007
C
C New subroutine TOP05 to correct Topside IRI-Bent model using 
C QF-factor obtained from h05top height at Ne=0.5*NmF2, NmF2 and hmF2  .
C
CREM        SUBROUTINE TOPH05(covi,AMLAT,TIME,HMAX,HT05,sg)
         SUBROUTINE TOPH05(rz,AMLAT,TIME,HMAX,HT05,sg)
C
C
C//  Gulyaeva T.L. (2003) Variations in the half-width of the topside ionosphere 
C//     according to the observations by space ionosondes ISIS 1,ISIS 2, and IK19.
C//     International J. of Geomagnetism and Aeronomy, 4(3), 201-207.
C//  Gulyaeva T.L., Titheridge J.E. (2006) Advanced specification of electron density 
C//     and temperature in the IRI ionosphere-plasmasphere model. 
C//     Adv. Space Res. 38(11), 2587-2595, doi:10.1016/j.asr.2005.08.045.
C
C  Implementation of empirical RAT=(h05top-hmF2)/hmF2 derived from ISIS and IK19
C  topside electron density profiles to obtain half peak density topside height
C  h05top  from the Chebishev polinomial coefficients given for 
C  (1) 4 levels of solar activity: Rz= 0,  50, 100, 150 
C      solar radio flux          covi=60, 106, 152, 198 (F10.7 option)
C  (2) 10 selected grids of geomagnetic latitude (N=S):0,10,20,30,40,50,60,70,80,90
C  (3) 5 selected grids of local time: 0, 6, 12, 18, 24.
C  (4) 4 seasonal grids: 1 equinox(sg=90deg), 2 summer (sg=180), 
C                        3 equinox (sg=270), 4 winter(sg=360)
C
      DIMENSION CVLEV(4)   
      COMMON     /BLOCK1/HMF2,XNMF2,XHMF1,F1REG         
     &          /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU,HTOP,QF0
crem  DATA CVLEV/60.,106.,152.,198./
      DATA CVLEV/0.,50.,100.,150./      ! Rz levels
      LOGICAL F1REG

      ABMLAT=ABS(AMLAT)
crem      IR=IFIX((covi-60.)/46.)+1 
      IR=IFIX(rz/50.)+1 
      M1=IFIX(ABMLAT/10.)+1
         L1=IFIX(TIME/6.)+1
       M2=M1+1
          IF(M1.EQ.10) M2=10
       L2=L1+1
          IF(L1.EQ.5) L2=5
C
C INTERPOLATE RAT FOR GIVEN RZI
C
C Call Chebishev approximation to interpolate for given ABMLAT, HRLT
      call CHEBISH(CVLEV(IR),TIME,ABMLAT,XX,SG)
        IF (IR.EQ.4) THEN
             RAT05=XX
            GOTO 10
      ENDIF
      call CHEBISH(CVLEV(IR+1),TIME,ABMLAT,YY,sg)
crem  RAT05=XX+(YY-XX)*(COVI-CVLEV(IR))/46.
      RAT05=XX+(YY-XX)*(RZ-CVLEV(IR))/50.
   10 HT05=HMAX*(1.+RAT05)
      RETURN
      END
C********************************************************************
      SUBROUTINE CHEBISH(COVS,HOURLT,ABMLAT,RATCH,SG)
C CHEBISHEV POLINOMIALS FOR ABMLAT(10),HOURLT(5)
C CR((C0...C5),(LT=0,6,...24),(SG=season grids=90,180,270,360)
C                    (COV=60,106,152,198)
C Ref. T.L.Gulyaeva, J.E.Titheridge. Advanced specification of electron density and
C temperature in the IRI ionosphere-plasmasphere model. Adv. Space Res., 2005. 

c      REAL UK(0:10),CR(0:5,5,3,4),YI(5),YY(5,3)
      REAL BR(6,5,3,4),YI(5),YY(5,3)
      REAL PL1(5),PL2(5),PL3(5),CL(0:3)
C  
      DATA rad/0.01745329/
      DATA PL1/-2.,-1.,0.,1.,2./
      DATA PL2/2.,-1.,-2.,-1.,2./
      DATA PL3/-1.,2.,0.,-2.,1./
      DATA BR/
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=60 (mlat/10=0,1,...,9)
C Equinox   B0MLAT:
     *  -1.5859,  3.5789, -3.7884, 2.7094,-1.2962,.6759
     *, -5.3449, 12.8379,-12.0165, 5.9746,-1.9084,.7669
     *,-12.8000, 35.3084,-38.0043,19.6004,-4.4974,.6975
     *,  5.8282,-13.3538,  9.1674,-0.9593,-0.8909,.6062
     *, -1.5859,  3.5789, -3.7884, 2.7094,-1.2962,.6759
C Summer B0MLAT    
     *, -7.1103, 21.0389,-24.5539,14.1607,-3.8537,.7266
     *,  5.5333,-10.6242,  4.8751, 1.0587,-1.0821,.7527
     *,-15.4487, 42.9269,-45.0314,21.4718,-4.2116,.6026
     *, -6.6436, 16.4533,-15.5142, 6.8287,-1.2871,.4976
     *, -7.1103, 21.0389,-24.5539,14.1607,-3.8537,.7266
C Winter   B0MLAT                                                                       
     *, 14.9103,-35.2337, 27.3078,-6.5362,-0.6265,.7509
     *,  2.3846, -5.8840,  3.7023, 0.8525,-1.2663,.7086
     *, -9.8846, 26.6649,-27.0173,12.6959,-2.6536,.6295
     *,  1.7692, -2.3578, -0.7945, 2.2477,-0.9691,.5719
     *, 14.9103,-35.2337, 27.3078,-6.5362,-0.6265,.7509
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=106 (mlat=0,10,...,90)
C Equinox   B1MLAT
     *, -4.1218, 10.6136,-11.4922, 6.0470,-1.3620,.5563
     *,  0.9077,  2.9562, -8.9880, 6.8680,-1.9621,.7737
     *,-16.2744, 42.8047,-43.7009,20.7965,-4.0697,.6619
     *,-17.3038, 44.3336,-40.9249,15.9042,-2.1554,.4796
     *, -4.1218, 10.6136,-11.4922, 6.0470,-1.3620,.5563
C Summer B1MLAT  
     *, -4.9692, 16.5753,-21.3543,12.7061,-3.1758,.6446
     *,  1.9000, -2.8167, -0.9962, 3.0687,-1.3454,.6859
     *,  7.6769,-14.8343,  6.7030, 1.5578,-1.0626,.4291
     *,  5.4833,-10.6322,  4.7571, 1.2178,-0.8223,.4615
     *, -4.9692, 16.5753,-21.3543,12.7061,-3.1758,.6446
C Winter B1MLAT  
     *, -4.7282, 13.4491,-15.6931, 8.8388,-1.9732,.5874
     *,  5.6756,-14.8458, 11.8927,-2.2632,-0.6122,.6948
     *,-14.2872, 40.0829,-41.2716,18.1696,-2.7203,.4916
     *,-13.6128, 33.4657,-29.7231,11.0972,-1.2884,.5034
     *, -4.7282, 13.4491,-15.6931, 8.8388,-1.9732,.5874
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=152 (mlat=0,10,...,90)
C Equinox   B2MLAT
     *, -3.3282, 10.4296,-12.4722, 6.7623,-1.5172,.4931
     *, -8.9744, 20.1311,-17.4552, 7.6518,-1.7371,.6702
     *, 12.0462,-27.8932, 20.6241,-4.5781, 0.0814,.3501
     *,-17.0551, 42.3258,-37.1874,13.3608,-1.4804,.4216
     *, -3.3282, 10.4296,-12.4722, 6.7623,-1.5172,.4931
C Summer B2MLAT  
     *,  7.3077,-17.1579, 11.6872,-0.7405,-1.0298,.5754
     *, 19.2641,-45.1886, 34.3297,-8.1879,-0.1875,.6562
     *,  6.0987,-11.0903,  4.3569, 1.4001,-0.7309,.3885
     *,  5.9295,-13.9205, 10.2347,-2.2818, 0.0853,.3915
     *,  7.3077,-17.1579, 11.6872,-0.7405,-1.0298,.5754
C Winter B2MLAT  
     *, -1.6821,  8.6010,-13.6570, 8.6307,-1.9846,.5635
     *,  5.4679,-12.3750,  7.5620, 0.5394,-1.4415,.6659
     *, -8.0821, 21.9288,-21.8597, 9.3455,-1.4644,.3599
     *, -8.3000, 19.3076,-16.3295, 6.1619,-0.9144,.3846
     *, -1.6821,  8.6010,-13.6570, 8.6307,-1.9846,.5635
C Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=198 (mlat=0,10,...,90)
C Equinox   B3MLAT
     *,-16.4051, 28.2880,-16.0982, 4.6328,-1.0405,.5486
     *, 13.0846,-34.8291, 30.0074,-8.6402, 0.1529,.6165
     *, 19.7474,-42.7116, 28.9430,-6.0487, 0.1492,.3748
     *, 16.2795,-36.6982, 26.5094,-6.3492, 0.2926,.3946
     *,-16.4051, 28.2880,-16.0982, 4.6328,-1.0405,.5486
C Summer B3MLAT
     *,  4.6410,-13.7931, 11.6548,-1.9248,-0.7246,.5264
     *, -2.4090,  3.1805, -2.8423, 2.8861,-0.9937,.5234
     *,  6.3410,-13.9643,  8.2461,-0.0186,-0.7009,.3582
     *,  9.0987,-20.8618, 14.7262,-2.8798,-0.0512,.3662
     *,  4.6410,-13.7931, 11.6548,-1.9248,-0.7246,.5264
C Winter B3MLAT
     *, -4.6526, 12.1878,-14.4047, 8.5226,-2.0493,.5903
     *,  3.9821, -6.9477,  0.8382, 3.4898,-1.5694,.6283
     *, -7.0474, 17.3974,-17.3465, 8.3671,-1.5708,.3759
     *,  4.2782, -9.9880,  5.9834, 0.0975,-0.4900,.3842
     *, -4.6526, 12.1878,-14.4047, 8.5226,-2.0493,.5903/

C  DATA UL/-2.,-1.,0.,1.,2./

      do k=0,3
      cl(k)=0.
      enddo
C
c,,   IR=IFIX((covs-60.)/46.)+1  
      IR=IFIX(covs/50.)+1  
C Given geomagnetic latitude parameter:
      xi=abmlat/100.
      DO LS=1,3
          DO LL=1,5
      B1=BR(6,LL,LS,IR)      
      B2=BR(5,LL,LS,IR)        
      B3=BR(4,LL,LS,IR)       
      B4=BR(3,LL,LS,IR)        
      B5=BR(2,LL,LS,IR)        
      B6=BR(1,LL,LS,IR)        
       HLT=(LL-1)*6.0

      YY(LL,LS)=B1+xi*(B2+xi*(B3+xi*(B4+xi*(B5+xi*B6))))
      ENDDO
      ENDDO            ! end of season/day cycle
C Apply seasonal interpolation
       do i=1,5
      p0=(2.*YY(i,1)+YY(i,2)+YY(i,3))/4.
      p1=(YY(i,3)-YY(i,2))/2.
      p2=(YY(i,2)+YY(i,3)-2.*YY(i,1))/4.
      YI(i)=p0+p1*cos(sg*rad)+p2*cos(2.*sg*rad)
       enddo
      DO K=1,5
      CL(0)=CL(0)+YI(K)
      CL(1)=CL(1)+YI(K)*PL1(K)
      CL(2)=CL(2)+YI(K)*PL2(K)
      CL(3)=CL(3)+YI(K)*PL3(K)
      ENDDO
      CL(0)=CL(0)/5.
      CL(1)=CL(1)/10.
      CL(2)=CL(2)/14.
      CL(3)=CL(3)/12.
      ULL=(HOURLT-12.)/6.
      ZA=CL(0)-2.*CL(2)
      RATCH=ZA+ULL*(CL(1)-3.4*CL(3)+ULL*(CL(2)+ULL*CL(3)))

      RETURN
      END   
C 
      real function fmlt(ut,xlat,xlong,day)
C Calculats geomagnetic local time in hours for given
C universal time UT(sec)
C XLAT/XLONG=geodetic Latitude/Longitude (deg.)
C DAY= day of year
C D. Bilitza, COSPAR-1980, paper 7.3.1
C
      rad=57.29578
      flat=xlat/rad
      flong=xlong/rad
      delta=0.409207*sin((day-80.0)*3.14159/184.0)
      beta=(180.-ut/240.0)/rad
      px=cos(flat)*cos(flong)
      py=cos(flat)*sin(flong)
      pz=sin(flat)
      sx=cos(delta)*cos(beta)
      sy=cos(delta)*sin(beta)
      sz=sin(delta)
      p1=0.35117*px-0.91483*py-0.19937*pz
      p2=0.93358*px+0.35837*py
      s1=0.35117*sx-0.91483*sy-0.19937*sz
      s2=0.93358*sx+0.35837*sy
      thetp=atan2(p2,p1)
      thets=atan2(s2,s1)
      fmlt=(thetp-thets+3.141593)/.2617994
      afmlt=anint(fmlt*10.)
      fmlt=afmlt/10.
      if(fmlt.ge.24.) fmlt=fmlt-24.
      if(fmlt.lt.0.) fmlt=fmlt+24.
      return
      end
C
C ---------------------------------------------------------------------
C
      SUBROUTINE CONVER(rga,rgo,rgma)

C     This subroutine converts a geographic latitude and longitude
C     location to a corrected geomagnetic latitude.
C
C     INPUT: 
C       geographic latitude   -90. to +90.
C       geographic longitude  0. to 360. positive east from Greenwich.
C
C     OUTPUT:
C       corrected geomagnetic latitude -90. to +90.


      DIMENSION CORMAG(20,91)      
      DATA ((CORMAG(i,j),i=1,20),j=1,31)/
     +163.68,163.68,163.68,163.68,163.68,163.68,
     +163.68,163.68,163.68,163.68,163.68,163.68,163.68,163.68,
     +163.68,163.68,163.68,163.68,163.68,163.68,162.60,163.12,
     +163.64,164.18,164.54,164.90,165.16,165.66,166.00,165.86,
     +165.20,164.38,163.66,162.94,162.42,162.00,161.70,161.70,
     +161.80,162.14,161.20,162.18,163.26,164.44,165.62,166.60,
     +167.42,167.80,167.38,166.82,166.00,164.66,163.26,162.16,
     +161.18,160.40,159.94,159.80,159.98,160.44,159.80,161.14,
     +162.70,164.50,166.26,167.90,169.18,169.72,169.36,168.24,
     +166.70,164.80,162.90,161.18,159.74,158.60,157.94,157.80,
     +157.98,158.72,158.40,160.10,162.02,164.28,166.64,169.00,
     +170.80,171.72,171.06,169.46,167.10,164.64,162.18,160.02,
     +158.20,156.80,156.04,155.80,156.16,157.02,157.00,158.96,
     +161.24,163.86,166.72,169.80,172.42,173.72,172.82,170.34,
     +167.30,164.22,161.34,158.74,156.60,155.00,154.08,153.90,
     +154.36,155.36,155.50,157.72,160.36,163.32,166.60,170.20,
     +173.70,175.64,174.18,170.80,167.10,163.56,160.24,157.36,
     +154.96,153.10,152.08,151.92,152.46,153.76,154.10,156.52,
     +159.36,162.52,166.24,170.30,174.62,177.48,175.04,170.82,
     +166.60,162.70,159.02,155.88,153.22,151.20,150.08,149.92,
     +150.64,152.20,152.80,155.32,158.28,161.70,165.58,170.00,
     +174.84,178.46,175.18,170.38,165.80,161.64,157.80,154.38,
     +151.52,149.30,148.18,148.02,148.92,150.60,151.40,154.08,
     +157.18,160.68,164.78,169.40,174.34,177.44,174.28,169.44,
     +164.70,160.34,156.30,152.78,149.72,147.40,146.18,146.04,
     +147.12,149.04,150.10,152.88,156.00,159.58,163.78,168.50,
     +173.28,175.60,172.86,168.14,163.40,158.98,154.88,151.10,
     +147.98,145.50,144.18,144.14,145.40,147.48,148.80,151.68,
     +154.88,158.48,162.68,167.40,171.76,173.60,171.12,166.68,
     +162.00,157.48,153.28,149.50,146.18,143.50,142.18,142.24,
     +143.68,145.98,147.50,150.54,153.68,157.28,161.42,166.10,
     +170.10,171.48,169.22,164.98,160.40,155.88,151.68,147.80,
     +144.34,141.60,140.18,140.26,141.98,144.62,146.30,149.34,
     +152.48,155.98,160.08,164.60,168.34,169.38,167.20,163.18,
     +158.60,154.18,149.98,146.02,142.54,139.70,138.18,138.46,
     +140.26,143.16,145.10,148.14,151.18,154.60,158.68,163.10,
     +166.48,167.28,165.18,161.32,156.90,152.48,148.28,144.32,
     +140.74,137.80,136.22,136.48,138.64,141.76,143.90,146.98,
     +149.98,153.30,157.24,161.40,164.52,165.16,162.86,159.42,
     +155.00,150.68,146.48,142.52,138.94,135.90,134.22,134.68,
     +137.02,140.40,142.70,145.84,148.76,151.92,155.74,159.70,
     +162.52,162.96,160.98,157.42,153.10,148.84,144.68,140.82,
     +137.20,134.00,132.32,132.80,135.42,139.10,141.60,144.74,
     +147.46,150.52,154.20,158.00,160.46,160.76,158.86,155.36,
     +151.20,146.94,142.88,139.02,135.40,132.10,130.32,131.00,
     +133.80,137.74,140.50,143.58,146.24,149.12,152.60,156.20,
     +158.40,158.66,156.76,153.36,149.30,145.04,141.08,137.30,
     +133.60,130.30,128.42,129.12,132.28,136.44,139.30,142.48,
     +144.94,147.64,150.48,154.30,156.34,156.36,154.56,151.26,
     +147.30,143.14,139.20,135.50,131.90,128.40,126.52,127.32,
     +130.76,135.18,138.20,141.28,143.72,146.24,149.26,152.40,
     +154.24,154.16,152.36,149.16,145.30,141.24,137.30,133.70,
     +130.10,126.60,124.62,125.54,129.16,133.92,137.10,140.18,
     +142.42,144.66,147.62,150.50,152.18,151.96,150.16,147.10,
     +143.30,139.24,135.50,131.90,128.36,124.80,122.72,123.74,
     +127.64,132.62,135.90,139.02,141.12,143.18,145.92,148.60,
     +149.98,149.76,148.04,145.00,141.20,137.30,133.60,130.10,
     +126.60,123.00,120.86,121.96,126.12,131.36,134.80,137.88,
     +139.80,141.68,144.08,146.60,147.88,147.56,145.84,142.90,
     +139.20,135.30,131.70,128.28,124.86,121.30,118.96,120.18,
     +124.70,130.16,133.60,136.72,138.48,140.10,142.38,144.60,
     +145.72,145.34,143.64,140.80,137.10,133.30,129.72,126.48,
     +123.10,119.50,117.16,118.48,123.18,128.86,132.40,135.42,
     +137.08,138.50,140.54,142.60,143.52,143.06,141.44,138.70,
     +135.10,131.30,127.82,124.58,121.40,117.70,115.26,116.70,
     +121.66,127.60,131.20,134.22,135.66,136.82,138.70,140.60,
     +141.36,140.86,139.24,136.50,133.00,129.30,125.92,122.78,
     +119.60,116.00,113.40,114.92,120.16,126.30,130.00,132.92,
     +134.24,135.14,136.80,138.60,139.16,138.64,137.12,134.40,
     +130.90,127.20,123.92,120.96,117.90,114.20,111.56,113.12,
     +118.64,124.90,128.70,131.56,132.74,133.44,134.90,136.50,
     +137.00,136.36,134.82,132.30,128.70,125.16,121.94,119.06,
     +116.10,112.50,109.70,111.42,117.14,123.60,127.30,130.16,
     +131.22,131.66,133.00,134.50,134.80,134.14,132.62,130.14,
     +126.60,123.06,119.94,117.16,114.30,110.70,107.80,109.64,
     +115.62,122.24,125.90,128.76,129.62,129.96,131.06,132.40,
     +132.60,131.86,130.42,128.00,124.50,120.96,117.96,115.26,
     +112.54,108.90,105.94,107.86,114.02,120.84/

      DATA ((CORMAG(i,j),i=1,20),j=32,61)/
     +124.05,126.79,
     +127.55,127.83,128.90,130.21,130.41,129.71,128.33,125.96,
     +122.49,118.96,115.97,113.26,110.52,106.89,104.01,106.00,
     +112.21,119.06,122.19,124.82,125.48,125.69,126.73,128.03,
     +128.22,127.55,126.23,123.92,120.47,116.97,113.97,111.26,
     +108.50,104.89,102.08,104.14,110.41,117.29,120.34,122.85,
     +123.40,123.56,124.57,125.84,126.03,125.40,124.14,121.88,
     +118.46,114.97,111.98,109.26,106.48,102.88,100.15,102.28,
     +108.60,115.51,118.49,120.88,121.33,121.42,122.40,123.65,
     +123.84,123.24,122.04,119.83,116.45,112.97,109.98,107.26,
     +104.46,100.87,098.22,100.42,106.79,113.74,116.63,118.91,
     +119.26,119.29,120.24,121.47,121.65,121.09,119.95,117.79,
     +114.43,110.98,107.99,105.26,102.44,098.87,096.29,098.56,
     +104.98,111.96,114.78,116.94,117.19,117.15,118.07,119.28,
     +119.46,118.93,117.86,115.75,112.42,108.98,106.00,103.26,
     +100.42,096.86,094.36,096.70,103.18,110.19,112.93,114.97,
     +115.12,115.02,115.91,117.09,117.27,116.78,115.76,113.71,
     +110.41,106.98,104.00,101.26,098.40,094.85,092.43,094.84,
     +101.37,108.41,111.07,113.00,113.04,112.88,113.74,114.91,
     +115.08,114.62,113.67,111.67,108.39,104.99,102.01,099.26,
     +096.38,092.85,090.51,092.97,099.56,106.64,109.22,111.03,
     +110.97,110.75,111.58,112.72,112.89,112.47,111.57,109.63,
     +106.38,102.99,100.01,097.26,094.36,090.84,088.58,091.11,
     +097.75,104.86,107.37,109.06,108.90,108.61,109.41,110.53,
     +110.70,110.31,109.48,107.59,104.37,100.99,098.02,095.26,
     +092.34,088.83,086.65,089.25,095.95,103.09,105.51,107.09,
     +106.83,106.48,107.25,108.35,108.51,108.16,107.39,105.55,
     +102.35,099.00,096.03,093.26,090.32,086.83,084.72,087.39,
     +094.14,101.31,103.66,105.12,104.76,104.34,105.08,106.16,
     +106.32,106.00,105.29,103.50,100.34,097.00,094.03,091.26,
     +088.30,084.82,082.79,085.53,092.33,099.54,101.81,103.15,
     +102.68,102.21,102.92,103.97,104.13,103.85,103.20,101.46,
     +098.33,095.00,092.04,089.26,086.28,082.81,080.86,083.67,
     +090.52,097.76,099.95,101.18,100.61,100.07,100.75,101.79,
     +101.94,101.69,101.10,099.42,096.31,093.01,090.04,087.26,
     +084.26,080.81,078.93,081.81,088.72,095.99,098.10,099.21,
     +098.54,097.94,098.59,099.60,099.75,099.54,099.01,097.38,
     +094.30,091.01,088.05,085.26,082.24,078.80,077.00,079.95,
     +086.91,094.21,096.25,097.24,096.47,095.81,096.43,097.41,
     +097.56,097.39,096.92,095.34,092.29,089.01,086.06,083.26,
     +080.22,076.79,075.07,078.09,085.10,092.43,094.39,095.27,
     +094.40,093.67,094.26,095.23,095.37,095.23,094.82,093.30,
     +090.27,087.02,084.06,081.26,078.20,074.79,073.14,076.23,
     +083.30,090.66,092.54,093.30,092.32,091.54,092.10,093.04,
     +093.18,093.08,092.73,091.26,088.26,085.02,082.07,079.26,
     +076.18,072.78,071.21,074.37,081.49,088.88,090.69,091.33,
     +090.25,089.40,089.93,090.85,090.99,090.92,090.63,089.21,
     +086.25,083.02,080.07,077.26,074.16,070.77,069.28,072.51,
     +079.68,087.11,088.83,089.36,088.18,087.27,087.77,088.67,
     +088.80,088.77,088.54,087.17,084.23,081.03,078.08,075.26,
     +072.14,068.77,067.35,070.65,077.87,085.33,086.98,087.39,
     +086.11,085.13,085.60,086.48,086.61,086.61,086.45,085.13,
     +082.22,079.03,076.09,073.26,070.12,066.76,065.42,068.79,
     +076.07,083.56,085.13,085.42,084.04,083.00,083.44,084.29,
     +084.42,084.46,084.35,083.09,080.21,077.03,074.09,071.26,
     +068.10,064.75,063.49,066.93,074.26,081.78,083.27,083.45,
     +081.96,080.86,081.27,082.11,082.23,082.30,082.26,081.05,
     +078.19,075.04,072.10,069.26,066.08,062.75,061.57,065.06,
     +072.45,080.01,081.42,081.48,079.89,078.73,079.11,079.92,
     +080.04,080.15,080.16,079.01,076.18,073.04,070.10,067.26,
     +064.06,060.74,059.64,063.20,070.64,078.23,079.57,079.51,
     +077.82,076.59,076.94,077.73,077.85,077.99,078.07,076.97,
     +074.17,071.04,068.11,065.26,062.04,058.73,057.71,061.34,
     +068.84,076.46,077.71,077.54,075.75,074.46,074.78,075.55,
     +075.66,075.84,075.98,074.93,072.15,069.05,066.12,063.26,
     +060.02,056.73,055.78,059.48,067.03,074.68,075.86,075.57,
     +073.68,072.32,072.61,073.36,073.47,073.68,073.88,072.88,
     +070.14,067.05,064.12,061.26,058.00,054.72,053.85,057.62,
     +065.22,072.91,074.01,073.60,071.60,070.19,070.45,071.17,
     +071.28,071.53,071.79,070.84,068.13,065.05,062.13,059.26,
     +055.98,052.71,051.92,055.76,063.41,071.13,072.15,071.63,
     +069.53,068.05,068.28,068.99,069.09,069.37,069.69,068.80,
     +066.11,063.06,060.13,057.26,053.96,050.71,049.99,053.90,
     +061.61,069.36,070.30,069.66,067.46,065.92,066.12,066.80,
     +066.90,067.22,067.60,066.76,064.10,061.06,058.14,055.26,
     +051.94,048.70,048.06,052.04,059.80,067.58/

      DATA ((CORMAG(i,j),i=1,20),j=62,91)/
     +067.70,067.06,
     +065.08,063.72,063.98,064.60,064.80,065.12,065.60,064.86,
     +062.40,059.26,056.24,053.18,049.84,046.60,046.12,050.12,
     +057.52,064.80,064.90,064.42,062.70,061.62,061.78,062.40,
     +062.60,063.04,063.58,063.00,060.60,057.46,054.42,051.18,
     +047.70,044.60,044.22,048.02,055.06,061.92,062.10,061.72,
     +060.32,059.50,059.68,060.20,060.46,060.94,061.58,061.00,
     +058.70,055.66,052.52,049.18,045.60,042.50,042.22,046.00,
     +052.60,058.98,059.20,059.18,058.12,057.32,057.48,058.00,
     +058.30,058.84,059.48,059.04,056.90,053.86,050.62,047.10,
     +043.50,040.50,040.28,043.98,050.22,056.18,056.40,056.64,
     +055.84,055.20,055.38,055.80,056.16,056.84,057.48,057.04,
     +055.10,052.06,048.70,045.10,041.40,038.40,038.28,041.88,
     +047.94,053.44,053.70,054.14,053.56,053.10,053.24,053.70,
     +054.06,054.74,055.38,055.14,053.20,050.26,046.80,043.10,
     +039.34,036.40,036.38,039.96,045.56,050.84,051.10,051.70,
     +051.36,051.00,051.14,051.50,051.96,052.64,053.38,053.08,
     +051.30,048.36,044.90,041.02,037.24,034.40,034.38,037.86,
     +043.28,048.20,048.50,049.26,049.18,048.90,049.04,049.40,
     +049.86,050.64,051.28,051.08,049.40,046.46,042.98,039.02,
     +035.14,032.40,032.48,035.72,041.00,045.70,046.00,046.96,
     +046.98,046.80,046.94,047.30,047.76,048.54,049.28,049.08,
     +047.40,044.56,041.08,037.02,033.14,030.40,030.58,033.84,
     +038.72,043.20,043.50,044.62,044.80,044.80,044.94,045.20,
     +045.76,046.54,047.18,046.98,045.50,042.66,039.08,035.02,
     +031.14,028.40,028.58,031.82,036.52,040.80,041.20,042.32,
     +042.54,042.70,042.84,043.20,043.66,044.44,045.08,044.98,
     +043.50,040.76,037.08,033.04,029.04,026.40,026.68,029.82,
     +034.34,038.40,038.80,040.12,040.60,040.70,040.84,041.10,
     +041.62,042.34,042.98,042.88,041.50,038.76,035.18,031.04,
     +027.14,024.50,024.78,027.70,032.14,036.06,036.50,037.88,
     +038.50,038.68,038.84,039.10,039.56,040.34,040.88,040.82,
     +039.40,036.76,033.18,029.12,025.14,022.50,022.88,025.90,
     +029.96,033.86,034.30,035.68,036.42,036.68,036.84,037.10,
     +037.56,038.24,038.88,038.72,037.40,034.76,031.18,027.12,
     +023.14,020.60,020.98,023.90,027.88,031.66,032.10,033.58,
     +034.32,034.68,034.84,035.10,035.56,036.24,036.78,036.62,
     +035.30,032.72,029.18,025.14,021.24,018.70,019.08,021.90,
     +025.88,029.42,029.90,031.48,032.32,032.68,032.84,033.10,
     +033.56,034.22,034.68,034.42,033.20,030.72,027.28,023.22,
     +019.34,016.80,017.24,020.00,023.78,027.32,027.70,029.38,
     +030.24,030.68,030.94,031.20,031.66,032.22,032.58,032.32,
     +031.10,028.62,025.28,021.32,017.48,015.00,015.38,018.18,
     +021.80,025.22,025.70,027.28,028.24,028.78,029.04,029.30,
     +029.66,030.22,030.50,030.22,029.00,026.62,023.30,019.42,
     +015.64,013.10,013.54,016.28,019.80,023.12,023.60,025.24,
     +026.24,026.78,027.14,027.40,027.76,028.22,028.40,028.12,
     +026.80,024.52,021.30,017.52,013.78,011.30,011.74,014.48,
     +017.90,021.12,021.60,023.24,024.34,024.88,025.24,025.50,
     +025.86,026.22,026.40,025.98,024.70,022.48,019.40,015.72,
     +012.04,009.50,009.94,012.58,016.02,019.12,019.60,021.24,
     +022.34,022.98,023.34,023.70,024.00,024.30,024.40,023.88,
     +022.60,020.48,017.52,014.00,010.34,007.80,008.18,010.88,
     +014.22,017.18,017.60,019.34,020.44,021.16,021.54,021.90,
     +022.16,022.40,022.32,021.78,020.60,018.48,015.62,012.20,
     +008.68,006.00,006.44,009.18,012.42,015.28,015.80,017.44,
     +018.54,019.26,019.74,020.10,020.30,020.50,020.32,019.72,
     +018.50,016.54,013.84,010.68,007.14,004.40,004.74,007.58,
     +010.74,013.48,014.00,015.54,016.74,017.46,017.94,018.30,
     +018.50,018.58,018.32,017.72,016.50,014.64,012.24,009.18,
     +005.84,002.90,003.30,006.16,009.14,011.84,012.30,013.78,
     +014.94,015.66,016.24,016.50,016.70,016.70,016.42,015.78,
     +014.60,012.90,010.66,007.86,004.88,001.60,001.72,004.96,
     +007.84,010.24,010.70,012.14,013.24,013.96,014.44,014.80,
     +014.90,014.88,014.52,013.92,012.80,011.30,009.28,006.94,
     +004.32,001.80,001.94,004.34,006.78,008.94,009.40,010.58,
     +011.64,012.36,012.74,013.10,013.20,013.08,012.72,012.12,
     +011.10,009.86,008.30,006.50,004.60,003.10,003.16,004.50,
     +006.20,007.90,008.40,009.42,010.14,010.76,011.14,011.40,
     +011.40,011.38,011.02,010.46,009.70,008.72,007.64,006.46,
     +005.42,004.60,004.70,005.34,006.24,007.36,007.90,008.46,
     +008.92,009.28,009.54,009.70,009.70,009.68,009.42,009.06,
     +008.60,008.08,007.56,007.02,006.56,006.30,006.30,006.52,
     +006.96,007.38,008.15,008.15,008.15,008.15,008.15,008.15,
     +008.15,008.15,008.15,008.15,008.15,008.15,008.15,008.15,
     +008.15,008.15,008.15,008.15,008.15,008.15/

C     Data Input      
      rlan = rga
      rlo = rgo      
      
C     From "normal" geographic latitude 
C     to angle from South Pole.  
      rla = rlan + 90

      IF (rlo .EQ. 360) THEN
         rlo = 0
        END IF

C     PROXIMITY

C     coefficients of the latitudinal points
      LA1 = (INT(rla/2)+1)
      LA2 = LA1 + 1
      if(la2.gt.91) la2=91

C     coefficients of the longitudinal points
      LO1 = (INT(rlo/18)+1)
corr      LO2 = LO1 + 1
      LO2 = MOD(LO1,20) + 1 

C     Four points of Geomagnetic Coordinates
      gm1 = CORMAG(LO1,LA1)
      gm2 = CORMAG(LO1,LA2) 
      gm3 = CORMAG(LO2,LA1)
      gm4 = CORMAG(LO2,LA2)

C     latitudinal points
      X1 = ABS(rla - (INT(rla)))                        
      X2 = 2. - X1

C     longitudinal points
      Y1 = ABS(rlo - (INT(rlo)))
      Y2 = 18. - Y1
      
C     X AND Y VALUES
      x = X1 / (X1 + X2)
      y = Y1 / (Y1 + Y2)

C     INTERPOLATION
      gmla = gm1 * (1 - x) * (1 - y) + gm2 * (1 - y) * (x) + gm3 * (y)
     1 * (1 - x) + gm4 * (x) * (y)

C     OUTPUT OF THE PROGRAM
C     From corrected geomagnetic latitude from North Pole
C     to "normal"  geomagnetic latitude.       
      rgma = 90. - gmla

      END
c
C--------------------------------------------------
      real function peakh(foE,foF2,M3000)
      real MF,M3000
      sqM=M3000*M3000
      MF=M3000*sqrt((0.0196*sqM+1.)/(1.2967*sqM-1.0))
      If(foE.ge.1.0E-30) then
         ratio=foF2/foE
         ratio=djoin(ratio,1.75,20.0,ratio-1.75)
         dM=0.253/(ratio-1.215)-0.012
      else
         dM=-0.012
      endif
      peakh=1490.0*MF/(M3000+dM)-176.0
      return
      end
C-----------------------------------------------------------------
         real function djoin(f1,f2,alpha,x)
      real f1,f2,alpha,x,ee,fexp
      ee=fexp(alpha*x)
      djoin=(f1*ee+f2)/(ee+1.0)
      return
      end
C-----------------------------------------------------------------
       real function fexp(a)
      real a
      if(a.gt.80.0) then
         fexp=5.5406E34
         return
      endif
      if(a.lt.-80.0D0) then
         fexp=1.8049E-35
         return
      endif
      fexp=exp(a)
      return
      end
C-----------------------------------------------------------------

