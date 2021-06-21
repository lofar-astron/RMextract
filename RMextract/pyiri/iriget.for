      SUBROUTINE IRI_GET_SUB(path,JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
c iriget.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
c
c test program for the iri_sub subroutine
c     
      LOGICAL JF(50)
      DIMENSION OARR(100),OUTF(20,1000)
      CHARACTER path*200
      call set_datapath(path)
      call read_ig_rz
      call readapf107
      
c calling IRI subroutine
c 
      call iri_sub(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
      
      RETURN
      end


      SUBROUTINE IRI_GET_TEC(path,JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,tec,tecb,tect)
c iriget.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
c
c test program for the iri_tec subroutine
c     
      LOGICAL JF(50)
      CHARACTER path*200
      call set_datapath(path)
      call read_ig_rz
      call readapf107
      
c calling IRI subroutine
c 

      call irit13(ALATI,ALONG,jmag,jf,iyyyy,mmdd,dhour,heibeg,heiend,
     &                          tec,tecb,tect)
      RETURN
      end
