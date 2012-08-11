 
      SUBROUTINE READPHOTO

      INCLUDE '../INCLUDECHEM/parNZ.inc'
      DIMENSION XL(10)
      INCLUDE '../INCLUDECHEM/comQBLOK.inc'
      INCLUDE '../INCLUDECHEM/comFLUXPHOTO.inc'

C ***** READ THE PHOTOLYSIS DATAFILE *****
c-as  This file contains the Sun flux from Atmospheric Ozone 1985
c-as   Vol. 1 World Meteorogical Organization. Global Ozone 
c-as   Research and Monitoring Project. Report No. 16. 
c
      READ(62,100)
 100  FORMAT(/)
      DO 41 L=1,108
  41  READ(62,101) WAVL(L),WAVU(L),FLUX(L),SO31(L),SO32(L),SMCHO(L),
     &  SMCOM(L)

 101  FORMAT(8X,F6.1,1X,F6.1,2X,E9.2,4(3X,E8.1))

      READ(62,102)
 102  FORMAT(////)
      DO 42 L=1,17
  42  READ(62,103) KA(L),(A(L,K),K=1,9)
 103  FORMAT(4X,I1,1X,9(1X,E13.6))
C
      READ(62,102)
      DO 43 L=1,17
  43  READ(62,104) KB(L),(B(L,K),K=1,5)
 104  FORMAT(4X,I1,1X,5(1X,E13.6))
C
      READ(62,102)
      DO 44 L=1,35
  44  READ(62,105) SO2HZ(L),SH2O(L),SCO2(L),SHO2(L),SN2O(L),SHCL(L),
     2  SCCL2F2(L),SCHCLF2(L),SCH3CL(L),SMCCL3(L)
 105  FORMAT(5X,10(E8.1,1X))
C
      READ(62,106)
 106  FORMAT(//)
      DO 45 L=1,68
  45  READ(62,107) SH2O2(L),SHNO3(L),SHNO4(L),SPAN(L),SPPN(L),
     2  SMOOH(L),SNO2(L),SH2CO(L),RHCO(L),RH2(L)
 107  FORMAT(6X,8E8.1,2(1X,F4.2))
C
      READ(62,106)
      DO 46 L=1,68
  46  READ(62,108) SCL2(L),SHOCL(L),SCLNO(L),SCLONO(L),SCLONO2(L),
     2  SCCL4(L),SCLO2(L),SCLO3(L),SCCL3F(L)
 108  FORMAT(6X,9E8.1)
C
      READ(62,109)
 109  FORMAT(//////)
      DO 58 L=1,17
      READ(62,110) NK(L),(ALPHA(L,K),K=1,4)
 110  FORMAT(6X,I1,1X,4F13.5)
  58  READ(62,111) (BETA(L,K),K=1,4)
 111  FORMAT(8X,4E13.4/)
C
      READ(62,210)
C   (Sulfur species cross sections)
 210  FORMAT(//)
C-AP
      DO 67 L=1,68
  67  READ(62,113) SSO2(L),SSO21(L),SSO23(L),SSO(L),SH2S(L)
 113  FORMAT(15X,5(E9.3,4X))
C-AP
      READ(62,212)
 212  FORMAT(//)
      DO 151 L=1,68
 151  READ(62,213) SN2O5(L),SCL2O2(L)
 213  FORMAT(6X,E9.2,2X,E9.2)
      close(62)      

********************
c   READ absorption coefficients for the far UV
      READ(73,100)
      do L=1,10
      read(73,*)X,(SIGMA(J,L),J=1,3),SSO2A(L),SSOA(L),SCH4(L)
      end do
      close(73) 

      RETURN
      END




