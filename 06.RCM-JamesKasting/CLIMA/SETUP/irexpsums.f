
      SUBROUTINE IREXPSUMS

c-mm  This subroutine will READ in the exponential sum datafiles in the
c-mm  ir for both h2o and co2 and output the appropriate value for each
c-mm  wavelength interval.  This will crash if our values exceed the
c-mm  boundaries set forth in the data, but this shouldn't happen under
c-mm  normal Martian conditions.
c-mm  Indices for xkappa and kappa are as follows
c-mm  Index #1:  Temp   1=less than T   2=more than T
c-mm  Index #2:  Pres   1=less than p   2=more than p
c-mm  Index #3:  lam    wavelength bin
c-mm  Index #4:  i      gauss point
c-mm  Index #5:  Spec   1=co2  2=h2o  3=ch4

      INTEGER LINE,i,lam,junk,Tind,pind,Tinddum,j,k,l,m
      REAL zeroed
      INCLUDE '../INCLUDECLIM/comIRDATA.inc'
     
      do ii=1,3
       read(36,*)
      enddo
      do 9400 i=1,8
         do 9401 j=1,12
            do 9402 k=1,55
               do 9403 l=1,8
                  do 9404 m=1,3
                     xkappa(i,j,k,l,m)=0.0
 9404             continue
 9403          continue
 9402       continue
 9401    continue
 9400 continue
      READ(36,9950)
        do 9951 i=1,8
          READ(36,9952) weight(i,1)
 9951   continue
        do 9700 i=1,6
           READ(36,*)
 9700   continue

        do 9500 Tind=1,7
           do 9502 pind=1,11
              do 9957 lam=9,48

c-mm  If zeroed=0.0, that implies there is only one set of numbers
c-mm  for that wavelength bin.  Otherwise, there are two and we
c-mm  READ from the second set.  xkappa is the particular kappa
c-mm  (1-8) for the given indices.  This will be later interpolated

          READ(36,9956) zeroed
          READ(36,*)
          do 9955 i=1,8
             if (zeroed.ne.(0.0)) then
                READ(36,9958) xkappa(Tind,pind,lam,i,1)
             else
                READ(36,9960) xkappa(Tind,pind,lam,i,1)
             endif
 9955     continue
        READ(36,*)
 9957   continue
        READ(36,*)
 9502   continue
 9500   continue

c-mm   H2O below this line----------------------
        READ(36,9950)
        do 9901 i=1,8
           READ(36,9952) weight(i,2)
 9901       continue
        do 9701 i=1,6
           READ(36,*)
 9701      continue
        do 9504 Tind=1,7
           do 9506 pind=1,11
              do 9508 lam=1,55
 
c-mm  If zeroed=0.0, that implies there is only one set of numbers
c-mm  for that wavelength bin.  Otherwise, there are two and we
c-mm  READ from the second set.  xkappa is the particular kappa
c-mm  (1-8) for the given indices.  This will be later interpolated
 
          READ(36,9956) zeroed
          READ(36,*)
          do 9510 i=1,8
             if (zeroed.ne.(0.0)) then
                READ(36,9958) xkappa(Tind,pind,lam,i,2)
             else
                READ(36,9960) xkappa(Tind,pind,lam,i,2)
             endif
 9510     continue
        READ(36,*)
 9508 continue
        READ(36,*)
 9506 continue
 9504 continue

c-mm  CH4 below this line--------------------------------
      weight(1,3)=0.08566225
      weight(2,3)=0.18038079
      weight(3,3)=0.23395697
      weight(4,3)=0.23395697
      weight(5,3)=0.18038079
      weight(6,3)=0.08566225
      do ii=1,3 
       READ(36,*)
      enddo
      do 9600 Tinddum=1,3
         Tind=Tinddum*3-Tinddum/3
         do 9601 pind=7,12
            do 9602 lam=1,38
               READ(36,*) junk,(xkappa(Tind,pind,lam,ii,3), ii=1,6)
 9602          continue
              READ(36,*)
              READ(36,*)
              READ(36,*)
 9601          continue
 9600          continue

        close(36)

        do 9570 i=7,12
        do 9571 j=1,38
        do 9572 k=1,6
          xkappa(1,i,j,k,3)=xkappa(3,i,j,k,3)
 9572 continue
 9571 continue
 9570 continue
 
 9950   format(/)
 9952   format(27X,E16.10)
 9956   format(83X,E9.3)
 9958   format(56X,E12.6)
 9960   format(21X,E12.6)
 
        RETURN
        END
