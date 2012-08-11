
      FUNCTION TBDY(A0,AI,CN,CM,T,D)

      B0 = A0*(300./T)**CN
      BI = AI*(300./T)**CM
      Y = ALOG10(B0*D/BI)
      X = 1./(1. + Y**2)
      TBDY = B0*D/(1. + B0*D/BI) * 0.6**X
      RETURN
      END
