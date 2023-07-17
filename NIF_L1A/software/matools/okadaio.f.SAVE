      subroutine  okadaio(alp, xrec, source, uout)
C
C*********************************************************
C*****                                               ***** 
C*****    subrouttine to pass data to SRECTF         ***** 
C*****                    BY  Y.OKADA                *****
C*****                                               *****
C********************************************************* 
      real*8  ALP,X,Y,DEP,AL1,AL2,AW1,AW2,SD,CD,DISL1,DISL2,DISL3
      real*8  U1,U2,U3,U11,U12,U21,U22,U31,U32
      real*8  xrec(2), source(10), uout(9)
C
C                                                         
C***** INPUT                                                          
C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)
C*****   X,Y     : COORDINATE OF STATION                            
C*****   DEP     : SOURCE DEPTH                                    
C*****   AL1,AL2 : FAULT LENGTH RANGE                             
C*****   AW1,AW2 : FAULT WIDTH RANGE                             
C*****   SD,CD   : SIN,COS OF DIP-ANGLE                         
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION     
C                                                                     
C***** OUTPUT                                                        
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL     )  
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /     
C*****   U31,U32         : TILT                 UNIT OF X,Y,,,AW )   
C                                                                   
      x = xrec(1)
      y = xrec(2)
      DEP  = source(1)
      AL1  = source(2)
      AL2  = source(3)
      AW1  = source(4)
      AW2  = source(5)
      SD   = source(6)
      CD   = source(7)
      DISL1= source(8)
      DISL2= source(9)
      DISL3= source(10)
C
      call  SRECTF(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,          
     *                   SD,CD,DISL1,DISL2,DISL3, 
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32) 
C
      uout(1) = U1 
      uout(2) = U2 
      uout(3) = U3 
      uout(4) = U11  
      uout(5) = U12  
      uout(6) = U21  
      uout(7) = U22  
      uout(8) = U31  
      uout(9) = U32  
C
      RETURN 
      END   
      SUBROUTINE  SRECTF(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,                   01310003
     *                   SD,CD,DISL1,DISL2,DISL3,                       01311003
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              01320000
      IMPLICIT REAL*8 (A-H,O-Z)                                         01330000
C                                                                       01340000
C*********************************************************              01350000
C*****                                               *****              01360000
C*****    SURFACE DISPLACEMENT,STRAIN,TILT           *****              01370000
C*****    DUE TO RECTANGULAR FAULT IN A HALF-SPACE   *****              01380000
C*****              CODED BY  Y.OKADA ... JAN 1985   *****              01390000
C*****                                               *****              01400000
C*********************************************************              01410000
C                                                                       01420000
C***** INPUT                                                            01430000
C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)   01431003
C*****   X,Y     : COORDINATE OF STATION                                01450003
C*****   DEP     : SOURCE DEPTH                                         01460003
C*****   AL1,AL2 : FAULT LENGTH RANGE                                   01470003
C*****   AW1,AW2 : FAULT WIDTH RANGE                                    01471003
C*****   SD,CD   : SIN,COS OF DIP-ANGLE                                 01480003
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)01490000
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION      01500000
C                                                                       01510000
C***** OUTPUT                                                           01520000
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL     )      01530000
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /          01540000
C*****   U31,U32         : TILT                 UNIT OF X,Y,,,AW )      01550000
C                                                                       01560000
C***** SUBROUTINE USED...SRECTG                                         01570000
C                                                                       01580000
      DIMENSION  U(9),DU(9)                                             01590000
      DATA  F0, F1 / 0.D0, 1.D0 /                                       01600000
C-----                                                                  01610000
      P = Y*CD + DEP*SD                                                 01620000
      Q = Y*SD - DEP*CD                                                 01630000
      DO 1111  I=1,9                                                    01640000
 1111 U(I)=F0                                                           01650000
C-----                                                                  01660000
      DO 5555  K=1,2                                                    01670000
       IF(K.EQ.1)  ET=P-AW1                                             01680003
       IF(K.EQ.2)  ET=P-AW2                                             01690003
       DO 4444  J=1,2                                                   01700000
        IF(J.EQ.1)  XI=X-AL1                                            01710003
        IF(J.EQ.2)  XI=X-AL2                                            01720003
        JK=J+K                                                          01730000
        IF(JK.NE.3)  SIGN= F1                                           01740000
        IF(JK.EQ.3)  SIGN=-F1                                           01750000
        CALL SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3,                01760000
     *           DU(1),DU(2),DU(3),DU(4),DU(5),DU(6),DU(7),DU(8),DU(9)) 01770000
        DO 3333  I=1,9                                                  01780000
         U(I)=U(I)+SIGN*DU(I)                                           01790000
 3333   CONTINUE                                                        01800000
 4444  CONTINUE                                                         01810000
 5555 CONTINUE                                                          01820000
      U1 =U(1)                                                          01830000
      U2 =U(2)                                                          01840000
      U3 =U(3)                                                          01850000
      U11=U(4)                                                          01860000
      U12=U(5)                                                          01870000
      U21=U(6)                                                          01880000
      U22=U(7)                                                          01890000
      U31=U(8)                                                          01900000
      U32=U(9)                                                          01910000
      RETURN                                                            01920000
      END                                                               01930000
      SUBROUTINE  SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3,           01940000
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              01950000
      IMPLICIT REAL*8 (A-H,O-Z)                                         01960000
C                                                                       01970000
C*********************************************************************  01980000
C*****                                                           *****  01990000
C*****  INDEFINITE INTEGRAL OF SURFACE DISPLACEMENT,STRAIN,TILT  *****  02000000
C*****  DUE TO RECTANGULAR FAULT IN A HALF-SPACE                 *****  02010000
C*****                          CODED BY  Y.OKADA ... JAN 1985   *****  02020000
C*****                                                           *****  02030000
C*********************************************************************  02040000
C                                                                       02050000
C***** INPUT                                                            02060000
C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)   02061003
C*****   XI,ET,Q : FAULT COORDINATE                                     02080000
C*****   SD,CD   : SIN,COS OF DIP-ANGLE                                 02090000
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)02100000
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION      02110000
C                                                                       02120000
C***** OUTPUT                                                           02130000
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL    )       02140000
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /          02150000
C*****   U31,U32         : TILT                 UNIT OF XI,ET,Q )       02160000
C                                                                       02170000
      DATA  F0,F1,F2/ 0.D0, 1.D0, 2.D0 /                                02180000
      PI2=6.283185307179586D0                                           02190000
C-----                                                                  02200000
      XI2=XI*XI                                                         02210000
      ET2=ET*ET                                                         02220000
      Q2=Q*Q                                                            02230000
      R2=XI2+ET2+Q2                                                     02240000
      R =DSQRT(R2)                                                      02250000
      R3=R*R2                                                           02260000
      D =ET*SD-Q*CD                                                     02270000
      Y =ET*CD+Q*SD                                                     02280000
      RET=R+ET                                                          02290000
      IF(RET.LT.F0)  RET=F0                                             02300000
      RD =R+D                                                           02310000
      RRD=F1/(R*RD)                                                     02320000
C-----                                                                  02330000
      IF( Q .NE.F0)  TT = DATAN( XI*ET/(Q*R) )                          02340000
      IF( Q .EQ.F0)  TT = F0                                            02350000
      IF(RET.NE.F0)  RE = F1/RET                                        02360000
      IF(RET.EQ.F0)  RE = F0                                            02370000
      IF(RET.NE.F0)  DLE= DLOG(RET)                                     02380000
      IF(RET.EQ.F0)  DLE=-DLOG(R-ET)                                    02390000
      RRX=F1/(R*(R+XI))                                                 02400000
      RRE=RE/R                                                          02410000
      AXI=(F2*R+XI)*RRX*RRX/R                                           02420000
      AET=(F2*R+ET)*RRE*RRE/R                                           02430000
      IF(CD.EQ.F0)  GO TO 20                                            02440000
C==============================                                         02450000
C=====   INCLINED FAULT   =====                                         02460000
C==============================                                         02470000
      TD=SD/CD                                                          02480000
      X =DSQRT(XI2+Q2)                                                  02490000
      IF(XI.EQ.F0)  A5=F0                                               02500000
      IF(XI.NE.F0)                                                      02510000
     *A5= ALP*F2/CD*DATAN( (ET*(X+Q*CD)+X*(R+X)*SD) / (XI*(R+X)*CD) )   02520000
      A4= ALP/CD*( DLOG(RD) - SD*DLE )                                  02530000
      A3= ALP*(Y/RD/CD - DLE) + TD*A4                                   02540000
      A1=-ALP/CD*XI/RD        - TD*A5                                   02550000
      C1= ALP/CD*XI*(RRD - SD*RRE)                                      02560000
      C3= ALP/CD*(Q*RRE - Y*RRD)                                        02570000
      B1= ALP/CD*(XI2*RRD - F1)/RD - TD*C3                              02580000
      B2= ALP/CD*XI*Y*RRD/RD       - TD*C1                              02590000
      GO TO 30                                                          02600000
C==============================                                         02610000
C=====   VERTICAL FAULT   =====                                         02620000
C==============================                                         02630000
   20 RD2=RD*RD                                                         02640000
      A1=-ALP/F2*XI*Q/RD2                                               02650000
      A3= ALP/F2*( ET/RD + Y*Q/RD2 - DLE )                              02660000
      A4=-ALP*Q/RD                                                      02670000
      A5=-ALP*XI*SD/RD                                                  02680000
      B1= ALP/F2*  Q  /RD2*(F2*XI2*RRD - F1)                            02690000
      B2= ALP/F2*XI*SD/RD2*(F2*Q2 *RRD - F1)                            02700000
      C1= ALP*XI*Q*RRD/RD                                               02710000
      C3= ALP*SD/RD*(XI2*RRD - F1)                                      02720000
C-----                                                                  02730000
   30 A2=-ALP*DLE - A3                                                  02740000
      B3=-ALP*XI*RRE - B2                                               02750000
      B4=-ALP*( CD/R + Q*SD*RRE ) - B1                                  02760000
      C2= ALP*(-SD/R + Q*CD*RRE ) - C3                                  02770000
C-----                                                                  02780000
      U1 =F0                                                            02790000
      U2 =F0                                                            02800000
      U3 =F0                                                            02810000
      U11=F0                                                            02820000
      U12=F0                                                            02830000
      U21=F0                                                            02840000
      U22=F0                                                            02850000
      U31=F0                                                            02860000
      U32=F0                                                            02870000
C======================================                                 02880000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 02890000
C======================================                                 02900000
      IF(DISL1.EQ.F0)  GO TO 200                                        02910000
      UN=DISL1/PI2                                                      02920000
      REQ=RRE*Q                                                         02930000
      U1 =U1 - UN*( REQ*XI +   TT    + A1*SD )                          02940000
      U2 =U2 - UN*( REQ*Y  + Q*CD*RE + A2*SD )                          02950000
      U3 =U3 - UN*( REQ*D  + Q*SD*RE + A4*SD )                          02960000
      U11=U11+ UN*( XI2*Q*AET - B1*SD )                                 02970000
      U12=U12+ UN*( XI2*XI*( D/(ET2+Q2)/R3 - AET*SD ) - B2*SD )         02980000
      U21=U21+ UN*( XI*Q/R3*CD + (XI*Q2*AET - B2)*SD )                  02990000
      U22=U22+ UN*( Y *Q/R3*CD + (Q*SD*(Q2*AET-F2*RRE)                  03000000
     *                            -(XI2+ET2)/R3*CD - B4)*SD )           03010000
      U31=U31+ UN*(-XI*Q2*AET*CD + (XI*Q/R3 - C1)*SD )                  03020000
      U32=U32+ UN*( D*Q/R3*CD + (XI2*Q*AET*CD - SD/R + Y*Q/R3 - C2)*SD )03030000
C===================================                                    03040000
C=====  DIP-SLIP CONTRIBUTION  =====                                    03050000
C===================================                                    03060000
  200 IF(DISL2.EQ.F0)  GO TO 300                                        03070000
      UN=DISL2/PI2                                                      03080000
      SDCD=SD*CD                                                        03090000
      U1 =U1 - UN*( Q/R             - A3*SDCD )                         03100000
      U2 =U2 - UN*( Y*Q*RRX + CD*TT - A1*SDCD )                         03110000
      U3 =U3 - UN*( D*Q*RRX + SD*TT - A5*SDCD )                         03120000
      U11=U11+ UN*( XI*Q/R3            + B3*SDCD )                      03130000
      U12=U12+ UN*( Y *Q/R3 - SD/R     + B1*SDCD )                      03140000
      U21=U21+ UN*( Y *Q/R3 + Q*CD*RRE + B1*SDCD )                      03150000
      U22=U22+ UN*( Y*Y*Q*AXI - (F2*Y*RRX + XI*CD*RRE)*SD + B2*SDCD )   03160000
      U31=U31+ UN*( D *Q/R3 + Q*SD*RRE + C3*SDCD )                      03170000
      U32=U32+ UN*( Y*D*Q*AXI - (F2*D*RRX + XI*SD*RRE)*SD + C1*SDCD )   03180000
C========================================                               03190000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               03200000
C========================================                               03210000
  300 IF(DISL3.EQ.F0)  GO TO 900                                        03220000
      UN=DISL3/PI2                                                      03230000
      SDSD=SD*SD                                                        03240000
      U1 =U1 + UN*( Q2*RRE                       - A3*SDSD )            03250000
      U2 =U2 + UN*(-D*Q*RRX - SD*(XI*Q*RRE - TT) - A1*SDSD )            03260000
      U3 =U3 + UN*( Y*Q*RRX + CD*(XI*Q*RRE - TT) - A5*SDSD )            03270000
      U11=U11- UN*( XI*Q2*AET             + B3*SDSD )                   03280000
      U12=U12- UN*(-D*Q/R3 - XI2*Q*AET*SD + B1*SDSD )                   03290000
      U21=U21- UN*( Q2*(CD/R3 + Q*AET*SD) + B1*SDSD )                   03300000
      U22=U22- UN*((Y*CD-D*SD)*Q2*AXI - F2*Q*SD*CD*RRX                  03310000
     *                      - (XI*Q2*AET - B2)*SDSD )                   03320000
      U31=U31- UN*( Q2*(SD/R3 - Q*AET*CD) + C3*SDSD )                   03330000
      U32=U32- UN*((Y*SD+D*CD)*Q2*AXI + XI*Q2*AET*SD*CD                 03340000
     *                       - (F2*Q*RRX - C1)*SDSD )                   03350000
C-----                                                                  03360000
  900 RETURN                                                            03370000
      END                                                               03380000
