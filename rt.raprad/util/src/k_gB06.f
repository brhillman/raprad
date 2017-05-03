      PARAMETER (MG = 16)
      REAL KA(5,13,MG)
      DIMENSION SELFREF(10,MG)

c      COMMON /HVRSNB/ HVRKG(NBANDS)
      COMMON /K6/ KA ,SELFREF

c      CHARACTER*8 HVRKG

c      DATA HVRKG(6)  / '2.4' /

C     The array KA contains absorption coefs at the 16 chosen g-values 
C     for a range of pressure levels > ~100mb and temperatures.  The first
C     index in the array, JT, which runs from 1 to 5, corresponds to 
C     different temperatures.  More specifically, JT = 3 means that the 
C     data are for the corresponding TREF for this  pressure level, 
C     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
C     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
C     index, JP, runs from 1 to 13 and refers to the corresponding 
C     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
C     The third index, IG, goes from 1 to 16, and tells us which 
C     g-interval the absorption coefficients are for.
      DATA (KA(JT, 1, 1),JT=1,5) /
     &2.2568E-05,2.5521E-05,2.4392E-05,2.3686E-05,2.3850E-05/
      DATA (KA(JT, 2, 1),JT=1,5) /
     &1.8872E-05,1.8763E-05,1.8106E-05,1.7702E-05,1.7681E-05/
      DATA (KA(JT, 3, 1),JT=1,5) /
     &1.2996E-05,1.3086E-05,1.3378E-05,1.3260E-05,1.1771E-05/
      DATA (KA(JT, 4, 1),JT=1,5) /
     &9.6276E-06,9.9774E-06,9.6694E-06,8.4861E-06,8.9160E-06/
      DATA (KA(JT, 5, 1),JT=1,5) /
     &7.1690E-06,7.2984E-06,5.8013E-06,6.2225E-06,6.6628E-06/
      DATA (KA(JT, 6, 1),JT=1,5) /
     &5.6038E-06,4.4095E-06,4.4722E-06,4.4395E-06,5.4630E-06/
      DATA (KA(JT, 7, 1),JT=1,5) /
     &3.4775E-06,3.0437E-06,3.3876E-06,3.8636E-06,4.9257E-06/
      DATA (KA(JT, 8, 1),JT=1,5) /
     &2.7310E-06,3.2224E-06,3.6552E-06,4.3975E-06,5.5359E-06/
      DATA (KA(JT, 9, 1),JT=1,5) /
     &6.6325E-06,7.7076E-06,8.9070E-06,1.0179E-05,1.1909E-05/
      DATA (KA(JT,10, 1),JT=1,5) /
     &2.4079E-05,2.6989E-05,3.0753E-05,3.4229E-05,3.8060E-05/
      DATA (KA(JT,11, 1),JT=1,5) /
     &3.3970E-05,3.7965E-05,4.2993E-05,4.7630E-05,5.2629E-05/
      DATA (KA(JT,12, 1),JT=1,5) /
     &3.5907E-05,4.0095E-05,4.5267E-05,5.0045E-05,5.5178E-05/
      DATA (KA(JT,13, 1),JT=1,5) /
     &3.0600E-05,3.4162E-05,3.8530E-05,4.2579E-05,4.6937E-05/
      DATA (KA(JT, 1, 2),JT=1,5) /
     &2.9639E-05,3.3221E-05,3.2806E-05,2.8964E-05,2.6726E-05/
      DATA (KA(JT, 2, 2),JT=1,5) /
     &2.3653E-05,2.4355E-05,2.4028E-05,2.0887E-05,2.2925E-05/
      DATA (KA(JT, 3, 2),JT=1,5) /
     &1.6734E-05,1.6216E-05,1.4570E-05,1.6223E-05,1.7762E-05/
      DATA (KA(JT, 4, 2),JT=1,5) /
     &1.0615E-05,1.0115E-05,1.1508E-05,1.1885E-05,1.3929E-05/
      DATA (KA(JT, 5, 2),JT=1,5) /
     &7.5406E-06,7.3085E-06,8.6006E-06,9.3888E-06,1.1719E-05/
      DATA (KA(JT, 6, 2),JT=1,5) /
     &4.8289E-06,6.2547E-06,5.8818E-06,7.9263E-06,1.0542E-05/
      DATA (KA(JT, 7, 2),JT=1,5) /
     &4.3622E-06,4.7605E-06,5.4656E-06,7.2224E-06,1.0212E-05/
      DATA (KA(JT, 8, 2),JT=1,5) /
     &4.2923E-06,4.6584E-06,6.2268E-06,8.3462E-06,1.1797E-05/
      DATA (KA(JT, 9, 2),JT=1,5) /
     &9.7259E-06,1.1243E-05,1.3832E-05,1.8290E-05,2.4869E-05/
      DATA (KA(JT,10, 2),JT=1,5) /
     &3.3596E-05,3.8056E-05,4.3895E-05,5.4699E-05,7.1889E-05/
      DATA (KA(JT,11, 2),JT=1,5) /
     &4.7132E-05,5.3142E-05,6.1177E-05,7.5188E-05,9.8046E-05/
      DATA (KA(JT,12, 2),JT=1,5) /
     &4.9755E-05,5.6067E-05,6.4461E-05,7.8769E-05,1.0238E-04/
      DATA (KA(JT,13, 2),JT=1,5) /
     &4.2448E-05,4.7895E-05,5.4990E-05,6.7099E-05,8.7152E-05/
      DATA (KA(JT, 1, 3),JT=1,5) /
     &2.4391E-05,2.6179E-05,3.0228E-05,3.9877E-05,4.1318E-05/
      DATA (KA(JT, 2, 3),JT=1,5) /
     &2.0506E-05,2.0015E-05,2.6885E-05,2.9284E-05,2.9837E-05/
      DATA (KA(JT, 3, 3),JT=1,5) /
     &1.3542E-05,1.5721E-05,1.9501E-05,1.8947E-05,2.4397E-05/
      DATA (KA(JT, 4, 3),JT=1,5) /
     &9.0261E-06,1.3633E-05,1.2485E-05,1.7278E-05,2.0007E-05/
      DATA (KA(JT, 5, 3),JT=1,5) /
     &8.7945E-06,9.7631E-06,1.1634E-05,1.3600E-05,1.9745E-05/
      DATA (KA(JT, 6, 3),JT=1,5) /
     &7.5275E-06,7.1532E-06,9.4521E-06,1.3168E-05,1.9311E-05/
      DATA (KA(JT, 7, 3),JT=1,5) /
     &5.6021E-06,6.6005E-06,8.8148E-06,1.3299E-05,1.9798E-05/
      DATA (KA(JT, 8, 3),JT=1,5) /
     &5.6076E-06,7.6899E-06,1.0237E-05,1.5421E-05,2.3499E-05/
      DATA (KA(JT, 9, 3),JT=1,5) /
     &1.3625E-05,1.7122E-05,2.2283E-05,3.1506E-05,4.8193E-05/
      DATA (KA(JT,10, 3),JT=1,5) /
     &4.6234E-05,5.4855E-05,6.8590E-05,8.9964E-05,1.3093E-04/
      DATA (KA(JT,11, 3),JT=1,5) /
     &6.5003E-05,7.6686E-05,9.5408E-05,1.2507E-04,1.8079E-04/
      DATA (KA(JT,12, 3),JT=1,5) /
     &6.8759E-05,8.1014E-05,1.0079E-04,1.3196E-04,1.9017E-04/
      DATA (KA(JT,13, 3),JT=1,5) /
     &5.8741E-05,6.9297E-05,8.6397E-05,1.1292E-04,1.6288E-04/
      DATA (KA(JT, 1, 4),JT=1,5) /
     &3.6952E-05,3.9167E-05,4.5764E-05,4.7701E-05,5.5871E-05/
      DATA (KA(JT, 2, 4),JT=1,5) /
     &2.6008E-05,3.6209E-05,2.8658E-05,3.8449E-05,3.8668E-05/
      DATA (KA(JT, 3, 4),JT=1,5) /
     &2.1596E-05,2.1370E-05,2.2422E-05,2.7123E-05,3.0425E-05/
      DATA (KA(JT, 4, 4),JT=1,5) /
     &1.4408E-05,1.3094E-05,1.9460E-05,2.1634E-05,3.3208E-05/
      DATA (KA(JT, 5, 4),JT=1,5) /
     &8.9120E-06,1.2874E-05,1.5169E-05,2.3591E-05,3.5032E-05/
      DATA (KA(JT, 6, 4),JT=1,5) /
     &7.0427E-06,1.1114E-05,1.4987E-05,2.2856E-05,3.4461E-05/
      DATA (KA(JT, 7, 4),JT=1,5) /
     &7.4919E-06,9.5389E-06,1.5038E-05,2.3474E-05,3.5967E-05/
      DATA (KA(JT, 8, 4),JT=1,5) /
     &8.6186E-06,1.1413E-05,1.7846E-05,2.8336E-05,4.4637E-05/
      DATA (KA(JT, 9, 4),JT=1,5) /
     &1.9661E-05,2.5521E-05,3.7417E-05,5.9323E-05,9.6439E-05/
      DATA (KA(JT,10, 4),JT=1,5) /
     &6.3994E-05,8.1146E-05,1.0850E-04,1.6659E-04,2.7179E-04/
      DATA (KA(JT,11, 4),JT=1,5) /
     &9.0362E-05,1.1392E-04,1.5220E-04,2.3565E-04,3.8761E-04/
      DATA (KA(JT,12, 4),JT=1,5) /
     &9.5893E-05,1.2069E-04,1.6199E-04,2.5340E-04,4.1839E-04/
      DATA (KA(JT,13, 4),JT=1,5) /
     &8.2144E-05,1.0344E-04,1.3998E-04,2.2132E-04,3.6483E-04/
      DATA (KA(JT, 1, 5),JT=1,5) /
     &3.3532E-05,4.1940E-05,5.7696E-05,6.0186E-05,6.3029E-05/
      DATA (KA(JT, 2, 5),JT=1,5) /
     &2.3830E-05,3.5587E-05,4.6056E-05,4.1505E-05,5.7960E-05/
      DATA (KA(JT, 3, 5),JT=1,5) /
     &2.1812E-05,2.5706E-05,2.8233E-05,3.5229E-05,5.8555E-05/
      DATA (KA(JT, 4, 5),JT=1,5) /
     &1.7315E-05,1.9695E-05,2.2530E-05,3.8597E-05,5.8973E-05/
      DATA (KA(JT, 5, 5),JT=1,5) /
     &1.1890E-05,1.3867E-05,2.3935E-05,3.9940E-05,5.8993E-05/
      DATA (KA(JT, 6, 5),JT=1,5) /
     &1.0677E-05,1.4192E-05,2.4839E-05,4.0149E-05,5.9735E-05/
      DATA (KA(JT, 7, 5),JT=1,5) /
     &9.2784E-06,1.5528E-05,2.6050E-05,4.1864E-05,6.4186E-05/
      DATA (KA(JT, 8, 5),JT=1,5) /
     &1.1359E-05,1.8672E-05,3.1354E-05,5.1149E-05,8.0918E-05/
      DATA (KA(JT, 9, 5),JT=1,5) /
     &2.6841E-05,4.0944E-05,6.8221E-05,1.1515E-04,1.8821E-04/
      DATA (KA(JT,10, 5),JT=1,5) /
     &8.7415E-05,1.2053E-04,1.9441E-04,3.3166E-04,5.5740E-04/
      DATA (KA(JT,11, 5),JT=1,5) /
     &1.2386E-04,1.7143E-04,2.7887E-04,4.7716E-04,8.0477E-04/
      DATA (KA(JT,12, 5),JT=1,5) /
     &1.3204E-04,1.8364E-04,3.0111E-04,5.1685E-04,8.7520E-04/
      DATA (KA(JT,13, 5),JT=1,5) /
     &1.1358E-04,1.5915E-04,2.6323E-04,4.5343E-04,7.7073E-04/
      DATA (KA(JT, 1, 6),JT=1,5) /
     &5.6626E-05,6.2952E-05,5.6997E-05,7.8387E-05,1.4023E-04/
      DATA (KA(JT, 2, 6),JT=1,5) /
     &4.2639E-05,4.0472E-05,5.1521E-05,7.2367E-05,1.2444E-04/
      DATA (KA(JT, 3, 6),JT=1,5) /
     &3.3191E-05,2.9821E-05,4.1035E-05,7.5835E-05,1.1292E-04/
      DATA (KA(JT, 4, 6),JT=1,5) /
     &1.7563E-05,2.2621E-05,4.3898E-05,7.3541E-05,1.0514E-04/
      DATA (KA(JT, 5, 6),JT=1,5) /
     &1.4086E-05,2.3545E-05,4.5696E-05,7.1090E-05,1.0366E-04/
      DATA (KA(JT, 6, 6),JT=1,5) /
     &1.2590E-05,2.4048E-05,4.4731E-05,6.9646E-05,1.0367E-04/
      DATA (KA(JT, 7, 6),JT=1,5) /
     &1.3429E-05,2.6924E-05,4.6412E-05,7.4051E-05,1.1343E-04/
      DATA (KA(JT, 8, 6),JT=1,5) /
     &1.8185E-05,3.3168E-05,5.7618E-05,9.5901E-05,1.5240E-04/
      DATA (KA(JT, 9, 6),JT=1,5) /
     &4.0915E-05,7.1523E-05,1.2899E-04,2.2341E-04,3.6849E-04/
      DATA (KA(JT,10, 6),JT=1,5) /
     &1.2433E-04,2.0393E-04,3.7307E-04,6.6855E-04,1.1425E-03/
      DATA (KA(JT,11, 6),JT=1,5) /
     &1.8043E-04,3.0126E-04,5.5504E-04,9.9966E-04,1.7159E-03/
      DATA (KA(JT,12, 6),JT=1,5) /
     &1.9679E-04,3.3338E-04,6.1765E-04,1.1153E-03,1.9175E-03/
      DATA (KA(JT,13, 6),JT=1,5) /
     &1.7363E-04,2.9831E-04,5.5435E-04,1.0025E-03,1.7267E-03/
      DATA (KA(JT, 1, 7),JT=1,5) /
     &5.4728E-05,6.3216E-05,9.7227E-05,1.9830E-04,3.0115E-04/
      DATA (KA(JT, 2, 7),JT=1,5) /
     &3.7424E-05,4.6891E-05,9.3922E-05,1.8575E-04,2.6393E-04/
      DATA (KA(JT, 3, 7),JT=1,5) /
     &3.2885E-05,5.4240E-05,1.1210E-04,1.7073E-04,2.3633E-04/
      DATA (KA(JT, 4, 7),JT=1,5) /
     &2.6560E-05,6.2341E-05,1.0718E-04,1.5691E-04,2.1887E-04/
      DATA (KA(JT, 5, 7),JT=1,5) /
     &2.8005E-05,6.1158E-05,9.7995E-05,1.4461E-04,2.0623E-04/
      DATA (KA(JT, 6, 7),JT=1,5) /
     &2.6304E-05,5.4858E-05,8.8237E-05,1.3524E-04,2.0068E-04/
      DATA (KA(JT, 7, 7),JT=1,5) /
     &2.8057E-05,5.2849E-05,8.8492E-05,1.4199E-04,2.1851E-04/
      DATA (KA(JT, 8, 7),JT=1,5) /
     &3.4390E-05,6.3412E-05,1.1163E-04,1.8592E-04,2.9376E-04/
      DATA (KA(JT, 9, 7),JT=1,5) /
     &7.2790E-05,1.4039E-04,2.5926E-04,4.5257E-04,7.4795E-04/
      DATA (KA(JT,10, 7),JT=1,5) /
     &2.1021E-04,4.1798E-04,8.0790E-04,1.4658E-03,2.5048E-03/
      DATA (KA(JT,11, 7),JT=1,5) /
     &3.2720E-04,6.6181E-04,1.2928E-03,2.3683E-03,4.0754E-03/
      DATA (KA(JT,12, 7),JT=1,5) /
     &3.7693E-04,7.7171E-04,1.5176E-03,2.7982E-03,4.8311E-03/
      DATA (KA(JT,13, 7),JT=1,5) /
     &3.5194E-04,7.2710E-04,1.4369E-03,2.6518E-03,4.5779E-03/
      DATA (KA(JT, 1, 8),JT=1,5) /
     &8.8622E-05,2.2200E-04,4.1760E-04,5.9916E-04,8.1306E-04/
      DATA (KA(JT, 2, 8),JT=1,5) /
     &7.7484E-05,1.9083E-04,3.2988E-04,4.7824E-04,6.5145E-04/
      DATA (KA(JT, 3, 8),JT=1,5) /
     &6.8249E-05,1.8888E-04,2.9496E-04,4.2064E-04,5.7727E-04/
      DATA (KA(JT, 4, 8),JT=1,5) /
     &1.0307E-04,1.8157E-04,2.7450E-04,3.9411E-04,5.4734E-04/
      DATA (KA(JT, 5, 8),JT=1,5) /
     &1.0013E-04,1.6578E-04,2.5412E-04,3.7165E-04,5.2629E-04/
      DATA (KA(JT, 6, 8),JT=1,5) /
     &8.7661E-05,1.4688E-04,2.3121E-04,3.4816E-04,5.0439E-04/
      DATA (KA(JT, 7, 8),JT=1,5) /
     &7.9171E-05,1.3672E-04,2.2306E-04,3.4738E-04,5.2076E-04/
      DATA (KA(JT, 8, 8),JT=1,5) /
     &8.3889E-05,1.5157E-04,2.5895E-04,4.2268E-04,6.6361E-04/
      DATA (KA(JT, 9, 8),JT=1,5) /
     &1.6488E-04,3.1976E-04,5.8317E-04,1.0082E-03,1.6627E-03/
      DATA (KA(JT,10, 8),JT=1,5) /
     &4.7725E-04,9.7517E-04,1.8746E-03,3.3981E-03,5.7914E-03/
      DATA (KA(JT,11, 8),JT=1,5) /
     &7.9688E-04,1.6560E-03,3.2130E-03,5.8470E-03,9.9511E-03/
      DATA (KA(JT,12, 8),JT=1,5) /
     &9.7776E-04,2.0671E-03,4.0593E-03,7.3885E-03,1.2564E-02/
      DATA (KA(JT,13, 8),JT=1,5) /
     &9.8124E-04,2.1024E-03,4.1582E-03,7.5870E-03,1.2921E-02/
      DATA (KA(JT, 1, 9),JT=1,5) /
     &7.4912E-04,1.1099E-03,1.5644E-03,2.1464E-03,2.8897E-03/
      DATA (KA(JT, 2, 9),JT=1,5) /
     &6.2597E-04,9.2912E-04,1.3142E-03,1.8181E-03,2.4773E-03/
      DATA (KA(JT, 3, 9),JT=1,5) /
     &4.7892E-04,7.1171E-04,1.0243E-03,1.4397E-03,1.9820E-03/
      DATA (KA(JT, 4, 9),JT=1,5) /
     &3.8453E-04,5.9140E-04,8.7342E-04,1.2581E-03,1.7584E-03/
      DATA (KA(JT, 5, 9),JT=1,5) /
     &3.6098E-04,5.7696E-04,8.8197E-04,1.2954E-03,1.8431E-03/
      DATA (KA(JT, 6, 9),JT=1,5) /
     &3.4542E-04,5.7060E-04,8.9328E-04,1.3402E-03,1.9357E-03/
      DATA (KA(JT, 7, 9),JT=1,5) /
     &3.3888E-04,5.8373E-04,9.4422E-04,1.4513E-03,2.1408E-03/
      DATA (KA(JT, 8, 9),JT=1,5) /
     &3.7707E-04,6.8408E-04,1.1538E-03,1.8357E-03,2.7831E-03/
      DATA (KA(JT, 9, 9),JT=1,5) /
     &7.0739E-04,1.3735E-03,2.4588E-03,4.1112E-03,6.4788E-03/
      DATA (KA(JT,10, 9),JT=1,5) /
     &1.8575E-03,3.8603E-03,7.2941E-03,1.2720E-02,2.0839E-02/
      DATA (KA(JT,11, 9),JT=1,5) /
     &2.9958E-03,6.2197E-03,1.1775E-02,2.0597E-02,3.3906E-02/
      DATA (KA(JT,12, 9),JT=1,5) /
     &3.6668E-03,7.6016E-03,1.4363E-02,2.5216E-02,4.1719E-02/
      DATA (KA(JT,13, 9),JT=1,5) /
     &3.7360E-03,7.7341E-03,1.4600E-02,2.5684E-02,4.2617E-02/
      DATA (KA(JT, 1,10),JT=1,5) /
     &2.2204E-03,3.1215E-03,4.1592E-03,5.5251E-03,7.2462E-03/
      DATA (KA(JT, 2,10),JT=1,5) /
     &1.8794E-03,2.6515E-03,3.6524E-03,4.9115E-03,6.4747E-03/
      DATA (KA(JT, 3,10),JT=1,5) /
     &1.4689E-03,2.0942E-03,2.9433E-03,4.0441E-03,5.4540E-03/
      DATA (KA(JT, 4,10),JT=1,5) /
     &1.1490E-03,1.6495E-03,2.3294E-03,3.2354E-03,4.4282E-03/
      DATA (KA(JT, 5,10),JT=1,5) /
     &8.6570E-04,1.3179E-03,1.9604E-03,2.8379E-03,4.0441E-03/
      DATA (KA(JT, 6,10),JT=1,5) /
     &8.5652E-04,1.4272E-03,2.2416E-03,3.3643E-03,4.8824E-03/
      DATA (KA(JT, 7,10),JT=1,5) /
     &9.3347E-04,1.6326E-03,2.6435E-03,4.0816E-03,6.0203E-03/
      DATA (KA(JT, 8,10),JT=1,5) /
     &1.1389E-03,2.0806E-03,3.5187E-03,5.6514E-03,8.6161E-03/
      DATA (KA(JT, 9,10),JT=1,5) /
     &2.3318E-03,4.5800E-03,8.2437E-03,1.3780E-02,2.1744E-02/
      DATA (KA(JT,10,10),JT=1,5) /
     &6.2326E-03,1.3165E-02,2.5232E-02,4.4615E-02,7.3492E-02/
      DATA (KA(JT,11,10),JT=1,5) /
     &9.8704E-03,2.0942E-02,4.0187E-02,7.0789E-02,1.1609E-01/
      DATA (KA(JT,12,10),JT=1,5) /
     &1.1999E-02,2.5297E-02,4.8337E-02,8.4705E-02,1.3826E-01/
      DATA (KA(JT,13,10),JT=1,5) /
     &1.1789E-02,2.4642E-02,4.6747E-02,8.1645E-02,1.3284E-01/
      DATA (KA(JT, 1,11),JT=1,5) /
     &3.1261E-03,4.2813E-03,5.7811E-03,7.6623E-03,1.0021E-02/
      DATA (KA(JT, 2,11),JT=1,5) /
     &2.7730E-03,3.9256E-03,5.2937E-03,7.0385E-03,9.1561E-03/
      DATA (KA(JT, 3,11),JT=1,5) /
     &2.2211E-03,3.2332E-03,4.4680E-03,6.0516E-03,7.9945E-03/
      DATA (KA(JT, 4,11),JT=1,5) /
     &1.7260E-03,2.5271E-03,3.5974E-03,5.0208E-03,6.7957E-03/
      DATA (KA(JT, 5,11),JT=1,5) /
     &1.3640E-03,2.0174E-03,2.8988E-03,4.1037E-03,5.6676E-03/
      DATA (KA(JT, 6,11),JT=1,5) /
     &1.1158E-03,1.8312E-03,2.8771E-03,4.3425E-03,6.2982E-03/
      DATA (KA(JT, 7,11),JT=1,5) /
     &1.2321E-03,2.1507E-03,3.5466E-03,5.5748E-03,8.3259E-03/
      DATA (KA(JT, 8,11),JT=1,5) /
     &1.5691E-03,2.9089E-03,5.0376E-03,8.1424E-03,1.2407E-02/
      DATA (KA(JT, 9,11),JT=1,5) /
     &3.5241E-03,7.0356E-03,1.2741E-02,2.1461E-02,3.3839E-02/
      DATA (KA(JT,10,11),JT=1,5) /
     &1.0170E-02,2.1914E-02,4.2506E-02,7.5273E-02,1.2386E-01/
      DATA (KA(JT,11,11),JT=1,5) /
     &1.6549E-02,3.5671E-02,6.9003E-02,1.2202E-01,2.0089E-01/
      DATA (KA(JT,12,11),JT=1,5) /
     &2.0073E-02,4.2914E-02,8.2517E-02,1.4566E-01,2.3923E-01/
      DATA (KA(JT,13,11),JT=1,5) /
     &1.9754E-02,4.1795E-02,8.0135E-02,1.4074E-01,2.2980E-01/
      DATA (KA(JT, 1,12),JT=1,5) /
     &4.1564E-03,5.9925E-03,8.3120E-03,1.1209E-02,1.4597E-02/
      DATA (KA(JT, 2,12),JT=1,5) /
     &3.9567E-03,5.4909E-03,7.6176E-03,1.0371E-02,1.3856E-02/
      DATA (KA(JT, 3,12),JT=1,5) /
     &3.3974E-03,4.7965E-03,6.6252E-03,8.9930E-03,1.2050E-02/
      DATA (KA(JT, 4,12),JT=1,5) /
     &2.7271E-03,4.0009E-03,5.6537E-03,7.7402E-03,1.0333E-02/
      DATA (KA(JT, 5,12),JT=1,5) /
     &2.1241E-03,3.2132E-03,4.6676E-03,6.6066E-03,9.0309E-03/
      DATA (KA(JT, 6,12),JT=1,5) /
     &1.6495E-03,2.5482E-03,3.8834E-03,5.7278E-03,8.1979E-03/
      DATA (KA(JT, 7,12),JT=1,5) /
     &1.5776E-03,2.7768E-03,4.6314E-03,7.2948E-03,1.0930E-02/
      DATA (KA(JT, 8,12),JT=1,5) /
     &2.1394E-03,4.0088E-03,6.9630E-03,1.1391E-02,1.7592E-02/
      DATA (KA(JT, 9,12),JT=1,5) /
     &5.1508E-03,1.0419E-02,1.9252E-02,3.2952E-02,5.2981E-02/
      DATA (KA(JT,10,12),JT=1,5) /
     &1.6525E-02,3.6418E-02,7.1731E-02,1.2850E-01,2.1301E-01/
      DATA (KA(JT,11,12),JT=1,5) /
     &2.8843E-02,6.3252E-02,1.2384E-01,2.2076E-01,3.6471E-01/
      DATA (KA(JT,12,12),JT=1,5) /
     &3.6629E-02,7.9908E-02,1.5563E-01,2.7617E-01,4.5458E-01/
      DATA (KA(JT,13,12),JT=1,5) /
     &3.6750E-02,7.9425E-02,1.5346E-01,2.7137E-01,4.4636E-01/
      DATA (KA(JT, 1,13),JT=1,5) /
     &6.2842E-03,8.5901E-03,1.1770E-02,1.5848E-02,2.1231E-02/
      DATA (KA(JT, 2,13),JT=1,5) /
     &5.8832E-03,8.3961E-03,1.1688E-02,1.5821E-02,2.1081E-02/
      DATA (KA(JT, 3,13),JT=1,5) /
     &5.0693E-03,7.3487E-03,1.0502E-02,1.4654E-02,1.9671E-02/
      DATA (KA(JT, 4,13),JT=1,5) /
     &4.3064E-03,6.2164E-03,8.9638E-03,1.2692E-02,1.7434E-02/
      DATA (KA(JT, 5,13),JT=1,5) /
     &3.5860E-03,5.3245E-03,7.6791E-03,1.0860E-02,1.5079E-02/
      DATA (KA(JT, 6,13),JT=1,5) /
     &2.8314E-03,4.3212E-03,6.4657E-03,9.2062E-03,1.2862E-02/
      DATA (KA(JT, 7,13),JT=1,5) /
     &2.2377E-03,3.6343E-03,5.9331E-03,9.1864E-03,1.3701E-02/
      DATA (KA(JT, 8,13),JT=1,5) /
     &2.8263E-03,5.3372E-03,9.3113E-03,1.5184E-02,2.3750E-02/
      DATA (KA(JT, 9,13),JT=1,5) /
     &7.1632E-03,1.4853E-02,2.8075E-02,4.9370E-02,8.0768E-02/
      DATA (KA(JT,10,13),JT=1,5) /
     &2.5324E-02,5.6755E-02,1.1356E-01,2.0699E-01,3.4818E-01/
      DATA (KA(JT,11,13),JT=1,5) /
     &4.7551E-02,1.0638E-01,2.1201E-01,3.8439E-01,6.4430E-01/
      DATA (KA(JT,12,13),JT=1,5) /
     &6.5536E-02,1.4596E-01,2.8971E-01,5.2186E-01,8.7097E-01/
      DATA (KA(JT,13,13),JT=1,5) /
     &7.1040E-02,1.5739E-01,3.0996E-01,5.5635E-01,9.2532E-01/
      DATA (KA(JT, 1,14),JT=1,5) /
     &8.9073E-03,1.2577E-02,1.7104E-02,2.3298E-02,3.0760E-02/
      DATA (KA(JT, 2,14),JT=1,5) /
     &9.1449E-03,1.2712E-02,1.7315E-02,2.3477E-02,3.1268E-02/
      DATA (KA(JT, 3,14),JT=1,5) /
     &8.6458E-03,1.1907E-02,1.6470E-02,2.2100E-02,2.9699E-02/
      DATA (KA(JT, 4,14),JT=1,5) /
     &7.6642E-03,1.0832E-02,1.5035E-02,2.0201E-02,2.7549E-02/
      DATA (KA(JT, 5,14),JT=1,5) /
     &6.6532E-03,9.4973E-03,1.3436E-02,1.8486E-02,2.5410E-02/
      DATA (KA(JT, 6,14),JT=1,5) /
     &5.6659E-03,8.2188E-03,1.1655E-02,1.6380E-02,2.2833E-02/
      DATA (KA(JT, 7,14),JT=1,5) /
     &4.7388E-03,7.0553E-03,1.0065E-02,1.4371E-02,2.0476E-02/
      DATA (KA(JT, 8,14),JT=1,5) /
     &4.0074E-03,6.9317E-03,1.1650E-02,1.9270E-02,3.0109E-02/
      DATA (KA(JT, 9,14),JT=1,5) /
     &9.4290E-03,1.9762E-02,3.7876E-02,6.6863E-02,1.1024E-01/
      DATA (KA(JT,10,14),JT=1,5) /
     &3.5639E-02,8.0302E-02,1.6119E-01,2.9567E-01,5.0285E-01/
      DATA (KA(JT,11,14),JT=1,5) /
     &7.0994E-02,1.5991E-01,3.2149E-01,5.8889E-01,1.0009E+00/
      DATA (KA(JT,12,14),JT=1,5) /
     &1.0472E-01,2.3582E-01,4.7386E-01,8.6654E-01,1.4714E+00/
      DATA (KA(JT,13,14),JT=1,5) /
     &1.2287E-01,2.7626E-01,5.5413E-01,1.0143E+00,1.7166E+00/
      DATA (KA(JT, 1,15),JT=1,5) /
     &1.0915E-02,1.5406E-02,2.1787E-02,3.0490E-02,4.2541E-02/
      DATA (KA(JT, 2,15),JT=1,5) /
     &1.1315E-02,1.6166E-02,2.2798E-02,3.1923E-02,4.3992E-02/
      DATA (KA(JT, 3,15),JT=1,5) /
     &1.0941E-02,1.5701E-02,2.1967E-02,3.0980E-02,4.2137E-02/
      DATA (KA(JT, 4,15),JT=1,5) /
     &1.0488E-02,1.4945E-02,2.1284E-02,2.9684E-02,4.1028E-02/
      DATA (KA(JT, 5,15),JT=1,5) /
     &9.8738E-03,1.4407E-02,2.0449E-02,2.8357E-02,3.9999E-02/
      DATA (KA(JT, 6,15),JT=1,5) /
     &8.9244E-03,1.3460E-02,1.9534E-02,2.7000E-02,3.8027E-02/
      DATA (KA(JT, 7,15),JT=1,5) /
     &8.0056E-03,1.2241E-02,1.8433E-02,2.6084E-02,3.5910E-02/
      DATA (KA(JT, 8,15),JT=1,5) /
     &7.2557E-03,1.0997E-02,1.6847E-02,2.4941E-02,3.5979E-02/
      DATA (KA(JT, 9,15),JT=1,5) /
     &1.1206E-02,2.3492E-02,4.5074E-02,7.9749E-02,1.3179E-01/
      DATA (KA(JT,10,15),JT=1,5) /
     &4.3714E-02,9.8800E-02,1.9893E-01,3.6488E-01,6.2367E-01/
      DATA (KA(JT,11,15),JT=1,5) /
     &9.0564E-02,2.0463E-01,4.1227E-01,7.5790E-01,1.2930E+00/
      DATA (KA(JT,12,15),JT=1,5) /
     &1.4027E-01,3.1640E-01,6.3737E-01,1.1714E+00,1.9962E+00/
      DATA (KA(JT,13,15),JT=1,5) /
     &1.7360E-01,3.9198E-01,7.8920E-01,1.4509E+00,2.4685E+00/
      DATA (KA(JT, 1,16),JT=1,5) /
     &1.1824E-02,1.7312E-02,2.6044E-02,4.2271E-02,6.3915E-02/
      DATA (KA(JT, 2,16),JT=1,5) /
     &1.2333E-02,1.8642E-02,2.7212E-02,4.2155E-02,7.0740E-02/
      DATA (KA(JT, 3,16),JT=1,5) /
     &1.1914E-02,1.8039E-02,2.6483E-02,3.8087E-02,6.0410E-02/
      DATA (KA(JT, 4,16),JT=1,5) /
     &1.1617E-02,1.6736E-02,2.5551E-02,3.7243E-02,5.5742E-02/
      DATA (KA(JT, 5,16),JT=1,5) /
     &1.2207E-02,1.6058E-02,2.3899E-02,3.5249E-02,5.2521E-02/
      DATA (KA(JT, 6,16),JT=1,5) /
     &1.2552E-02,1.6221E-02,2.2453E-02,3.3451E-02,4.8959E-02/
      DATA (KA(JT, 7,16),JT=1,5) /
     &1.2682E-02,1.6583E-02,2.1720E-02,3.1647E-02,4.6528E-02/
      DATA (KA(JT, 8,16),JT=1,5) /
     &1.2800E-02,1.6916E-02,2.2044E-02,3.0239E-02,4.5143E-02/
      DATA (KA(JT, 9,16),JT=1,5) /
     &1.2606E-02,2.5371E-02,4.8585E-02,8.6042E-02,1.4214E-01/
      DATA (KA(JT,10,16),JT=1,5) /
     &4.7940E-02,1.0811E-01,2.1696E-01,3.9856E-01,6.7922E-01/
      DATA (KA(JT,11,16),JT=1,5) /
     &1.0038E-01,2.2692E-01,4.5634E-01,8.3917E-01,1.4277E+00/
      DATA (KA(JT,12,16),JT=1,5) /
     &1.5780E-01,3.5540E-01,7.1875E-01,1.3169E+00,2.2414E+00/
      DATA (KA(JT,13,16),JT=1,5) /
     &1.9907E-01,4.4847E-01,9.0510E-01,1.6605E+00,2.8238E+00/

C     The array SELFREF contains the coefficient of the water vapor
C     self-continuum (including the energy term).  The first index
C     refers to temperature in 7.2 degree increments.  For instance,
C     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
C     etc.  The second index runs over the g-channel (1 to 16).
      DATA (SELFREF(JT, 1),JT=1,10)  /
     & 7.71391E-02, 6.43175E-02, 5.36271E-02, 4.47136E-02, 3.72816E-02,
     & 3.10849E-02, 2.59182E-02, 2.16102E-02, 1.80183E-02, 1.50234E-02/
      DATA (SELFREF(JT, 2),JT=1,10)  /
     & 8.39275E-02, 7.04019E-02, 5.90561E-02, 4.95388E-02, 4.15552E-02,
     & 3.48583E-02, 2.92406E-02, 2.45283E-02, 2.05754E-02, 1.72595E-02/
      DATA (SELFREF(JT, 3),JT=1,10)  /
     & 8.77747E-02, 7.41411E-02, 6.26251E-02, 5.28978E-02, 4.46814E-02,
     & 3.77413E-02, 3.18791E-02, 2.69275E-02, 2.27449E-02, 1.92121E-02/
      DATA (SELFREF(JT, 4),JT=1,10)  /
     & 9.35855E-02, 7.94116E-02, 6.73843E-02, 5.71787E-02, 4.85187E-02,
     & 4.11704E-02, 3.49349E-02, 2.96439E-02, 2.51542E-02, 2.13445E-02/
      DATA (SELFREF(JT, 5),JT=1,10)  /
     & 9.31058E-02, 7.97629E-02, 6.83321E-02, 5.85395E-02, 5.01503E-02,
     & 4.29633E-02, 3.68062E-02, 3.15316E-02, 2.70128E-02, 2.31416E-02/
      DATA (SELFREF(JT, 6),JT=1,10)  /
     & 7.92454E-02, 7.03070E-02, 6.23769E-02, 5.53412E-02, 4.90991E-02,
     & 4.35611E-02, 3.86477E-02, 3.42885E-02, 3.04210E-02, 2.69898E-02/
      DATA (SELFREF(JT, 7),JT=1,10)  /
     & 7.29322E-02, 6.57740E-02, 5.93183E-02, 5.34963E-02, 4.82457E-02,
     & 4.35104E-02, 3.92399E-02, 3.53886E-02, 3.19152E-02, 2.87828E-02/
      DATA (SELFREF(JT, 8),JT=1,10)  /
     & 7.53680E-02, 6.63817E-02, 5.84669E-02, 5.14958E-02, 4.53558E-02,
     & 3.99480E-02, 3.51849E-02, 3.09897E-02, 2.72948E-02, 2.40404E-02/
      DATA (SELFREF(JT, 9),JT=1,10)  /
     & 7.17087E-02, 6.34545E-02, 5.61505E-02, 4.96871E-02, 4.39678E-02,
     & 3.89068E-02, 3.44283E-02, 3.04654E-02, 2.69586E-02, 2.38554E-02/
      DATA (SELFREF(JT,10),JT=1,10)  /
     & 7.20376E-02, 6.42071E-02, 5.72277E-02, 5.10070E-02, 4.54625E-02,
     & 4.05207E-02, 3.61161E-02, 3.21903E-02, 2.86912E-02, 2.55724E-02/
      DATA (SELFREF(JT,11),JT=1,10)  /
     & 8.56200E-02, 7.44248E-02, 6.46934E-02, 5.62344E-02, 4.88815E-02,
     & 4.24900E-02, 3.69342E-02, 3.21049E-02, 2.79070E-02, 2.42581E-02/
      DATA (SELFREF(JT,12),JT=1,10)  /
     & 9.30364E-02, 8.00633E-02, 6.88992E-02, 5.92918E-02, 5.10240E-02,
     & 4.39092E-02, 3.77864E-02, 3.25174E-02, 2.79832E-02, 2.40811E-02/
      DATA (SELFREF(JT,13),JT=1,10)  /
     & 9.58838E-02, 8.17150E-02, 6.96399E-02, 5.93492E-02, 5.05792E-02,
     & 4.31051E-02, 3.67354E-02, 3.13070E-02, 2.66808E-02, 2.27381E-02/
      DATA (SELFREF(JT,14),JT=1,10)  /
     & 8.92791E-02, 7.69435E-02, 6.63123E-02, 5.71500E-02, 4.92536E-02,
     & 4.24483E-02, 3.65833E-02, 3.15286E-02, 2.71723E-02, 2.34179E-02/
      DATA (SELFREF(JT,15),JT=1,10)  /
     & 9.09947E-02, 7.81455E-02, 6.71107E-02, 5.76340E-02, 4.94956E-02,
     & 4.25064E-02, 3.65041E-02, 3.13494E-02, 2.69226E-02, 2.31209E-02/
      DATA (SELFREF(JT,16),JT=1,10)  /
     & 7.99990E-02, 7.02289E-02, 6.16519E-02, 5.41224E-02, 4.75125E-02,
     & 4.17099E-02, 3.66159E-02, 3.21440E-02, 2.82183E-02, 2.47720E-02/

