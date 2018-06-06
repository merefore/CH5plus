subroutine walkercalcpot(local,nwalk,v,force)
	implicit none
	integer parmax,i
	double precision par(2)
	double precision, dimension(nwalk,10), intent(in) :: local
	integer, intent(in) :: nwalk
	double precision, dimension(nwalk), intent(out) :: v
	double precision, dimension(289), intent(out) :: force
	parmax = 110
	par(1) =           1.08594310
	par(2) =           1.84500000
!	force = 0.0d0
	force(3) =         3.19449869
	force(4) =       -10.34265977
	force(5) =     13365.56120756
	force(6) =     14522.52943848
	force(7) =      2865.66161462
	force(8) =       337.58839095
	force(9) =     40097.23679417
	force(10) =     13631.85433924
	force(11) =     -7644.30455936
	force(12) =       720.34838066
	force(13) =     -1197.83512846
	force(14) =     -1261.78859591
	force(15) =     -1575.53747450
	force(16) =    -1678.29016090
	force(17) =       2141.05957306
	force(18) =     -2684.39248481
	force(19) =        162.32893323
	force(20) =     -1979.38524595
	force(21) =        431.79424659
	force(22) =     -1616.07856310
	force(23) =       7517.55490718
	force(24) =        759.00378306
	force(25) =       4063.87998951
	force(26) =       3363.86975012
	force(27) =        334.88721174
	force(28) =        243.24530228
	force(29) =        969.18377042
	force(30) =       3381.93726426
	force(31) =     -2756.99088461
	force(32) =     -1160.72120437
	force(33) =       -766.85014511
	force(34) =        630.57394254
	force(35) =     -1138.55154168
	force(37) =       -676.41035587
	force(38) =       -145.04375785
	force(39) =       1170.90020653
	force(41) =       -101.26342738
	force(42) =        764.83623499
	force(43) =       -764.86622656
	force(44) =     -2408.85733254
	force(47) =       1203.45476220
	force(49) =     -1184.61431830
	force(51) =       -173.03874968
	force(52) =       -925.84454682
	force(53) =     -1090.90365667
	force(54) =     -1801.47555429
	force(55) =       2113.30077847
	force(56) =     30169.56814085
	force(57) =        326.35158276
	force(58) =     -1845.36983813
	force(59) =       -412.79540327
	force(60) =     19054.54416729
	force(61) =       2044.42657596
	force(63) =        153.02915468
	force(65) =       -557.74874364
	force(66) =       -875.22329191
	force(70) =         69.79639086
	force(76) =        629.04183401
	force(77) =        619.24816391
	force(78) =       1189.09371942
	force(85) =        936.50662942
	force(101) =     -1399.38536071
	force(102) =        956.47802196
	force(103) =        650.91787668
	force(104) =       -914.84850835
	force(105) =       -609.73751072
	force(106) =       -725.76801668
	force(107) =     -1011.87706341
	force(109) =     -1208.20197858
	force(113) =       -845.31633878
	force(114) =     -1770.99230049
	force(115) =     -1717.64535620
	force(117) =        584.01287825
	force(118) =       1696.96131573
	force(119) =     -3936.61587120
	force(120) =     -2584.37486858
	force(121) =       -722.50510494
	force(122) =       -171.38251217
	force(123) =        932.55923441
	force(128) =     10046.17195932
	force(129) =          7.22205526
	force(130) =     68366.01263817
	force(131) =       8997.36217189
	force(133) =     -3604.76155373
	force(134) =       5066.23917498
	force(135) =        834.35057355
	force(136) =    -15006.01079852
	force(137) =       1070.72787585
	force(138) =       1635.97703755
	force(139) =       1519.97675751
	force(140) =       1362.57981414
	force(141) =          2.52893341
	force(143) =     -2307.27237284
	force(144) =     -4057.09492824
	force(147) =        438.83273254
	force(151) =       -326.03316845
	force(154) =       1769.01587843
	force(165) =        599.81801957
	force(167) =     -1249.14473646
	force(168) =       1904.78023677
	force(170) =        996.31931610
	force(183) =     -1211.98148917
	force(184) =        -13.23461123
	force(243) =     -1433.66237498
	force(246) =       -168.35830718
	force(249) =       -496.33403249
	force(284) =       -912.33207568
	force(285) =     -1516.27178469
	force(288) =        -73.98447090
	force(289) =       1320.51318762
! 	local(1) = 1.0860100000  
!	local(2) = 1.0860100000  
!	local(3) = 1.0860100000  
!	local(4) = 1.0860100000  
!	local(5) = 109.4712206300  
!	local(6) = 109.4712206300  
!	local(7) = 109.4712206300  
!	local(8) = 109.4712206300  
!	local(9) = 109.4712206300  
!	local(10) = 109.4712206300

	do i = 1,nwalk
		call poten_xy4(local(i,:), parmax, par, force,v(i))
	enddo
	
end subroutine walkercalcpot


subroutine poten_xy4(local,parmax,par,force,f)
  implicit none

  integer, parameter :: ark         = selected_real_kind(25,32)
  integer, parameter :: ik          = selected_int_kind(8)

  integer parmax, ipar, term(289), k
  double precision local(10), force(289), par(2), dF(289)
  double precision f

  integer(ik)          :: N

  real(ark) :: y1, y2, y3, y4, y5, y6, y7, y8, y9, r1e, alphae, a0, deg, pi
  real(ark) :: s1,s2,s3
  real(ark) :: r1,r2,r3,r4,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34

  pi = 4.0d0 * datan2(1.0d0,1.0d0)

  deg = pi/180.0d0

  r1e    = par(1)
  a0     = par(2)

  r1  = local(1)
  r2  = local(2)
  r3  = local(3)
  r4  = local(4)

  alpha12 = local(5)*deg
  alpha13 = local(6)*deg
  alpha14 = local(7)*deg
  alpha23 = local(8)*deg
  alpha24 = local(9)*deg
  alpha34 = local(10)*deg
!  print*,r1,r2,r3,r4
!	print*,alpha12,alpha13
  y1=1.0_ark-exp(-a0*(r1-r1e))
  y2=1.0_ark-exp(-a0*(r2-r1e))
  y3=1.0_ark-exp(-a0*(r3-r1e))
  y4=1.0_ark-exp(-a0*(r4-r1e))
      
  y5=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
  y6=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
  y7=(alpha24-alpha13)/sqrt(2.0_ark)
  y8=(alpha23-alpha14)/sqrt(2.0_ark)
  y9=(alpha34-alpha12)/sqrt(2.0_ark)

!  print*,'y',y1,y2,y3,y4,y5,y6,y7,y8,y9
      dF(1) = 0._ark
      dF(2) = 0._ark
      dF(3) = 1.0_ark
      dF(4) = y2+y3+y4+y1
      dF(5) = y8**2+y7**2+y9**2
      dF(6) = y6**2+y5**2
      dF(7) = (-y7-y8-y9)*y1+(y7-y9+y8)*y2+(y8+y9-y7)*y3+(y9+y7-y8)*y4
      dF(8) = (y4+y3+y2)*y1+(y4+y3)*y2+y3*y4
      dF(9) = y2**2+y3**2+y4**2+y1**2
      dF(10) = y7*y8*y9
      dF(11) = (-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6
      dF(12) = y5**3-3._ark*y5*y6**2
      dF(13) = ((y8-2._ark*y9+y7)*y5+(sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y6)*y1+((-y8-2._ark*y9- &
          y7)*y5+(sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y6)*y2+((2._ark*y9+y7-y8)*y5+(-sqrt(3._ark)*y8- &
          sqrt(3._ark)*y7)*y6)*y3+((2._ark*y9+y8-y7)*y5+(sqrt(3._ark)*y7+sqrt(3._ark)*y8)*y6)*y4
      dF(14) = ((y9+y8)*y7+y8*y9)*y1+((y8-y9)*y7-y8*y9)*y2+((-y9-y8)*y7+y8*y9)*y3+((- &
          y8+y9)*y7-y8*y9)*y4
      dF(15) = (y8**2+y7**2+y9**2)*y1+(y8**2+y7**2+y9**2)*y2+(y8**2+y7**2+y9**2)*y3+ &
          (y8**2+y7**2+y9**2)*y4
      dF(16) = (y6**2+y5**2)*y1+(y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4
      dF(17) = (y3*y7+y4*y8+y2*y9)*y1+(-y3*y8-y4*y7)*y2-y3*y4*y9
      dF(18) = (y2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4)*y1+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+sqrt(3._ark)*y6/ &
          2._ark)*y4)*y2+y3*y4*y5
      dF(19) = ((y4+y3)*y2+y3*y4)*y1+y2*y3*y4
      dF(20) = (y9+y7+y8)*y1**2+(-y7+y9-y8)*y2**2+(-y8-y9+y7)*y3**2+(-y7-y9+y8)*y4**2
      dF(21) = (y4+y3+y2)*y1**2+(y3**2+y2**2+y4**2)*y1+(y4+y3)*y2**2+(y3**2+y4**2)*y2+ &
          y3**2*y4+y3*y4**2
      dF(22) = y3**3+y1**3+y4**3+y2**3
      dF(23) = (y9**2+y8**2)*y7**2+y8**2*y9**2
      dF(24) = y9**4+y8**4+y7**4
      dF(25) = -sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2
      dF(26) = (y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2
      dF(27) = y6**4+y5**4+2._ark*y5**2*y6**2
      dF(28) = (y3**2+y2**2+y4**2)*y1**2+(y3**2+y4**2)*y2**2+y3**2*y4**2
      dF(29) = (-y7-y8-y9)*y1**3+(y7-y9+y8)*y2**3+(y8+y9-y7)*y3**3+(y9+y7-y8)*y4**3
      dF(30) = y4**4+y3**4+y2**4+y1**4
      dF(31) = y3*y7*y8*y9+y2*y7*y8*y9+y4*y7*y8*y9+y1*y7*y8*y9
      dF(32) = ((y9+y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y1+((-y8+y9)*y7**2+ &
          (-y9**2-y8**2)*y7+y8**2*y9-y8*y9**2)*y2+((-y9-y8)*y7**2+(y9**2+y8**2)*y7- &
          y8**2*y9-y8*y9**2)*y3+((y8-y9)*y7**2+(-y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y4
      dF(33) = (y9**3+y8**3+y7**3)*y1+(y9**3-y8**3-y7**3)*y2+(-y8**3-y9**3+y7**3)*y3+ &
          (-y7**3-y9**3+y8**3)*y4
      dF(34) = ((-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y1+((-sqrt(3._ark)*y7**2/3._ark- &
          sqrt(3._ark)*y8**2/3._ark+2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y2+((- &
          sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+ &
          (y7**2-y8**2)*y6)*y3+((-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y4
      dF(35) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4
      dF(36) = (-sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(sqrt(3._ark)*y9/6._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y1+(-sqrt(3._ark)*y5**2*y9/2._ark+(y8- &
          y7)*y6*y5+(sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y6**2)*y2+ &
          (sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark- &
          sqrt(3._ark)*y9/6._ark)*y6**2)*y3+(sqrt(3._ark)*y5**2*y9/2._ark+(-y8-y7)*y6*y5+(- &
          sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y4
      dF(37) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3+((-y7-y9+y8)*y5**2+(-y7- &
          y9+y8)*y6**2)*y4
      dF(38) = (y5**3-3._ark*y5*y6**2)*y1+(y5**3-3._ark*y5*y6**2)*y2+(y5**3- &
          3._ark*y5*y6**2)*y3+(y5**3-3._ark*y5*y6**2)*y4
      dF(39) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1+(y4*y7**2+y3*y8**2)*y2+y3*y4*y9**2
      dF(40) = (y4*y7*y9+y3*y8*y9+y2*y7*y8)*y1+(-y4*y8*y9-y3*y7*y9)*y2-y3*y4*y7*y8
      dF(41) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1+((y7**2+ &
          y9**2)*y3+(y9**2+y8**2)*y4)*y2+(y8**2+y7**2)*y4*y3
      dF(42) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1+((y6**2+ &
          y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3
      dF(43) = ((sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark)*y2+(y5*y6+ &
          sqrt(3._ark)*y5**2/3._ark)*y3+(sqrt(3._ark)*y5**2/3._ark-y5*y6)*y4)*y1+ &
          ((sqrt(3._ark)*y5**2/3._ark-y5*y6)*y3+(y5*y6+sqrt(3._ark)*y5**2/3._ark)*y4)*y2+ &
          (sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark)*y4*y3
      dF(44) = (y2*y5*y9+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3+(-y5*y8/2._ark- &
          sqrt(3._ark)*y6*y8/2._ark)*y4)*y1+((y5*y8/2._ark+sqrt(3._ark)*y6*y8/2._ark)*y3+(y5*y7/ &
          2._ark-sqrt(3._ark)*y6*y7/2._ark)*y4)*y2-y3*y4*y5*y9
      dF(45) = (((y9+y7-y8)*y3+(y8+y9-y7)*y4)*y2+(y7-y9+y8)*y4*y3)*y1+(-y7-y8- &
          y9)*y4*y3*y2
      dF(46) = y1*y2*y3*y4
      dF(47) = (y3*y7+y4*y8+y2*y9)*y1**2+(y4**2*y8+y3**2*y7+y2**2*y9)*y1+(-y3*y8- &
          y4*y7)*y2**2+(-y3**2*y8-y4**2*y7)*y2-y3*y4**2*y9-y3**2*y4*y9
      dF(48) = ((-y8-y7)*y2+(-y9-y8)*y3+(-y9-y7)*y4)*y1**2+((y7+y8)*y2**2+(y9+ &
          y8)*y3**2+(y7+y9)*y4**2)*y1+((y7-y9)*y3+(y8-y9)*y4)*y2**2+((y9-y7)*y3**2+(-y8+ &
          y9)*y4**2)*y2+(y8-y7)*y4*y3**2+(y7-y8)*y4**2*y3
      dF(49) = (y2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4)*y1**2+(y2**2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y1+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+ &
          sqrt(3._ark)*y6/2._ark)*y4)*y2**2+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark+ &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y2+y3*y4**2*y5+y3**2*y4*y5
      dF(50) = ((y4+y3)*y2+y3*y4)*y1**2+((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+ &
          y3**2*y4)*y1+y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2
      dF(51) = (y4+y3+y2)*y1**3+(y4**3+y3**3+y2**3)*y1+(y4+y3)*y2**3+(y4**3+y3**3)*y2+ &
          y3*y4**3+y3**3*y4
      dF(52) = ((y9+y8)*y7+y8*y9)*y1**2+((y8-y9)*y7-y8*y9)*y2**2+((-y9-y8)*y7+ &
          y8*y9)*y3**2+((-y8+y9)*y7-y8*y9)*y4**2
      dF(53) = (y8**2+y7**2+y9**2)*y1**2+(y8**2+y7**2+y9**2)*y2**2+(y8**2+y7**2+ &
          y9**2)*y3**2+(y8**2+y7**2+y9**2)*y4**2
      dF(54) = (y6**2+y5**2)*y1**2+(y6**2+y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+ &
          y5**2)*y4**2
      dF(55) = ((y9-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y1**2+((y9+y7/2._ark+y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**2+((-y7/2._ark-y9+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6)*y3**2+((y7/2._ark-y8/2._ark-y9)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y4**2
      dF(56) = y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7
      dF(57) = (2._ark/3._ark*sqrt(3._ark)*y9**4-sqrt(3._ark)*y8**4/3._ark-sqrt(3._ark)*y7**4/ &
          3._ark)*y5+(-y8**4+y7**4)*y6
      dF(58) = sqrt(3._ark)*y5**3*y9**2+(y7**2-y8**2)*y6*y5**2+(-sqrt(3._ark)*y9**2/3._ark- &
         4._ark/3._ark*sqrt(3._ark)*y7**2-4._ark/3._ark*sqrt(3._ark)*y8**2)*y6**2*y5+(y7**2- &
         y8**2)*y6**3
      dF(59) = ((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+ &
          (sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6
      dF(60) = y5**2*y7*y8*y9+y6**2*y7*y8*y9
      dF(61) = (y8**2+y7**2+y9**2)*y5**3+(-3._ark*y8**2-3._ark*y7**2-3._ark*y9**2)*y6**2*y5
      dF(62) = -3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2
      s1 = (((y9/2._ark-y8)*y7**2+(y9**2/2._ark-y8**2)*y7+y8*y9**2/2._ark+y8**2*y9/2._ark)*y5+ &
          (-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y9/2._ark)*y6)*y1+(((y8+y9/2._ark)*y7**2+(-y9**2/2._ark+y8**2)*y7- &
          y8*y9**2/2._ark+y8**2*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/ &
          2._ark+sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y6)*y2
      dF(63) = s1+(((y8-y9/2._ark)*y7**2+(y9**2/2._ark-y8**2)*y7-y8**2*y9/2._ark-y8*y9**2/ &
          2._ark)*y5+(-sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8**2*y9/2._ark)*y6)*y3+(((-y9/2._ark-y8)*y7**2+ &
          (-y9**2/2._ark+y8**2)*y7-y8**2*y9/2._ark+y8*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9/ &
          2._ark-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6)*y4
      dF(64) = ((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(y8*y9- &
          y7*y9)*y6*y5+((sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y1+((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y2+((-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/2._ark)*y5**2+(y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark-2._ark/3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y3+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9- &
          y7*y9)*y6*y5+((-2._ark/3._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/6._ark)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y4
      dF(65) = (-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y1+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y2+(-sqrt(3._ark)*y5**2*y9**2/ &
          2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/3._ark- &
          sqrt(3._ark)*y7**2/3._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2- &
          y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/ &
          3._ark)*y6**2)*y4
      dF(66) = (((y9+y8)*y7+y8*y9)*y5**2+((y9+y8)*y7+y8*y9)*y6**2)*y1+(((y8-y9)*y7- &
          y8*y9)*y5**2+((y8-y9)*y7-y8*y9)*y6**2)*y2+(((-y9-y8)*y7+y8*y9)*y5**2+((-y9- &
          y8)*y7+y8*y9)*y6**2)*y3+(((-y8+y9)*y7-y8*y9)*y5**2+((-y8+y9)*y7-y8*y9)*y6**2)*y4
      dF(67) = ((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y1+((y8**2+y7**2+ &
          y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y2+((y8**2+y7**2+y9**2)*y5**2+(y8**2+ &
          y7**2+y9**2)*y6**2)*y3+((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y4
      dF(68) = ((sqrt(3._ark)*y7+sqrt(3._ark)*y8)*y5**3+(y8-y7)*y6*y5**2+(-5._ark/ &
          3._ark*sqrt(3._ark)*y7-5._ark/3._ark*sqrt(3._ark)*y8-8._ark/3._ark*sqrt(3._ark)*y9)*y6**2*y5+ &
          (y8-y7)*y6**3)*y1+((-sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y5**3+(y7-y8)*y6*y5**2+(-8._ark/ &
          3._ark*sqrt(3._ark)*y9+5._ark/3._ark*sqrt(3._ark)*y7+5._ark/3._ark*sqrt(3._ark)*y8)*y6**2*y5+ &
          (y7-y8)*y6**3)*y2+((sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y5**3+(-y8-y7)*y6*y5**2+(8._ark/ &
          3._ark*sqrt(3._ark)*y9+5._ark/3._ark*sqrt(3._ark)*y8-5._ark/3._ark*sqrt(3._ark)*y7)*y6**2*y5+(- &
          y8-y7)*y6**3)*y3+((sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y5**3+(y7+y8)*y6*y5**2+(5._ark/ &
          3._ark*sqrt(3._ark)*y7-5._ark/3._ark*sqrt(3._ark)*y8+8._ark/3._ark*sqrt(3._ark)*y9)*y6**2*y5+ &
          (y7+y8)*y6**3)*y4
      dF(69) = ((y9+y7+y8)*y5**3+(-3._ark*y8-3._ark*y9-3._ark*y7)*y6**2*y5)*y1+((-y7+y9- &
          y8)*y5**3+(3._ark*y7+3._ark*y8-3._ark*y9)*y6**2*y5)*y2+((-y8-y9+y7)*y5**3+(-3._ark*y7+ &
          3._ark*y8+3._ark*y9)*y6**2*y5)*y3+((-y7-y9+y8)*y5**3+(3._ark*y7+3._ark*y9- &
          3._ark*y8)*y6**2*y5)*y4
      dF(70) = (y6**4+y5**4+2._ark*y5**2*y6**2)*y1+(y6**4+y5**4+2._ark*y5**2*y6**2)*y2+ &
          (y6**4+y5**4+2._ark*y5**2*y6**2)*y3+(y6**4+y5**4+2._ark*y5**2*y6**2)*y4
      dF(71) = (y4*y7*y8*y9+y3*y7*y8*y9+y2*y7*y8*y9)*y1+(y4*y7*y8*y9+y3*y7*y8*y9)*y2+ &
          y3*y4*y7*y8*y9
      dF(72) = ((-y8**2*y9-y7**2*y9)*y2+(-y9**2-y8**2)*y7*y3+(-y8*y9**2- &
          y7**2*y8)*y4)*y1+((y7**2*y8+y8*y9**2)*y3+(y9**2+y8**2)*y7*y4)*y2+(y7**2*y9+ &
          y8**2*y9)*y4*y3
      dF(73) = (-y2*y9**3-y4*y8**3-y3*y7**3)*y1+(y4*y7**3+y3*y8**3)*y2+y3*y4*y9**3
      dF(74) = (((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark-y7**2/ &
          2._ark)*y6)*y2+(-sqrt(3._ark)*y5*y9**2/2._ark+(y9**2/2._ark+y8**2)*y6)*y3+(- &
          sqrt(3._ark)*y5*y9**2/2._ark+(-y7**2-y9**2/2._ark)*y6)*y4)*y1+((-sqrt(3._ark)*y5*y9**2/ &
          2._ark+(-y7**2-y9**2/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y9**2/2._ark+(y9**2/2._ark+ &
          y8**2)*y6)*y4)*y2+((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark- &
          y7**2/2._ark)*y6)*y4*y3
      dF(75) = (2._ark*y2*y5*y7*y8+(-y5*y8*y9+sqrt(3._ark)*y6*y8*y9)*y3+(- &
          sqrt(3._ark)*y6*y7*y9-y5*y7*y9)*y4)*y1+((y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y3+(- &
          sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y4)*y2-2._ark*y3*y4*y5*y7*y8
      dF(76) = (((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y7**2/ &
          2._ark)*y6)*y2+((-y9**2/2._ark+y8**2)*y5-sqrt(3._ark)*y6*y9**2/2._ark)*y3+((y7**2-y9**2/ &
          2._ark)*y5+sqrt(3._ark)*y6*y9**2/2._ark)*y4)*y1+(((y7**2-y9**2/2._ark)*y5+ &
          sqrt(3._ark)*y6*y9**2/2._ark)*y3+((-y9**2/2._ark+y8**2)*y5-sqrt(3._ark)*y6*y9**2/ &
          2._ark)*y4)*y2+((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y8**2/2._ark- &
          sqrt(3._ark)*y7**2/2._ark)*y6)*y4*y3
      dF(77) = (-2._ark*y2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3+(y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y4)*y1+((y5*y8**2+sqrt(3._ark)*y6*y8**2)*y3+(- &
          sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4)*y2-2._ark*y3*y4*y5*y9**2
      dF(78) = ((-sqrt(3._ark)*y6**2*y9/6._ark+sqrt(3._ark)*y5**2*y9/2._ark)*y2+(-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/3._ark)*y3+(sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y4)*y1+((- &
          sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y3+(-sqrt(3._ark)*y6**2*y7/3._ark+ &
          y5*y6*y7)*y4)*y2+(sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y4*y3
      dF(79) = ((-y5**3/3._ark+y5*y6**2)*y2+(-y5**3/3._ark+y5*y6**2)*y3+(-y5**3/3._ark+ &
          y5*y6**2)*y4)*y1+((-y5**3/3._ark+y5*y6**2)*y3+(-y5**3/3._ark+y5*y6**2)*y4)*y2+(- &
          y5**3/3._ark+y5*y6**2)*y4*y3
      dF(80) = ((-y9*y6**2-y9*y5**2)*y2+(-y5**2*y7-y7*y6**2)*y3+(-y5**2*y8- &
          y8*y6**2)*y4)*y1+((y8*y6**2+y5**2*y8)*y3+(y7*y6**2+y5**2*y7)*y4)*y2+(y9*y6**2+ &
          y9*y5**2)*y4*y3
      dF(81) = ((5._ark/9._ark*sqrt(3._ark)*y5**3+sqrt(3._ark)*y5*y6**2)*y2+(y6**3+y5**2*y6- &
          4._ark/9._ark*sqrt(3._ark)*y5**3)*y3+(-y6**3-y5**2*y6-4._ark/ &
          9._ark*sqrt(3._ark)*y5**3)*y4)*y1+((-y6**3-y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3)*y3+ &
          (y6**3+y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3)*y4)*y2+(5._ark/9._ark*sqrt(3._ark)*y5**3+ &
          sqrt(3._ark)*y5*y6**2)*y4*y3
      dF(82) = (-y3*y8*y9-y2*y7*y8-y4*y7*y9)*y1**2+(-y4**2*y7*y9-y2**2*y7*y8- &
          y3**2*y8*y9)*y1+(y3*y7*y9+y4*y8*y9)*y2**2+(y4**2*y8*y9+y3**2*y7*y9)*y2+ &
          y3**2*y4*y7*y8+y3*y4**2*y7*y8
      dF(83) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1**2+((y8**2+ &
          y7**2)*y2**2+(y9**2+y8**2)*y3**2+(y7**2+y9**2)*y4**2)*y1+((y7**2+y9**2)*y3+ &
          (y9**2+y8**2)*y4)*y2**2+((y7**2+y9**2)*y3**2+(y9**2+y8**2)*y4**2)*y2+(y8**2+ &
          y7**2)*y4*y3**2+(y8**2+y7**2)*y4**2*y3
      dF(84) = ((-y8*y9-y7*y9)*y2+(-y9-y8)*y7*y3+(-y7*y8-y8*y9)*y4)*y1**2+((y8*y9+ &
          y7*y9)*y2**2+(y9+y8)*y7*y3**2+(y8*y9+y7*y8)*y4**2)*y1+((-y7*y8+y8*y9)*y3+(-y8+ &
          y9)*y7*y4)*y2**2+((-y8*y9+y7*y8)*y3**2+(y8-y9)*y7*y4**2)*y2+(-y8*y9+ &
          y7*y9)*y4*y3**2+(y8*y9-y7*y9)*y4**2*y3
      dF(85) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1**2+(y4**2*y8**2+y3**2*y7**2+ &
          y2**2*y9**2)*y1+(y4*y7**2+y3*y8**2)*y2**2+(y4**2*y7**2+y3**2*y8**2)*y2+ &
          y3*y4**2*y9**2+y3**2*y4*y9**2
      s1 = (((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/2._ark)*y6)*y2+(- &
          sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark- &
          y7)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/ &
          2._ark)*y6)*y2**2+(sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark-y8)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y9/2._ark+(y7+y9/2._ark)*y6)*y4**2)*y1+((-sqrt(3._ark)*y5*y9/2._ark+(-y9/ &
         2._ark+y7)*y6)*y3+(-sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark-y8)*y6)*y4)*y2**2
      dF(86) = s1+((sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/2._ark)*y6)*y3**2+(sqrt(3._ark)*y5*y9/ &
          2._ark+(y8-y9/2._ark)*y6)*y4**2)*y2+((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(- &
          y8/2._ark-y7/2._ark)*y6)*y4*y3**2+((-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y7/ &
          2._ark+y8/2._ark)*y6)*y4**2*y3
      dF(87) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1**2+((y6**2+ &
          y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y1+((y6**2+y5**2)*y3+ &
          (y6**2+y5**2)*y4)*y2**2+((y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2+(y6**2+ &
          y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3
      s1 = (((-y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y2+ &
          ((y8-y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3+((-y9/2._ark+y7)*y5+sqrt(3._ark)*y6*y9/ &
          2._ark)*y4)*y1**2+(((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**2+((y9/2._ark-y8)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3**2+((-y7+y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y1+(((-y9/2._ark-y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3+ &
          ((-y9/2._ark-y8)*y5-sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**2
      dF(88) = s1+(((y7+y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3**2+((y8+y9/2._ark)*y5+ &
          sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y2+((y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3**2+((y7/2._ark-y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4**2*y3
      dF(89) = ((((-y8+y9)*y7-y8*y9)*y3+((-y9-y8)*y7+y8*y9)*y4)*y2+((y8-y9)*y7- &
          y8*y9)*y4*y3)*y1+((y9+y8)*y7+y8*y9)*y4*y3*y2
      dF(90) = (((y8**2+y7**2+y9**2)*y3+(y8**2+y7**2+y9**2)*y4)*y2+(y8**2+y7**2+ &
          y9**2)*y4*y3)*y1+(y8**2+y7**2+y9**2)*y4*y3*y2
      dF(91) = (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3)*y1+(y6**2+ &
          y5**2)*y4*y3*y2
      dF(92) = ((-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2+(-y5*y6- &
          sqrt(3._ark)*y5**2/3._ark)*y3+(y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4)*y1**2+((- &
          sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2**2+(-y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3**2+(y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4**2)*y1+((y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4)*y2**2+((y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3**2+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4**2)*y2+(-sqrt(3._ark)*y6**2/2._ark+ &
          sqrt(3._ark)*y5**2/6._ark)*y4*y3**2+(-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/ &
          6._ark)*y4**2*y3
      dF(93) = ((y3*y8+y4*y7)*y2+y3*y4*y9)*y1**2+((-y3*y7-y4*y8)*y2**2+(-y4**2*y9- &
          y3**2*y9)*y2-y3*y4**2*y7-y3**2*y4*y8)*y1+y2**2*y3*y4*y9+(y3*y4**2*y8+ &
          y3**2*y4*y7)*y2
      dF(94) = (((-sqrt(3._ark)*y5/3._ark-y6)*y3+(-sqrt(3._ark)*y5/3._ark+y6)*y4)*y2+2._ark/ &
          3._ark*sqrt(3._ark)*y3*y4*y5)*y1**2+(((-sqrt(3._ark)*y5/3._ark+y6)*y3+(-sqrt(3._ark)*y5/ &
          3._ark-y6)*y4)*y2**2+(2._ark/3._ark*sqrt(3._ark)*y4**2*y5+2._ark/ &
          3._ark*sqrt(3._ark)*y3**2*y5)*y2+(-sqrt(3._ark)*y5/3._ark-y6)*y4*y3**2+(-sqrt(3._ark)*y5/ &
          3._ark+y6)*y4**2*y3)*y1+2._ark/3._ark*sqrt(3._ark)*y2**2*y3*y4*y5+((-sqrt(3._ark)*y5/3._ark+ &
          y6)*y4*y3**2+(-sqrt(3._ark)*y5/3._ark-y6)*y4**2*y3)*y2
      dF(95) = ((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+y3**2*y4)*y1**2+((y3**2+ &
          y4**2)*y2**2+y3**2*y4**2)*y1+(y3*y4**2+y3**2*y4)*y2**2+y2*y3**2*y4**2
      dF(96) = (-y4*y8-y3*y7-y2*y9)*y1**3+(-y4**3*y8-y3**3*y7-y2**3*y9)*y1+(y3*y8+ &
          y4*y7)*y2**3+(y3**3*y8+y4**3*y7)*y2+y3**3*y4*y9+y3*y4**3*y9
      dF(97) = ((y7+y8)*y2+(y9+y8)*y3+(y7+y9)*y4)*y1**3+((-y8-y7)*y2**3+(-y9- &
          y8)*y3**3+(-y9-y7)*y4**3)*y1+((y9-y7)*y3+(-y8+y9)*y4)*y2**3+((y7-y9)*y3**3+(y8- &
          y9)*y4**3)*y2+(y7-y8)*y4*y3**3+(y8-y7)*y4**3*y3
      dF(98) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4)*y1**3+(-2._ark/3._ark*sqrt(3._ark)*y2**3*y5+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y3**3+(y6+sqrt(3._ark)*y5/3._ark)*y4**3)*y1+((y6+sqrt(3._ark)*y5/ &
          3._ark)*y3+(-y6+sqrt(3._ark)*y5/3._ark)*y4)*y2**3+((y6+sqrt(3._ark)*y5/3._ark)*y3**3+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**3)*y2-2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y5-2._ark/ &
          3._ark*sqrt(3._ark)*y3*y4**3*y5
      dF(99) = ((y4+y3)*y2+y3*y4)*y1**3+((y4+y3)*y2**3+(y4**3+y3**3)*y2+y3**3*y4+ &
          y3*y4**3)*y1+y2**3*y3*y4+(y3**3*y4+y3*y4**3)*y2
      dF(100) = (y4+y3+y2)*y1**4+(y4**4+y2**4+y3**4)*y1+(y4+y3)*y2**4+(y4**4+ &
          y3**4)*y2+y3*y4**4+y3**4*y4
      dF(101) = y2**2*y7*y8*y9+y3**2*y7*y8*y9+y1**2*y7*y8*y9+y4**2*y7*y8*y9
      dF(102) = ((-y9-y8)*y7**2+(-y9**2-y8**2)*y7-y8*y9**2-y8**2*y9)*y1**2+((y8- &
          y9)*y7**2+(y9**2+y8**2)*y7-y8**2*y9+y8*y9**2)*y2**2+((y9+y8)*y7**2+(-y9**2- &
          y8**2)*y7+y8**2*y9+y8*y9**2)*y3**2+((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+ &
          y8**2*y9)*y4**2
      dF(103) = (-y8**3-y7**3-y9**3)*y1**2+(-y9**3+y8**3+y7**3)*y2**2+(y8**3+y9**3- &
          y7**3)*y3**2+(y7**3+y9**3-y8**3)*y4**2
      dF(104) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1**2+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2**2+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3**2+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**2
      dF(105) = ((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y1**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y2**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y3**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4**2
      dF(106) = ((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5**2+(y8-y7)*y6*y5+(- &
          sqrt(3._ark)*y8/6._ark-2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y6**2)*y1**2+ &
          ((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5**2+(y7-y8)*y6*y5+(-2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y6**2)*y2**2+((- &
          sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5**2+(-y8-y7)*y6*y5+(-sqrt(3._ark)*y7/ &
          6._ark+2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/6._ark)*y6**2)*y3**2+((sqrt(3._ark)*y7/ &
          2._ark-sqrt(3._ark)*y8/2._ark)*y5**2+(y7+y8)*y6*y5+(2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark)*y6**2)*y4**2
      dF(107) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1**2+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2**2+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3**2+((-y7-y9+y8)*y5**2+ &
          (-y7-y9+y8)*y6**2)*y4**2
      dF(108) = (y5**3-3._ark*y5*y6**2)*y1**2+(y5**3-3._ark*y5*y6**2)*y2**2+(y5**3- &
          3._ark*y5*y6**2)*y3**2+(y5**3-3._ark*y5*y6**2)*y4**2
      dF(109) = (2._ark/3._ark*sqrt(3._ark)*y2*y5*y9+(-sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3+(- &
          sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4)*y1**2+(2._ark/3._ark*sqrt(3._ark)*y2**2*y5*y9+(- &
          sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3**2+(-sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4**2)*y1+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4)*y2**2+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3**2+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4**2)*y2- &
          2._ark/3._ark*sqrt(3._ark)*y3*y4**2*y5*y9-2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5*y9
      dF(110) = (-y3**2*y7-y4**2*y8-y2**2*y9)*y1**2+(y4**2*y7+y3**2*y8)*y2**2+ &
          y3**2*y4**2*y9
      dF(111) = (-2._ark/3._ark*sqrt(3._ark)*y2**2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3**2+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**2)*y1**2+((y6+sqrt(3._ark)*y5/3._ark)*y3**2+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**2)*y2**2-2._ark/3._ark*sqrt(3._ark)*y3**2*y4**2*y5
      dF(112) = ((y9+y8)*y7+y8*y9)*y1**3+((y8-y9)*y7-y8*y9)*y2**3+((-y9-y8)*y7+ &
          y8*y9)*y3**3+((-y8+y9)*y7-y8*y9)*y4**3
      dF(113) = (y8**2+y7**2+y9**2)*y1**3+(y8**2+y7**2+y9**2)*y2**3+(y8**2+y7**2+ &
          y9**2)*y3**3+(y8**2+y7**2+y9**2)*y4**3
      dF(114) = ((-2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5+ &
          (y8-y7)*y6)*y1**3+((-sqrt(3._ark)*y7/3._ark-2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y8/ &
          3._ark)*y5+(y7-y8)*y6)*y2**3+((sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(-y8-y7)*y6)*y3**3+((2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y5+(y7+y8)*y6)*y4**3
      dF(115) = (y6**2+y5**2)*y1**3+(y6**2+y5**2)*y2**3+(y6**2+y5**2)*y3**3+(y6**2+ &
          y5**2)*y4**3
      dF(116) = (y3**2+y2**2+y4**2)*y1**3+(y4**3+y3**3+y2**3)*y1**2+(y3**2+ &
          y4**2)*y2**3+(y4**3+y3**3)*y2**2+y3**2*y4**3+y3**3*y4**2
      dF(117) = (-y7-y8-y9)*y1**4+(y7-y9+y8)*y2**4+(y8+y9-y7)*y3**4+(y9+y7-y8)*y4**4
      dF(118) = y4**5+y3**5+y2**5+y1**5
      dF(119) = (y7**2*y8*y9+(y8*y9**2+y8**2*y9)*y7)*y1+(-y7**2*y8*y9+(y8*y9**2- &
          y8**2*y9)*y7)*y2+(y7**2*y8*y9+(-y8**2*y9-y8*y9**2)*y7)*y3+(-y7**2*y8*y9+(- &
          y8*y9**2+y8**2*y9)*y7)*y4
      dF(120) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y1+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y2+((y9**2+y8**2)*y7**2+y8**2*y9**2)*y3+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y4
      dF(121) = ((y9+y8)*y7**3+(y9**3+y8**3)*y7+y8**3*y9+y8*y9**3)*y1+((y8-y9)*y7**3+ &
          (y8**3-y9**3)*y7-y8*y9**3-y8**3*y9)*y2+((-y9-y8)*y7**3+(-y9**3-y8**3)*y7+ &
          y8**3*y9+y8*y9**3)*y3+((-y8+y9)*y7**3+(y9**3-y8**3)*y7-y8*y9**3-y8**3*y9)*y4
      dF(122) = (y9**4+y8**4+y7**4)*y1+(y9**4+y8**4+y7**4)*y2+(y9**4+y8**4+y7**4)*y3+ &
          (y9**4+y8**4+y7**4)*y4
      s1 = ((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5+((y8+y9/2._ark)*y7**2+(-y9**2/2._ark-y8**2)*y7- &
          y8**2*y9/2._ark+y8*y9**2/2._ark)*y6)*y1+((-sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y5+ &
          ((y9/2._ark-y8)*y7**2+(y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/2._ark)*y6)*y2
      dF(123) = s1+((sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5+((-y9/2._ark-y8)*y7**2+(- &
          y9**2/2._ark-y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y6)*y3+((sqrt(3._ark)*y8**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y5+((y8-y9/2._ark)*y7**2+(y9**2/2._ark+y8**2)*y7+y8**2*y9/2._ark+y8*y9**2/ &
          2._ark)*y6)*y4
      dF(124) = ((-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3)*y5+(-y8**3+y7**3)*y6)*y1+((sqrt(3._ark)*y7**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3+y8**3)*y6)*y2+((-2._ark/ &
          3._ark*sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark)*y5+(y8**3+ &
          y7**3)*y6)*y3+((-2._ark/3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3-y8**3)*y6)*y4
      dF(125) = ((((-sqrt(3._ark)*y8/3._ark-2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          3._ark)*y5+(-y8-y7)*y6)*y3+((-sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(y7+y8)*y6)*y4)*y2+((sqrt(3._ark)*y8/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/3._ark)*y5+(y8-y7)*y6)*y4*y3)*y1+((2._ark/ &
          3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5+(y7- &
          y8)*y6)*y4*y3*y2
      dF(126) = (((y7+y9)*y3+(y9+y8)*y4)*y2+(y7+y8)*y4*y3)*y1**2+(((-y8+y9)*y3+(y9- &
          y7)*y4)*y2**2+((y7-y8)*y3**2+(y8-y7)*y4**2)*y2+(y7-y9)*y4*y3**2+(y8- &
          y9)*y4**2*y3)*y1+(-y8-y7)*y4*y3*y2**2+((-y9-y8)*y4*y3**2+(-y9-y7)*y4**2*y3)*y2
      dF(127) = y1**2*y2*y3*y4+(y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2)*y1
      dF(128) = (y9**2+y8**2)*y7**4+(y8**4+y9**4)*y7**2+y8**4*y9**2+y8**2*y9**4
      dF(129) = y7**6+y9**6+y8**6
      dF(130) = y7**2*y8**2*y9**2
      dF(131) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y5**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y6**2
      dF(132) = -6._ark*y5**2*y6**4+9._ark*y5**4*y6**2+y6**6
      dF(133) = (y7**3*y8*y9+(-2._ark*y8*y9**3+y8**3*y9)*y7)*y5+(sqrt(3._ark)*y7*y8**3*y9- &
          sqrt(3._ark)*y7**3*y8*y9)*y6
      dF(134) = ((2._ark/3._ark*sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2/6._ark)*y7**2+ &
          sqrt(3._ark)*y8**2*y9**2/6._ark)*y5**2+(-y8**2*y9**2+y7**2*y9**2)*y6*y5+ &
          (sqrt(3._ark)*y7**2*y9**2/2._ark+sqrt(3._ark)*y8**2*y9**2/2._ark)*y6**2
      dF(135) = -sqrt(3._ark)*y5**2*y9**4/2._ark+(-y8**4+y7**4)*y6*y5+(-sqrt(3._ark)*y7**4/ &
          3._ark-sqrt(3._ark)*y8**4/3._ark+sqrt(3._ark)*y9**4/6._ark)*y6**2
      dF(136) = -y5**3*y7*y8*y9/3._ark+y5*y6**2*y7*y8*y9
      dF(137) = (y9**4+y8**4+y7**4)*y5**2+(y9**4+y8**4+y7**4)*y6**2
      dF(138) = 3._ark/4._ark*y5**4*y9**2+(-y9**2/2._ark+y7**2+y8**2)*y6**2*y5**2+(-2._ark/ &
          3._ark*sqrt(3._ark)*y7**2+2._ark/3._ark*sqrt(3._ark)*y8**2)*y6**3*y5+(y8**2/3._ark+y9**2/ &
          12._ark+y7**2/3._ark)*y6**4
      dF(139) = -sqrt(3._ark)*y5**4*y9**2/4._ark+(y7**2-y8**2)*y6*y5**3- &
          sqrt(3._ark)*y5**2*y6**2*y9**2/2._ark+(-y8**2/3._ark+y7**2/3._ark)*y6**3*y5+(-2._ark/ &
          9._ark*sqrt(3._ark)*y8**2+7._ark/36._ark*sqrt(3._ark)*y9**2-2._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y6**4
      dF(140) = (-y9**2/2._ark+y7**2+y8**2)*y5**4+3._ark*y5**2*y6**2*y9**2+(4._ark/ &
          3._ark*sqrt(3._ark)*y7**2-4._ark/3._ark*sqrt(3._ark)*y8**2)*y6**3*y5+(y8**2/3._ark+5._ark/ &
          6._ark*y9**2+y7**2/3._ark)*y6**4
      dF(141) = 9._ark*y5**2*y6**4+y5**6-6._ark*y5**4*y6**2
      dF(142) = ((y9**2+y8**2)*y7**3+(y9**3+y8**3)*y7**2+y8**3*y9**2+y8**2*y9**3)*y1+ &
          ((-y9**2-y8**2)*y7**3+(y9**3-y8**3)*y7**2+y8**2*y9**3-y8**3*y9**2)*y2+((y9**2+ &
          y8**2)*y7**3+(-y9**3-y8**3)*y7**2-y8**2*y9**3-y8**3*y9**2)*y3+((-y9**2- &
          y8**2)*y7**3+(y8**3-y9**3)*y7**2+y8**3*y9**2-y8**2*y9**3)*y4
      dF(143) = ((y8*y9**2+y8**2*y9)*y7**2+y7*y8**2*y9**2)*y1+((-y8*y9**2+ &
          y8**2*y9)*y7**2-y7*y8**2*y9**2)*y2+((-y8**2*y9-y8*y9**2)*y7**2+ &
          y7*y8**2*y9**2)*y3+((y8*y9**2-y8**2*y9)*y7**2-y7*y8**2*y9**2)*y4
      dF(144) = (y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7)*y1+(y7**3*y8*y9+(y8**3*y9+ &
          y8*y9**3)*y7)*y2+(y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7)*y3+(y7**3*y8*y9+(y8**3*y9+ &
          y8*y9**3)*y7)*y4
      dF(145) = ((y9+y8)*y7**4+(y8**4+y9**4)*y7+y8**4*y9+y8*y9**4)*y1+((-y8+y9)*y7**4+ &
          (-y8**4-y9**4)*y7+y8**4*y9-y8*y9**4)*y2+((-y9-y8)*y7**4+(y8**4+y9**4)*y7- &
          y8**4*y9-y8*y9**4)*y3+((y8-y9)*y7**4+(-y8**4-y9**4)*y7+y8*y9**4-y8**4*y9)*y4
      dF(146) = (-y7**5-y8**5-y9**5)*y1+(y7**5-y9**5+y8**5)*y2+(-y7**5+y8**5+ &
          y9**5)*y3+(y7**5+y9**5-y8**5)*y4
      s1 = ((-sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y8*y9**3/ &
          2._ark-sqrt(3._ark)*y8**3*y9/2._ark)*y5+((y8+y9/2._ark)*y7**3+(-y9**3/2._ark-y8**3)*y7- &
          y8**3*y9/2._ark+y8*y9**3/2._ark)*y6)*y1+((sqrt(3._ark)*y8**3*y9/2._ark+ &
          sqrt(3._ark)*y7**3*y9/2._ark-sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8*y9**3/2._ark)*y5+ &
          ((y8-y9/2._ark)*y7**3+(y9**3/2._ark-y8**3)*y7+y8**3*y9/2._ark-y8*y9**3/2._ark)*y6)*y2
      dF(147) = s1+((-sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8**3*y9/2._ark+ &
          sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y8*y9**3/2._ark)*y5+((-y9/2._ark-y8)*y7**3+ &
          (y8**3+y9**3/2._ark)*y7-y8**3*y9/2._ark+y8*y9**3/2._ark)*y6)*y3+((sqrt(3._ark)*y8**3*y9/ &
          2._ark-sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8*y9**3/ &
          2._ark)*y5+((y9/2._ark-y8)*y7**3+(-y9**3/2._ark+y8**3)*y7+y8**3*y9/2._ark-y8*y9**3/ &
          2._ark)*y6)*y4
      dF(148) = ((y7**2*y8*y9/2._ark+(y8**2*y9/2._ark-y8*y9**2)*y7)*y5+ &
          (sqrt(3._ark)*y7*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y1+((-y7**2*y8*y9/ &
          2._ark+(-y8**2*y9/2._ark-y8*y9**2)*y7)*y5+(-sqrt(3._ark)*y7*y8**2*y9/2._ark+ &
          sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y2+((y7**2*y8*y9/2._ark+(-y8**2*y9/2._ark+ &
          y8*y9**2)*y7)*y5+(-sqrt(3._ark)*y7*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y8*y9/ &
          2._ark)*y6)*y3+((-y7**2*y8*y9/2._ark+(y8*y9**2+y8**2*y9/2._ark)*y7)*y5+ &
          (sqrt(3._ark)*y7*y8**2*y9/2._ark+sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y4
      dF(149) = ((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y1+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y2+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y3+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y4
      dF(150) = ((y9+y8)*y7+y8*y9)*y1**4+((y8-y9)*y7-y8*y9)*y2**4+((-y9-y8)*y7+ &
          y8*y9)*y3**4+((-y8+y9)*y7-y8*y9)*y4**4
      dF(151) = (((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+ &
          (sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y1+(((-y9**2/2._ark+ &
          y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark- &
          sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y2+(((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/ &
          2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y3+(((- &
          y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark- &
          sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y4
      s1 = (((y9/2._ark-y8)*y7**3+(y9**3/2._ark-y8**3)*y7+y8**3*y9/2._ark+y8*y9**3/2._ark)*y5+ &
          (-sqrt(3._ark)*y8*y9**3/2._ark+sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y7**3*y9/2._ark- &
          sqrt(3._ark)*y8**3*y9/2._ark)*y6)*y1+(((-y9/2._ark-y8)*y7**3+(-y9**3/2._ark-y8**3)*y7- &
          y8**3*y9/2._ark-y8*y9**3/2._ark)*y5+(-sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y8**3*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**3/2._ark-sqrt(3._ark)*y7**3*y9/2._ark)*y6)*y2
      dF(152) = s1+(((y8-y9/2._ark)*y7**3+(-y9**3/2._ark+y8**3)*y7+y8**3*y9/2._ark+y8*y9**3/ &
          2._ark)*y5+(-sqrt(3._ark)*y8*y9**3/2._ark-sqrt(3._ark)*y7*y9**3/2._ark- &
          sqrt(3._ark)*y8**3*y9/2._ark-sqrt(3._ark)*y7**3*y9/2._ark)*y6)*y3+(((y8+y9/2._ark)*y7**3+ &
          (y8**3+y9**3/2._ark)*y7-y8**3*y9/2._ark-y8*y9**3/2._ark)*y5+(sqrt(3._ark)*y7*y9**3/2._ark+ &
          sqrt(3._ark)*y8**3*y9/2._ark+sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y8*y9**3/ &
          2._ark)*y6)*y4
      s1 = (((-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y9**2/3._ark)*y7+sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/ &
          6._ark)*y5**2+(-y7**2*y8+(y9**2+y8**2)*y7-y8*y9**2)*y6*y5+(sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8**2*y9/2._ark)*y6**2)*y1+(((-sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/ &
          3._ark)*y7**2+(-sqrt(3._ark)*y9**2/3._ark-sqrt(3._ark)*y8**2/3._ark)*y7- &
          sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/6._ark)*y5**2+(y7**2*y8+(-y9**2- &
          y8**2)*y7+y8*y9**2)*y6*y5+(sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8**2*y9/ &
          2._ark)*y6**2)*y2
      dF(153) = s1+(((sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark)*y7**2+(sqrt(3._ark)*y8**2/ &
          3._ark+sqrt(3._ark)*y9**2/3._ark)*y7+sqrt(3._ark)*y8**2*y9/6._ark-sqrt(3._ark)*y8*y9**2/ &
          3._ark)*y5**2+(y7**2*y8+(y9**2+y8**2)*y7+y8*y9**2)*y6*y5+(-sqrt(3._ark)*y8**2*y9/ &
          2._ark-sqrt(3._ark)*y7**2*y9/2._ark)*y6**2)*y3+(((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y9/ &
          6._ark)*y7**2+(-sqrt(3._ark)*y9**2/3._ark-sqrt(3._ark)*y8**2/3._ark)*y7+ &
          sqrt(3._ark)*y8*y9**2/3._ark+sqrt(3._ark)*y8**2*y9/6._ark)*y5**2+(-y7**2*y8+(-y9**2- &
          y8**2)*y7-y8*y9**2)*y6*y5+(-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6**2)*y4
      s1 = ((sqrt(3._ark)*y8/6._ark-5._ark/24._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/6._ark)*y5**4+(- &
          sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark+7._ark/12._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(y7- &
          y8)*y6**3*y5+sqrt(3._ark)*y6**4*y9/8._ark)*y1+((-sqrt(3._ark)*y8/6._ark-5._ark/ &
          24._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y5**4+(7._ark/12._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/6._ark)*y6**2*y5**2+(y8-y7)*y6**3*y5+ &
          sqrt(3._ark)*y6**4*y9/8._ark)*y2
      dF(154) = s1+((-sqrt(3._ark)*y8/6._ark+5._ark/24._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          6._ark)*y5**4+(-7._ark/12._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/ &
          6._ark)*y6**2*y5**2+(y7+y8)*y6**3*y5-sqrt(3._ark)*y6**4*y9/8._ark)*y3+((sqrt(3._ark)*y8/ &
          6._ark-sqrt(3._ark)*y7/6._ark+5._ark/24._ark*sqrt(3._ark)*y9)*y5**4+(sqrt(3._ark)*y7/6._ark- &
          sqrt(3._ark)*y8/6._ark-7._ark/12._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(-y8-y7)*y6**3*y5- &
          sqrt(3._ark)*y6**4*y9/8._ark)*y4
      dF(155) = (((-y9-y8)*y7-y8*y9)*y5**2+((-y9-y8)*y7-y8*y9)*y6**2)*y1**2+(((-y8+ &
          y9)*y7+y8*y9)*y5**2+((-y8+y9)*y7+y8*y9)*y6**2)*y2**2+(((y9+y8)*y7-y8*y9)*y5**2+ &
          ((y9+y8)*y7-y8*y9)*y6**2)*y3**2+(((y8-y9)*y7+y8*y9)*y5**2+((y8-y9)*y7+ &
          y8*y9)*y6**2)*y4**2
      dF(156) = (y6**4+y5**4+2._ark*y5**2*y6**2)*y1**2+(y6**4+y5**4+ &
          2._ark*y5**2*y6**2)*y2**2+(y6**4+y5**4+2._ark*y5**2*y6**2)*y3**2+(y6**4+y5**4+ &
          2._ark*y5**2*y6**2)*y4**2
      dF(157) = ((-y7**3-y8**3)*y2+(-y9**3-y8**3)*y3+(-y7**3-y9**3)*y4)*y1**2+((y8**3+ &
          y7**3)*y2**2+(y9**3+y8**3)*y3**2+(y9**3+y7**3)*y4**2)*y1+((-y9**3+y7**3)*y3+ &
          (y8**3-y9**3)*y4)*y2**2+((-y7**3+y9**3)*y3**2+(y9**3-y8**3)*y4**2)*y2+(-y7**3+ &
          y8**3)*y4*y3**2+(-y8**3+y7**3)*y4**2*y3
      dF(158) = ((y8*y9**2+y7*y9**2)*y2+(y9+y8)*y7**2*y3+(y8**2*y9+ &
          y7*y8**2)*y4)*y1**2+((-y8*y9**2-y7*y9**2)*y2**2+(-y9-y8)*y7**2*y3**2+(-y8**2*y9- &
          y7*y8**2)*y4**2)*y1+((-y7*y8**2+y8**2*y9)*y3+(-y8+y9)*y7**2*y4)*y2**2+((- &
          y8**2*y9+y7*y8**2)*y3**2+(y8-y9)*y7**2*y4**2)*y2+(-y8*y9**2+y7*y9**2)*y4*y3**2+ &
          (-y7*y9**2+y8*y9**2)*y4**2*y3
      dF(159) = (y4*y7*y8*y9+y3*y7*y8*y9+y2*y7*y8*y9)*y1**2+(y4**2*y7*y8*y9+ &
          y3**2*y7*y8*y9+y2**2*y7*y8*y9)*y1+(y4*y7*y8*y9+y3*y7*y8*y9)*y2**2+ &
          (y4**2*y7*y8*y9+y3**2*y7*y8*y9)*y2+y3**2*y4*y7*y8*y9+y3*y4**2*y7*y8*y9
      dF(160) = ((y7*y8**2+y7**2*y8)*y2+(y8*y9**2+y8**2*y9)*y3+(y7*y9**2+ &
          y7**2*y9)*y4)*y1**2+((-y7**2*y8-y7*y8**2)*y2**2+(-y8**2*y9-y8*y9**2)*y3**2+(- &
          y7*y9**2-y7**2*y9)*y4**2)*y1+((-y7*y9**2+y7**2*y9)*y3+(-y8*y9**2+ &
          y8**2*y9)*y4)*y2**2+((y7*y9**2-y7**2*y9)*y3**2+(y8*y9**2-y8**2*y9)*y4**2)*y2+ &
          (y7*y8**2-y7**2*y8)*y4*y3**2+(-y7*y8**2+y7**2*y8)*y4**2*y3
      dF(161) = ((-y8**2*y9-y7**2*y9)*y2+(-y9**2-y8**2)*y7*y3+(-y8*y9**2- &
          y7**2*y8)*y4)*y1**2+((-y8**2*y9-y7**2*y9)*y2**2+(-y9**2-y8**2)*y7*y3**2+(- &
          y8*y9**2-y7**2*y8)*y4**2)*y1+((y7**2*y8+y8*y9**2)*y3+(y9**2+y8**2)*y7*y4)*y2**2+ &
          ((y7**2*y8+y8*y9**2)*y3**2+(y9**2+y8**2)*y7*y4**2)*y2+(y7**2*y9+ &
          y8**2*y9)*y4*y3**2+(y7**2*y9+y8**2*y9)*y4**2*y3
      dF(162) = (y4**2*y8**2+y3**2*y7**2+y2**2*y9**2)*y1**2+(y4**2*y7**2+ &
          y3**2*y8**2)*y2**2+y3**2*y4**2*y9**2
      dF(163) = ((y8**2+y7**2)*y2**2+(y9**2+y8**2)*y3**2+(y7**2+y9**2)*y4**2)*y1**2+ &
          ((y7**2+y9**2)*y3**2+(y9**2+y8**2)*y4**2)*y2**2+(y8**2+y7**2)*y4**2*y3**2
      dF(164) = ((-y9-y8)*y7**2+(-y9**2-y8**2)*y7-y8*y9**2-y8**2*y9)*y1**3+((y8- &
          y9)*y7**2+(y9**2+y8**2)*y7-y8**2*y9+y8*y9**2)*y2**3+((y9+y8)*y7**2+(-y9**2- &
          y8**2)*y7+y8**2*y9+y8*y9**2)*y3**3+((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+ &
          y8**2*y9)*y4**3
      dF(165) = (-y8**3-y7**3-y9**3)*y1**3+(-y9**3+y8**3+y7**3)*y2**3+(y8**3+y9**3- &
          y7**3)*y3**3+(y7**3+y9**3-y8**3)*y4**3
      dF(166) = y4**3*y7*y8*y9+y2**3*y7*y8*y9+y1**3*y7*y8*y9+y3**3*y7*y8*y9
      dF(167) = ((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y1**3+((sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y7**2/3._ark-2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y2**3+ &
          ((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(- &
          y7**2+y8**2)*y6)*y3**3+((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y4**3
      dF(168) = ((2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5+ &
          (y7-y8)*y6)*y1**4+((sqrt(3._ark)*y8/3._ark+2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          3._ark)*y5+(y8-y7)*y6)*y2**4+((-sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(y7+y8)*y6)*y3**4+((-sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/3._ark)*y5+(-y8-y7)*y6)*y4**4
      s1 = (((-sqrt(3._ark)*y9/3._ark-sqrt(3._ark)*y8/3._ark)*y7**2+(-sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y9**2/6._ark)*y7+sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y8**2*y9/ &
          3._ark)*y5**2+((-y9-y8)*y7**2+y7*y8**2+y8**2*y9)*y6*y5+(-sqrt(3._ark)*y8*y9**2/2._ark- &
          sqrt(3._ark)*y7*y9**2/2._ark)*y6**2)*y1+(((-sqrt(3._ark)*y9/3._ark+sqrt(3._ark)*y8/ &
          3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7- &
          sqrt(3._ark)*y8**2*y9/3._ark-sqrt(3._ark)*y8*y9**2/6._ark)*y5**2+((y8-y9)*y7**2- &
          y7*y8**2+y8**2*y9)*y6*y5+(sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/ &
          2._ark)*y6**2)*y2
      dF(169) = s1+(((sqrt(3._ark)*y9/3._ark+sqrt(3._ark)*y8/3._ark)*y7**2+(-sqrt(3._ark)*y8**2/ &
          3._ark+sqrt(3._ark)*y9**2/6._ark)*y7-sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y8**2*y9/ &
          3._ark)*y5**2+((y9+y8)*y7**2+y7*y8**2-y8**2*y9)*y6*y5+(sqrt(3._ark)*y8*y9**2/2._ark- &
          sqrt(3._ark)*y7*y9**2/2._ark)*y6**2)*y3+(((-sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y9/ &
          3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7+ &
          sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y8**2*y9/3._ark)*y5**2+((-y8+y9)*y7**2- &
          y7*y8**2-y8**2*y9)*y6*y5+(-sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7*y9**2/ &
          2._ark)*y6**2)*y4
      dF(170) = ((-sqrt(3._ark)*y9**3/6._ark+sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/ &
          3._ark)*y5**2+(-y8**3+y7**3)*y6*y5+sqrt(3._ark)*y6**2*y9**3/2._ark)*y1+((- &
          sqrt(3._ark)*y9**3/6._ark-sqrt(3._ark)*y8**3/3._ark-sqrt(3._ark)*y7**3/3._ark)*y5**2+(- &
          y7**3+y8**3)*y6*y5+sqrt(3._ark)*y6**2*y9**3/2._ark)*y2+((sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark+sqrt(3._ark)*y9**3/6._ark)*y5**2+(y8**3+y7**3)*y6*y5- &
          sqrt(3._ark)*y6**2*y9**3/2._ark)*y3+((-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark+ &
          sqrt(3._ark)*y9**3/6._ark)*y5**2+(-y7**3-y8**3)*y6*y5-sqrt(3._ark)*y6**2*y9**3/ &
          2._ark)*y4
      dF(171) = ((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y1+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y2+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y3+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y4
      dF(172) = (((y8/3._ark+y9/3._ark)*y7+y8*y9/3._ark)*y5**3+((-y9-y8)*y7- &
          y8*y9)*y6**2*y5)*y1+(((y8/3._ark-y9/3._ark)*y7-y8*y9/3._ark)*y5**3+((-y8+y9)*y7+ &
          y8*y9)*y6**2*y5)*y2+(((-y8/3._ark-y9/3._ark)*y7+y8*y9/3._ark)*y5**3+((y9+y8)*y7- &
          y8*y9)*y6**2*y5)*y3+(((y9/3._ark-y8/3._ark)*y7-y8*y9/3._ark)*y5**3+((y8-y9)*y7+ &
          y8*y9)*y6**2*y5)*y4
      dF(173) = ((y7**2*y9**2+y8**2*y9**2)*y2+(y9**2+y8**2)*y7**2*y3+(y8**2*y9**2+ &
          y7**2*y8**2)*y4)*y1+((y8**2*y9**2+y7**2*y8**2)*y3+(y9**2+y8**2)*y7**2*y4)*y2+ &
          (y7**2*y9**2+y8**2*y9**2)*y4*y3
      dF(174) = (((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(- &
          sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y2+((-y9**2-y8**2)*y6*y5+ &
          (sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y3+((y7**2+y9**2)*y6*y5+ &
          (sqrt(3._ark)*y7**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4)*y1+(((y7**2+ &
          y9**2)*y6*y5+(sqrt(3._ark)*y7**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y3+((-y9**2- &
          y8**2)*y6*y5+(sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4)*y2+ &
          ((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(-sqrt(3._ark)*y8**2/6._ark- &
          sqrt(3._ark)*y7**2/6._ark)*y6**2)*y4*y3
      dF(175) = (((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y2+((y9**2+y8**2)*y5**2+ &
          (y9**2+y8**2)*y6**2)*y3+((y7**2+y9**2)*y5**2+(y7**2+y9**2)*y6**2)*y4)*y1+ &
          (((y7**2+y9**2)*y5**2+(y7**2+y9**2)*y6**2)*y3+((y9**2+y8**2)*y5**2+(y9**2+ &
          y8**2)*y6**2)*y4)*y2+((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y4*y3
      s1 = (((-sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y2+ &
          (sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y3+(sqrt(3._ark)*y5*y7**2+ &
          (2._ark*y9**2+y7**2)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y5+(- &
          y7**2+y8**2)*y6)*y2**2+(sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y7**2+(2._ark*y9**2+y7**2)*y6)*y4**2)*y1+((sqrt(3._ark)*y5*y7**2+ &
          (2._ark*y9**2+y7**2)*y6)*y3+(sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y4)*y2**2
      dF(176) = s1+((sqrt(3._ark)*y5*y7**2+(2._ark*y9**2+y7**2)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y4**2)*y2+((-sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y4*y3**2+((-sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y4**2*y3
      dF(177) = (y4+y3+y2)*y1**5+(y2**5+y4**5+y3**5)*y1+(y4+y3)*y2**5+(y4**5+ &
          y3**5)*y2+y3**5*y4+y3*y4**5
      s2 = (((-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5**2+(y7-y8)*y6*y5)*y2+ &
          ((sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9-sqrt(3._ark)*y6**2*y8/ &
          2._ark)*y3+((sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y9/3._ark)*y5**2-y5*y6*y9- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4)*y1**2
      s3 = (((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5**2+(y8-y7)*y6*y5)*y2**2+((- &
          sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2-y5*y6*y9+sqrt(3._ark)*y6**2*y8/ &
          2._ark)*y3**2+((-sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9+ &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4**2)*y1+(((-sqrt(3._ark)*y9/3._ark-sqrt(3._ark)*y7/ &
          6._ark)*y5**2-y5*y6*y9+sqrt(3._ark)*y6**2*y7/2._ark)*y3+((-sqrt(3._ark)*y9/3._ark- &
          sqrt(3._ark)*y8/6._ark)*y5**2+y5*y6*y9+sqrt(3._ark)*y6**2*y8/2._ark)*y4)*y2**2
      s1 = s2+s3
      dF(178) = s1+(((sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y3**2+((sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2- &
          y5*y6*y9-sqrt(3._ark)*y6**2*y8/2._ark)*y4**2)*y2+((sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/ &
          3._ark)*y5**2+(y7+y8)*y6*y5)*y4*y3**2+((sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/ &
          3._ark)*y5**2+(-y8-y7)*y6*y5)*y4**2*y3
      s2 = (((-sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark)*y5**2+(y7-y8)*y6*y5+(- &
          sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6**2)*y2+((-2._ark/3._ark*sqrt(3._ark)*y9- &
          sqrt(3._ark)*y8/6._ark)*y5**2-y5*y6*y8-sqrt(3._ark)*y6**2*y8/2._ark)*y3+((-2._ark/ &
          3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y5**2+y5*y6*y7-sqrt(3._ark)*y6**2*y7/ &
          2._ark)*y4)*y1**2
      s3 = (((sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y5**2+(y8-y7)*y6*y5+ &
          (sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6**2)*y2**2+((2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y8/6._ark)*y5**2+y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark)*y3**2+ &
          ((sqrt(3._ark)*y7/6._ark+2._ark/3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4**2)*y1+(((sqrt(3._ark)*y7/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y7+sqrt(3._ark)*y6**2*y7/2._ark)*y3+((sqrt(3._ark)*y8/ &
          6._ark-2._ark/3._ark*sqrt(3._ark)*y9)*y5**2+y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark)*y4)*y2**2
      s1 = s2+s3
      dF(179) = s1+(((-sqrt(3._ark)*y7/6._ark+2._ark/3._ark*sqrt(3._ark)*y9)*y5**2+y5*y6*y7- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y3**2+((-sqrt(3._ark)*y8/6._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y8-sqrt(3._ark)*y6**2*y8/2._ark)*y4**2)*y2+((- &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/6._ark)*y5**2+(y7+y8)*y6*y5+(-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6**2)*y4*y3**2+((sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/ &
          6._ark)*y5**2+(-y8-y7)*y6*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6**2)*y4**2*y3
      dF(180) = ((y4*y7**2+y3*y8**2)*y2+y3*y4*y9**2)*y1**2+((y4*y8**2+y3*y7**2)*y2**2+ &
          (y4**2*y9**2+y3**2*y9**2)*y2+y3**2*y4*y8**2+y3*y4**2*y7**2)*y1+ &
          y2**2*y3*y4*y9**2+(y3*y4**2*y8**2+y3**2*y4*y7**2)*y2
      dF(181) = (((-y7*y8-y8*y9)*y3+(-y9-y8)*y7*y4)*y2+(-y8*y9-y7*y9)*y4*y3)*y1**2+ &
          (((-y8+y9)*y7*y3+(-y7*y8+y8*y9)*y4)*y2**2+((-y8*y9+y7*y9)*y3**2+(y8*y9- &
          y7*y9)*y4**2)*y2+(-y8*y9+y7*y8)*y4*y3**2+(y8-y9)*y7*y4**2*y3)*y1+(y8*y9+ &
          y7*y9)*y4*y3*y2**2+((y9+y8)*y7*y4*y3**2+(y8*y9+y7*y8)*y4**2*y3)*y2
      dF(182) = (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y1+(y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y2+ &
          (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y3+(y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y4
      dF(183) = (((y9+y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y5**2+((y9+ &
          y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y6**2)*y1+(((-y8+y9)*y7**2+(- &
          y9**2-y8**2)*y7+y8**2*y9-y8*y9**2)*y5**2+((-y8+y9)*y7**2+(-y9**2-y8**2)*y7+ &
          y8**2*y9-y8*y9**2)*y6**2)*y2+(((-y9-y8)*y7**2+(y9**2+y8**2)*y7-y8**2*y9- &
          y8*y9**2)*y5**2+((-y9-y8)*y7**2+(y9**2+y8**2)*y7-y8**2*y9-y8*y9**2)*y6**2)*y3+ &
          (((y8-y9)*y7**2+(-y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y5**2+((y8-y9)*y7**2+(- &
          y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y6**2)*y4
      dF(184) = ((y9**3+y8**3+y7**3)*y5**2+(y9**3+y8**3+y7**3)*y6**2)*y1+((y9**3- &
          y8**3-y7**3)*y5**2+(y9**3-y8**3-y7**3)*y6**2)*y2+((-y8**3-y9**3+y7**3)*y5**2+(- &
          y8**3-y9**3+y7**3)*y6**2)*y3+((-y7**3-y9**3+y8**3)*y5**2+(-y7**3-y9**3+ &
          y8**3)*y6**2)*y4
      dF(185) = (((4._ark/9._ark*sqrt(3._ark)*y9-5._ark/9._ark*sqrt(3._ark)*y8)*y7+4._ark/ &
          9._ark*sqrt(3._ark)*y8*y9)*y5**3+(-y8*y9+y7*y9)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y7*y8+ &
          (-y8*y9+y7*y9)*y6**3)*y1+(((-5._ark/9._ark*sqrt(3._ark)*y8-4._ark/ &
          9._ark*sqrt(3._ark)*y9)*y7-4._ark/9._ark*sqrt(3._ark)*y8*y9)*y5**3+(y8*y9-y7*y9)*y6*y5**2- &
          sqrt(3._ark)*y5*y6**2*y7*y8+(y8*y9-y7*y9)*y6**3)*y2+(((-4._ark/9._ark*sqrt(3._ark)*y9+ &
          5._ark/9._ark*sqrt(3._ark)*y8)*y7+4._ark/9._ark*sqrt(3._ark)*y8*y9)*y5**3+(-y8*y9- &
          y7*y9)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y7*y8+(-y8*y9-y7*y9)*y6**3)*y3+(((4._ark/ &
          9._ark*sqrt(3._ark)*y9+5._ark/9._ark*sqrt(3._ark)*y8)*y7-4._ark/ &
          9._ark*sqrt(3._ark)*y8*y9)*y5**3+(y8*y9+y7*y9)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y7*y8+ &
          (y8*y9+y7*y9)*y6**3)*y4
      dF(186) = ((-4._ark/9._ark*sqrt(3._ark)*y8**2+5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+ &
          (y7**2-y8**2)*y6**3)*y1+((-4._ark/9._ark*sqrt(3._ark)*y8**2+5._ark/ &
          9._ark*sqrt(3._ark)*y9**2-4._ark/9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+ &
          sqrt(3._ark)*y5*y6**2*y9**2+(y7**2-y8**2)*y6**3)*y2+((-4._ark/9._ark*sqrt(3._ark)*y8**2+ &
          5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2- &
          y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+(y7**2-y8**2)*y6**3)*y3+((-4._ark/ &
          9._ark*sqrt(3._ark)*y8**2+5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+ &
          (y7**2-y8**2)*y6**3)*y4
      dF(187) = ((y9+y7+y8)*y5**4+(2._ark*y8+2._ark*y9+2._ark*y7)*y6**2*y5**2+(y9+y7+ &
          y8)*y6**4)*y1+((-y7+y9-y8)*y5**4+(-2._ark*y7+2._ark*y9-2._ark*y8)*y6**2*y5**2+(-y7+y9- &
          y8)*y6**4)*y2+((-y8-y9+y7)*y5**4+(-2._ark*y8+2._ark*y7-2._ark*y9)*y6**2*y5**2+(-y8-y9+ &
          y7)*y6**4)*y3+((-y7-y9+y8)*y5**4+(-2._ark*y9+2._ark*y8-2._ark*y7)*y6**2*y5**2+(-y7-y9+ &
          y8)*y6**4)*y4
      s1 = ((sqrt(3._ark)*y9/24._ark+sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y5**4+(y7- &
          y8)*y6*y5**3+(sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/4._ark+sqrt(3._ark)*y7/ &
          2._ark)*y6**2*y5**2+3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y1+((-sqrt(3._ark)*y8/6._ark- &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/24._ark)*y5**4+(y8-y7)*y6*y5**3+(-sqrt(3._ark)*y8/ &
          2._ark-sqrt(3._ark)*y9/4._ark-sqrt(3._ark)*y7/2._ark)*y6**2*y5**2+3._ark/ &
          8._ark*sqrt(3._ark)*y6**4*y9)*y2
      dF(188) = s1+((-sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/24._ark+sqrt(3._ark)*y7/ &
          6._ark)*y5**4+(y7+y8)*y6*y5**3+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark+ &
          sqrt(3._ark)*y9/4._ark)*y6**2*y5**2-3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y3+ &
          ((sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/24._ark-sqrt(3._ark)*y7/6._ark)*y5**4+(-y8- &
          y7)*y6*y5**3+(sqrt(3._ark)*y9/4._ark-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6**2*y5**2-3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y4
      dF(189) = (-3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2)*y1+(-3._ark*y5*y6**4+y5**5- &
          2._ark*y5**3*y6**2)*y2+(-3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2)*y3+(-3._ark*y5*y6**4+ &
          y5**5-2._ark*y5**3*y6**2)*y4
      dF(190) = ((sqrt(3._ark)*y5**2*y9**2/2._ark-sqrt(3._ark)*y6**2*y9**2/6._ark)*y2+(- &
          y5*y6*y7**2+sqrt(3._ark)*y6**2*y7**2/3._ark)*y3+(y5*y6*y8**2+sqrt(3._ark)*y6**2*y8**2/ &
          3._ark)*y4)*y1+((y5*y6*y8**2+sqrt(3._ark)*y6**2*y8**2/3._ark)*y3+(-y5*y6*y7**2+ &
          sqrt(3._ark)*y6**2*y7**2/3._ark)*y4)*y2+(sqrt(3._ark)*y5**2*y9**2/2._ark- &
          sqrt(3._ark)*y6**2*y9**2/6._ark)*y4*y3
      dF(191) = ((sqrt(3._ark)*y5**2*y7*y8/2._ark-sqrt(3._ark)*y6**2*y7*y8/6._ark)*y2+ &
          (sqrt(3._ark)*y6**2*y8*y9/3._ark-y5*y6*y8*y9)*y3+(sqrt(3._ark)*y6**2*y7*y9/3._ark+ &
          y5*y6*y7*y9)*y4)*y1+((-sqrt(3._ark)*y6**2*y7*y9/3._ark-y5*y6*y7*y9)*y3+(- &
          sqrt(3._ark)*y6**2*y8*y9/3._ark+y5*y6*y8*y9)*y4)*y2+(sqrt(3._ark)*y6**2*y7*y8/6._ark- &
          sqrt(3._ark)*y5**2*y7*y8/2._ark)*y4*y3
      dF(192) = ((y8*y7*y5**2+y8*y7*y6**2)*y2+(y5**2*y8*y9+y9*y8*y6**2)*y3+ &
          (y5**2*y7*y9+y9*y7*y6**2)*y4)*y1+((-y9*y7*y6**2-y5**2*y7*y9)*y3+(-y5**2*y8*y9- &
          y9*y8*y6**2)*y4)*y2+(-y8*y7*y5**2-y8*y7*y6**2)*y4*y3
      dF(193) = ((2._ark/9._ark*sqrt(3._ark)*y6**4+2._ark*sqrt(3._ark)*y5**2*y6**2)*y2+(5._ark/ &
          3._ark*y5*y6**3+7._ark/18._ark*sqrt(3._ark)*y6**4-y5**3*y6+sqrt(3._ark)*y5**4/2._ark)*y3+ &
          (7._ark/18._ark*sqrt(3._ark)*y6**4+sqrt(3._ark)*y5**4/2._ark+y5**3*y6-5._ark/ &
          3._ark*y5*y6**3)*y4)*y1+((7._ark/18._ark*sqrt(3._ark)*y6**4+sqrt(3._ark)*y5**4/2._ark+ &
          y5**3*y6-5._ark/3._ark*y5*y6**3)*y3+(5._ark/3._ark*y5*y6**3+7._ark/18._ark*sqrt(3._ark)*y6**4- &
          y5**3*y6+sqrt(3._ark)*y5**4/2._ark)*y4)*y2+(2._ark/9._ark*sqrt(3._ark)*y6**4+ &
          2._ark*sqrt(3._ark)*y5**2*y6**2)*y4*y3
      dF(194) = (y4*y8**3+y2*y9**3+y3*y7**3)*y1**2+(y2**2*y9**3+y4**2*y8**3+ &
          y3**2*y7**3)*y1+(-y4*y7**3-y3*y8**3)*y2**2+(-y3**2*y8**3-y4**2*y7**3)*y2- &
          y3**2*y4*y9**3-y3*y4**2*y9**3
      s1 = (((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y7*y9/2._ark+y8*y9/ &
          2._ark)*y6)*y2+(-sqrt(3._ark)*y5*y7*y9/2._ark+(y8+y9/2._ark)*y7*y6)*y3+(- &
          sqrt(3._ark)*y5*y8*y9/2._ark+(-y8*y9/2._ark-y7*y8)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y8*y9/ &
          2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5+(y7*y9/2._ark-y8*y9/2._ark)*y6)*y2**2+ &
          (sqrt(3._ark)*y5*y7*y9/2._ark+(-y9/2._ark-y8)*y7*y6)*y3**2+(sqrt(3._ark)*y5*y8*y9/2._ark+ &
          (y7*y8+y8*y9/2._ark)*y6)*y4**2)*y1+((sqrt(3._ark)*y5*y8*y9/2._ark+(y8*y9/2._ark- &
          y7*y8)*y6)*y3+(sqrt(3._ark)*y5*y7*y9/2._ark+(y8-y9/2._ark)*y7*y6)*y4)*y2**2
      dF(195) = s1+((-sqrt(3._ark)*y5*y8*y9/2._ark+(y7*y8-y8*y9/2._ark)*y6)*y3**2+(- &
          sqrt(3._ark)*y5*y7*y9/2._ark+(y9/2._ark-y8)*y7*y6)*y4**2)*y2+((-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y5+(y8*y9/2._ark+y7*y9/2._ark)*y6)*y4*y3**2+((- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y8*y9/2._ark-y7*y9/ &
          2._ark)*y6)*y4**2*y3
      dF(196) = (((-y8-y7)*y5**2+(-y8-y7)*y6**2)*y2+((-y9-y8)*y5**2+(-y9- &
          y8)*y6**2)*y3+((-y9-y7)*y5**2+(-y9-y7)*y6**2)*y4)*y1**2+(((y7+y8)*y5**2+(y7+ &
          y8)*y6**2)*y2**2+((y9+y8)*y5**2+(y9+y8)*y6**2)*y3**2+((y7+y9)*y5**2+(y7+ &
          y9)*y6**2)*y4**2)*y1+(((y7-y9)*y5**2+(y7-y9)*y6**2)*y3+((y8-y9)*y5**2+(y8- &
          y9)*y6**2)*y4)*y2**2+(((y9-y7)*y5**2+(y9-y7)*y6**2)*y3**2+((-y8+y9)*y5**2+(-y8+ &
          y9)*y6**2)*y4**2)*y2+((y8-y7)*y5**2+(y8-y7)*y6**2)*y4*y3**2+((y7-y8)*y5**2+(y7- &
          y8)*y6**2)*y4**2*y3
      s1 = (((y8**2+y7**2)*y5+(sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y2+((y9**2- &
          2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y3+((-2._ark*y7**2+y9**2)*y5- &
          sqrt(3._ark)*y6*y9**2)*y4)*y1**2+(((y8**2+y7**2)*y5+(sqrt(3._ark)*y7**2- &
          sqrt(3._ark)*y8**2)*y6)*y2**2+((y9**2-2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y3**2+ &
          ((-2._ark*y7**2+y9**2)*y5-sqrt(3._ark)*y6*y9**2)*y4**2)*y1+(((-2._ark*y7**2+y9**2)*y5- &
          sqrt(3._ark)*y6*y9**2)*y3+((y9**2-2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y4)*y2**2
      dF(197) = s1+(((-2._ark*y7**2+y9**2)*y5-sqrt(3._ark)*y6*y9**2)*y3**2+((y9**2- &
          2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y4**2)*y2+((y8**2+y7**2)*y5+ &
          (sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y4*y3**2+((y8**2+y7**2)*y5+ &
          (sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y4**2*y3
      dF(198) = (-2._ark*y2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3+(y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y4)*y1**2+(-2._ark*y2**2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+ &
          y5*y7**2)*y3**2+(y5*y8**2+sqrt(3._ark)*y6*y8**2)*y4**2)*y1+((y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y3+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4)*y2**2+((y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y3**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4**2)*y2- &
          2._ark*y3*y4**2*y5*y9**2-2._ark*y3**2*y4*y5*y9**2
      s1 = (((y8*y9/2._ark+y7*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/ &
          2._ark)*y6)*y2+((y9/2._ark-y8)*y7*y5+sqrt(3._ark)*y6*y7*y9/2._ark)*y3+((y8*y9/2._ark- &
          y7*y8)*y5-sqrt(3._ark)*y6*y8*y9/2._ark)*y4)*y1**2+(((-y8*y9/2._ark-y7*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/2._ark)*y6)*y2**2+((y8-y9/2._ark)*y7*y5- &
          sqrt(3._ark)*y6*y7*y9/2._ark)*y3**2+((y7*y8-y8*y9/2._ark)*y5+sqrt(3._ark)*y6*y8*y9/ &
          2._ark)*y4**2)*y1+(((-y8*y9/2._ark-y7*y8)*y5+sqrt(3._ark)*y6*y8*y9/2._ark)*y3+((-y9/ &
          2._ark-y8)*y7*y5-sqrt(3._ark)*y6*y7*y9/2._ark)*y4)*y2**2
      dF(199) = s1+(((y7*y8+y8*y9/2._ark)*y5-sqrt(3._ark)*y6*y8*y9/2._ark)*y3**2+((y8+y9/ &
          2._ark)*y7*y5+sqrt(3._ark)*y6*y7*y9/2._ark)*y4**2)*y2+((-y7*y9/2._ark+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3**2+((y7*y9/2._ark-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**2*y3
      dF(200) = ((y9*y6**2+y9*y5**2)*y2+(y7*y6**2+y5**2*y7)*y3+(y8*y6**2+ &
          y5**2*y8)*y4)*y1**2+((y9*y6**2+y9*y5**2)*y2**2+(y7*y6**2+y5**2*y7)*y3**2+ &
          (y8*y6**2+y5**2*y8)*y4**2)*y1+((-y5**2*y8-y8*y6**2)*y3+(-y5**2*y7- &
          y7*y6**2)*y4)*y2**2+((-y5**2*y8-y8*y6**2)*y3**2+(-y5**2*y7-y7*y6**2)*y4**2)*y2+ &
          (-y9*y6**2-y9*y5**2)*y4*y3**2+(-y9*y6**2-y9*y5**2)*y4**2*y3
      dF(201) = ((y5**3-3._ark*y5*y6**2)*y2+(y5**3-3._ark*y5*y6**2)*y3+(y5**3- &
          3._ark*y5*y6**2)*y4)*y1**2+((y5**3-3._ark*y5*y6**2)*y2**2+(y5**3- &
          3._ark*y5*y6**2)*y3**2+(y5**3-3._ark*y5*y6**2)*y4**2)*y1+((y5**3-3._ark*y5*y6**2)*y3+ &
          (y5**3-3._ark*y5*y6**2)*y4)*y2**2+((y5**3-3._ark*y5*y6**2)*y3**2+(y5**3- &
          3._ark*y5*y6**2)*y4**2)*y2+(y5**3-3._ark*y5*y6**2)*y4*y3**2+(y5**3- &
          3._ark*y5*y6**2)*y4**2*y3
      dF(202) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1**3+((y8**2+ &
          y7**2)*y2**3+(y9**2+y8**2)*y3**3+(y7**2+y9**2)*y4**3)*y1+((y7**2+y9**2)*y3+ &
          (y9**2+y8**2)*y4)*y2**3+((y7**2+y9**2)*y3**3+(y9**2+y8**2)*y4**3)*y2+(y8**2+ &
          y7**2)*y4*y3**3+(y8**2+y7**2)*y4**3*y3
      dF(203) = ((-y8*y9-y7*y9)*y2+(-y9-y8)*y7*y3+(-y7*y8-y8*y9)*y4)*y1**3+((y8*y9+ &
          y7*y9)*y2**3+(y9+y8)*y7*y3**3+(y8*y9+y7*y8)*y4**3)*y1+((-y7*y8+y8*y9)*y3+(-y8+ &
          y9)*y7*y4)*y2**3+((-y8*y9+y7*y8)*y3**3+(y8-y9)*y7*y4**3)*y2+(-y8*y9+ &
          y7*y9)*y4*y3**3+(y8*y9-y7*y9)*y4**3*y3
      dF(204) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1**3+(y3**3*y7**2+y4**3*y8**2+ &
          y2**3*y9**2)*y1+(y4*y7**2+y3*y8**2)*y2**3+(y3**3*y8**2+y4**3*y7**2)*y2+ &
          y3**3*y4*y9**2+y3*y4**3*y9**2
      dF(205) = (y4*y7*y9+y3*y8*y9+y2*y7*y8)*y1**3+(y3**3*y8*y9+y4**3*y7*y9+ &
          y2**3*y7*y8)*y1+(-y4*y8*y9-y3*y7*y9)*y2**3+(-y3**3*y7*y9-y4**3*y8*y9)*y2- &
          y3**3*y4*y7*y8-y3*y4**3*y7*y8
      s1 = (((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/2._ark)*y6)*y2+(- &
          sqrt(3._ark)*y5*y8/2._ark+(y9+y8/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y7/2._ark+(-y9-y7/ &
          2._ark)*y6)*y4)*y1**3+(((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
          2._ark)*y6)*y2**3+(sqrt(3._ark)*y5*y8/2._ark+(-y8/2._ark-y9)*y6)*y3**3+ &
          (sqrt(3._ark)*y5*y7/2._ark+(y9+y7/2._ark)*y6)*y4**3)*y1+((sqrt(3._ark)*y5*y7/2._ark+(y7/ &
          2._ark-y9)*y6)*y3+(sqrt(3._ark)*y5*y8/2._ark+(-y8/2._ark+y9)*y6)*y4)*y2**3
      dF(206) = s1+((-sqrt(3._ark)*y5*y7/2._ark+(-y7/2._ark+y9)*y6)*y3**3+(- &
          sqrt(3._ark)*y5*y8/2._ark+(-y9+y8/2._ark)*y6)*y4**3)*y2+((sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark+y8/2._ark)*y6)*y4*y3**3+((-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/2._ark)*y6)*y4**3*y3
      s1 = (((y7/2._ark+y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y2+((- &
          y9+y8/2._ark)*y5+sqrt(3._ark)*y6*y8/2._ark)*y3+((y7/2._ark-y9)*y5-sqrt(3._ark)*y6*y7/ &
          2._ark)*y4)*y1**3+(((-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**3+((-y8/2._ark+y9)*y5-sqrt(3._ark)*y6*y8/2._ark)*y3**3+((-y7/2._ark+ &
          y9)*y5+sqrt(3._ark)*y6*y7/2._ark)*y4**3)*y1+(((-y9-y7/2._ark)*y5+sqrt(3._ark)*y6*y7/ &
          2._ark)*y3+((-y8/2._ark-y9)*y5-sqrt(3._ark)*y6*y8/2._ark)*y4)*y2**3
      dF(207) = s1+(((y9+y7/2._ark)*y5-sqrt(3._ark)*y6*y7/2._ark)*y3**3+((y9+y8/2._ark)*y5+ &
          sqrt(3._ark)*y6*y8/2._ark)*y4**3)*y2+((y7/2._ark-y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3**3+((y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4**3*y3
      dF(208) = ((-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y2+(-y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y1**3+((- &
          sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y2**3+(-y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3**3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4**3)*y1+((y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3+(-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y2**3+((y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3**3+(-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4**3)*y2+(-sqrt(3._ark)*y6**2/6._ark+ &
          sqrt(3._ark)*y5**2/2._ark)*y4*y3**3+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4**3*y3
      dF(209) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1**3+((y6**2+ &
          y5**2)*y2**3+(y6**2+y5**2)*y3**3+(y6**2+y5**2)*y4**3)*y1+((y6**2+y5**2)*y3+ &
          (y6**2+y5**2)*y4)*y2**3+((y6**2+y5**2)*y3**3+(y6**2+y5**2)*y4**3)*y2+(y6**2+ &
          y5**2)*y4*y3**3+(y6**2+y5**2)*y4**3*y3
      dF(210) = (y3*y7+y4*y8+y2*y9)*y1**4+(y2**4*y9+y4**4*y8+y3**4*y7)*y1+(-y3*y8- &
          y4*y7)*y2**4+(-y4**4*y7-y3**4*y8)*y2-y3**4*y4*y9-y3*y4**4*y9
      dF(211) = ((-y8-y7)*y2+(-y9-y8)*y3+(-y9-y7)*y4)*y1**4+((y7+y8)*y2**4+(y9+ &
          y8)*y3**4+(y7+y9)*y4**4)*y1+((y7-y9)*y3+(y8-y9)*y4)*y2**4+((y9-y7)*y3**4+(-y8+ &
          y9)*y4**4)*y2+(y8-y7)*y4*y3**4+(y7-y8)*y4**4*y3
      dF(212) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4)*y1**4+(-2._ark/3._ark*sqrt(3._ark)*y2**4*y5+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y3**4+(y6+sqrt(3._ark)*y5/3._ark)*y4**4)*y1+((y6+sqrt(3._ark)*y5/ &
          3._ark)*y3+(-y6+sqrt(3._ark)*y5/3._ark)*y4)*y2**4+((y6+sqrt(3._ark)*y5/3._ark)*y3**4+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**4)*y2-2._ark/3._ark*sqrt(3._ark)*y3*y4**4*y5-2._ark/ &
          3._ark*sqrt(3._ark)*y3**4*y4*y5
      dF(213) = ((y8**4+y7**4)*y2+(y8**4+y9**4)*y3+(y9**4+y7**4)*y4)*y1+((y9**4+ &
          y7**4)*y3+(y8**4+y9**4)*y4)*y2+(y8**4+y7**4)*y4*y3
      dF(214) = ((y7**3*y8+y7*y8**3)*y2+(y8**3*y9+y8*y9**3)*y3+(y7*y9**3+ &
          y7**3*y9)*y4)*y1+((-y7*y9**3-y7**3*y9)*y3+(-y8**3*y9-y8*y9**3)*y4)*y2+(- &
          y7**3*y8-y7*y8**3)*y4*y3
      dF(215) = (y4*y7**2*y9**2+y2*y7**2*y8**2+y3*y8**2*y9**2)*y1+(y3*y7**2*y9**2+ &
          y4*y8**2*y9**2)*y2+y3*y4*y7**2*y8**2
      dF(216) = (y3*y7**2*y8*y9+y4*y7*y8**2*y9+y2*y7*y8*y9**2)*y1+(-y4*y7**2*y8*y9- &
          y3*y7*y8**2*y9)*y2-y3*y4*y7*y8*y9**2
      dF(217) = (y4*y8**4+y3*y7**4+y2*y9**4)*y1+(y3*y8**4+y4*y7**4)*y2+y3*y4*y9**4
      dF(218) = (8._ark/3._ark*sqrt(3._ark)*y2*y5*y6**2*y9+(-sqrt(3._ark)*y5**3*y7+ &
          y5**2*y6*y7+5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7+y6**3*y7)*y3+(-y5**2*y6*y8+5._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2*y8-sqrt(3._ark)*y5**3*y8-y6**3*y8)*y4)*y1+((y6**3*y8+ &
          y5**2*y6*y8-5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8+sqrt(3._ark)*y5**3*y8)*y3+(- &
          y5**2*y6*y7-y6**3*y7+sqrt(3._ark)*y5**3*y7-5._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2*y7)*y4)*y2-8._ark/3._ark*sqrt(3._ark)*y3*y4*y5*y6**2*y9
      dF(219) = ((y5**2*y9**2+y6**2*y9**2)*y2+(y6**2*y7**2+y5**2*y7**2)*y3+ &
          (y6**2*y8**2+y5**2*y8**2)*y4)*y1+((y6**2*y8**2+y5**2*y8**2)*y3+(y6**2*y7**2+ &
          y5**2*y7**2)*y4)*y2+(y5**2*y9**2+y6**2*y9**2)*y4*y3
      dF(220) = ((4._ark/3._ark*y6**4+4._ark*y5**2*y6**2)*y2+(4._ark/3._ark*sqrt(3._ark)*y5*y6**3+ &
          y5**2*y6**2+5._ark/6._ark*y6**4+3._ark/2._ark*y5**4)*y3+(3._ark/2._ark*y5**4-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4)*y4)*y1+((3._ark/2._ark*y5**4- &
          4._ark/3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4)*y3+(4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4+3._ark/2._ark*y5**4)*y4)*y2+ &
          (4._ark/3._ark*y6**4+4._ark*y5**2*y6**2)*y4*y3
      dF(221) = ((y4*y7*y8*y9+y3*y7*y8*y9)*y2+y3*y4*y7*y8*y9)*y1+y2*y3*y4*y7*y8*y9
      dF(222) = (((((-y9/2._ark-y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark- &
          sqrt(3._ark)*y7*y9/2._ark)*y6)*y3+(((y9/2._ark-y8)*y7-y8*y9/2._ark)*y5+ &
          (sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4)*y2+(((y8+y9/2._ark)*y7+ &
          y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3)*y1+ &
          (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/ &
          2._ark)*y6)*y4*y3*y2
      dF(223) = ((((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y3+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4)*y2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4*y3)*y1+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4*y3*y2
      dF(224) = ((((y9+y7-y8)*y5**2+(y9+y7-y8)*y6**2)*y3+((y8+y9-y7)*y5**2+(y8+y9- &
          y7)*y6**2)*y4)*y2+((y7-y9+y8)*y5**2+(y7-y9+y8)*y6**2)*y4*y3)*y1+((-y7-y8- &
          y9)*y5**2+(-y7-y8-y9)*y6**2)*y4*y3*y2
      dF(225) = (((y5**3-3._ark*y5*y6**2)*y3+(y5**3-3._ark*y5*y6**2)*y4)*y2+(y5**3- &
          3._ark*y5*y6**2)*y4*y3)*y1+(y5**3-3._ark*y5*y6**2)*y4*y3*y2
      dF(226) = (((-sqrt(3._ark)*y6*y8-y5*y8)*y3+(sqrt(3._ark)*y6*y7-y5*y7)*y4)*y2+ &
          2._ark*y3*y4*y5*y9)*y1**2+(((-sqrt(3._ark)*y6*y7+y5*y7)*y3+(sqrt(3._ark)*y6*y8+ &
          y5*y8)*y4)*y2**2+(-2._ark*y3**2*y5*y9-2._ark*y4**2*y5*y9)*y2+(sqrt(3._ark)*y6*y8+ &
          y5*y8)*y4*y3**2+(-sqrt(3._ark)*y6*y7+y5*y7)*y4**2*y3)*y1+2._ark*y2**2*y3*y4*y5*y9+ &
          ((sqrt(3._ark)*y6*y7-y5*y7)*y4*y3**2+(-sqrt(3._ark)*y6*y8-y5*y8)*y4**2*y3)*y2
      dF(227) = (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3)*y1**2+ &
          (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2**2+((y6**2+y5**2)*y3**2+(y6**2+ &
          y5**2)*y4**2)*y2+(y6**2+y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3)*y1+(y6**2+ &
          y5**2)*y4*y3*y2**2+((y6**2+y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3)*y2
      dF(228) = ((sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/3._ark)*y2+(y6**3-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y3+(-y6**3-y5**2*y6-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2)*y4)*y1**2+((sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/ &
          3._ark)*y2**2+(y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y3**2+(-y6**3- &
          y5**2*y6-4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y4**2)*y1+((-y6**3-y5**2*y6-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2)*y3+(y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+ &
          y5**2*y6)*y4)*y2**2+((-y6**3-y5**2*y6-4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y3**2+ &
          (y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y4**2)*y2+(sqrt(3._ark)*y5**3- &
          sqrt(3._ark)*y5*y6**2/3._ark)*y4*y3**2+(sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/ &
          3._ark)*y4**2*y3
      dF(229) = ((-y8**2*y9+y7**2*y9)*y6*y2+((-sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/ &
          2._ark)*y7*y5+(y9**2/2._ark-y8**2/2._ark)*y7*y6)*y3+((-sqrt(3._ark)*y7**2*y8/2._ark+ &
          sqrt(3._ark)*y8*y9**2/2._ark)*y5+(-y8*y9**2/2._ark+y7**2*y8/2._ark)*y6)*y4)*y1+(((- &
          sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7**2*y8/2._ark)*y5+(-y7**2*y8/2._ark+y8*y9**2/ &
          2._ark)*y6)*y3+((sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y9**2/2._ark)*y7*y5+(y8**2/2._ark- &
          y9**2/2._ark)*y7*y6)*y4)*y2+(-y7**2*y9+y8**2*y9)*y6*y4*y3
      dF(230) = (y2*y5*y9**3+(sqrt(3._ark)*y6*y7**3/2._ark-y5*y7**3/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y8**3/2._ark-y5*y8**3/2._ark)*y4)*y1+((y5*y8**3/2._ark+ &
          sqrt(3._ark)*y6*y8**3/2._ark)*y3+(-sqrt(3._ark)*y6*y7**3/2._ark+y5*y7**3/2._ark)*y4)*y2- &
          y3*y4*y5*y9**3
      dF(231) = (y2*y5*y7*y8*y9+(sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/2._ark)*y3+(- &
          y5*y7*y8*y9/2._ark-sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y4)*y1+((-y5*y7*y8*y9/2._ark- &
          sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y3+(sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/ &
          2._ark)*y4)*y2+y3*y4*y5*y7*y8*y9
      dF(232) = ((y7**2*y9+y8**2*y9)*y5*y2+((-y8**2/2._ark-y9**2/2._ark)*y7*y5+ &
          (sqrt(3._ark)*y9**2/2._ark+sqrt(3._ark)*y8**2/2._ark)*y7*y6)*y3+((-y8*y9**2/2._ark- &
          y7**2*y8/2._ark)*y5+(-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y8/ &
          2._ark)*y6)*y4)*y1+(((y8*y9**2/2._ark+y7**2*y8/2._ark)*y5+(sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y8/2._ark)*y6)*y3+((y8**2/2._ark+y9**2/2._ark)*y7*y5+(- &
          sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y7*y6)*y4)*y2+(-y8**2*y9- &
          y7**2*y9)*y5*y4*y3
      dF(233) = (((-sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(y7**2- &
          y8**2)*y6*y5+(-sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y9**2+(-sqrt(3._ark)*y9**2/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y8**2)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y9**2+(- &
          sqrt(3._ark)*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2)*y4)*y1+((- &
          sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y9**2+(-sqrt(3._ark)*y9**2/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y7**2)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y9**2+(- &
          sqrt(3._ark)*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y8**2)*y6**2)*y4)*y2+((- &
          sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(y7**2-y8**2)*y6*y5+(- &
          sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y4*y3
      dF(234) = ((-3._ark*y5*y6**2*y9+y5**3*y9)*y2+(y5**3*y7-3._ark*y5*y6**2*y7)*y3+ &
          (y5**3*y8-3._ark*y5*y6**2*y8)*y4)*y1+((-y5**3*y8+3._ark*y5*y6**2*y8)*y3+(-y5**3*y7+ &
          3._ark*y5*y6**2*y7)*y4)*y2+(-y5**3*y9+3._ark*y5*y6**2*y9)*y4*y3
      dF(235) = ((-5._ark/3._ark*y6**4-6._ark*y5**2*y6**2+y5**4)*y2+(-8._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3-2._ark/3._ark*y6**4-2._ark*y5**4)*y3+(-2._ark*y5**4-2._ark/ &
          3._ark*y6**4+8._ark/3._ark*sqrt(3._ark)*y5*y6**3)*y4)*y1+((-2._ark*y5**4-2._ark/3._ark*y6**4+ &
          8._ark/3._ark*sqrt(3._ark)*y5*y6**3)*y3+(-8._ark/3._ark*sqrt(3._ark)*y5*y6**3-2._ark/ &
          3._ark*y6**4-2._ark*y5**4)*y4)*y2+(-5._ark/3._ark*y6**4-6._ark*y5**2*y6**2+y5**4)*y4*y3
      dF(236) = (((y7**3+y9**3-y8**3)*y3+(y8**3+y9**3-y7**3)*y4)*y2+(-y9**3+y8**3+ &
          y7**3)*y4*y3)*y1+(-y8**3-y7**3-y9**3)*y4*y3*y2
      dF(237) = ((((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+y8**2*y9)*y3+((y9+ &
          y8)*y7**2+(-y9**2-y8**2)*y7+y8**2*y9+y8*y9**2)*y4)*y2+((y8-y9)*y7**2+(y9**2+ &
          y8**2)*y7-y8**2*y9+y8*y9**2)*y4*y3)*y1+((-y9-y8)*y7**2+(-y9**2-y8**2)*y7- &
          y8*y9**2-y8**2*y9)*y4*y3*y2
      dF(238) = (((-sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y9/6._ark+ &
          sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9/2._ark+(- &
          y8-y7)*y6*y5+(sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/ &
          3._ark)*y6**2)*y4)*y2+(sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(-sqrt(3._ark)*y7/ &
          3._ark-sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y9/6._ark)*y6**2)*y4*y3)*y1+ &
          (sqrt(3._ark)*y5**2*y9/2._ark+(y8-y7)*y6*y5+(sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark- &
          sqrt(3._ark)*y9/6._ark)*y6**2)*y4*y3*y2
      dF(239) = (y8**2+y7**2+y9**2)*y4*y3*y2*y1
      dF(240) = (y6**2+y5**2)*y4*y3*y2*y1
      dF(241) = (y9**4+y8**4+y7**4)*y1**2+(y9**4+y8**4+y7**4)*y2**2+(y9**4+y8**4+ &
          y7**4)*y3**2+(y9**4+y8**4+y7**4)*y4**2
      dF(242) = (y7**2*y8*y9+(y8*y9**2+y8**2*y9)*y7)*y1**2+(-y7**2*y8*y9+(y8*y9**2- &
          y8**2*y9)*y7)*y2**2+(y7**2*y8*y9+(-y8**2*y9-y8*y9**2)*y7)*y3**2+(-y7**2*y8*y9+(- &
          y8*y9**2+y8**2*y9)*y7)*y4**2
      dF(243) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y1**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y2**2+((y9**2+y8**2)*y7**2+y8**2*y9**2)*y3**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y4**2
      dF(244) = ((y9+y8)*y7**3+(y9**3+y8**3)*y7+y8**3*y9+y8*y9**3)*y1**2+((y8- &
          y9)*y7**3+(y8**3-y9**3)*y7-y8*y9**3-y8**3*y9)*y2**2+((-y9-y8)*y7**3+(-y9**3- &
          y8**3)*y7+y8**3*y9+y8*y9**3)*y3**2+((-y8+y9)*y7**3+(y9**3-y8**3)*y7-y8*y9**3- &
          y8**3*y9)*y4**2
      s1 = ((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5+((y8+y9/2._ark)*y7**2+(-y9**2/2._ark-y8**2)*y7- &
          y8**2*y9/2._ark+y8*y9**2/2._ark)*y6)*y1**2+((-sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y5+ &
          ((y9/2._ark-y8)*y7**2+(y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/2._ark)*y6)*y2**2
      dF(245) = s1+((sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5+((-y9/2._ark-y8)*y7**2+(- &
          y9**2/2._ark-y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y6)*y3**2+ &
          ((sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y9/2._ark)*y5+((y8-y9/2._ark)*y7**2+(y9**2/2._ark+y8**2)*y7+y8**2*y9/ &
          2._ark+y8*y9**2/2._ark)*y6)*y4**2
      dF(246) = ((-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3)*y5+(-y8**3+y7**3)*y6)*y1**2+((sqrt(3._ark)*y7**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3+y8**3)*y6)*y2**2+((- &
          2._ark/3._ark*sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark)*y5+ &
          (y8**3+y7**3)*y6)*y3**2+((-2._ark/3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3-y8**3)*y6)*y4**2
      s1 = (((y8-y9/2._ark)*y7**2+(-y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/ &
          2._ark)*y5+(-sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y1**2+(((-y9/2._ark- &
          y8)*y7**2+(y9**2/2._ark-y8**2)*y7+y8*y9**2/2._ark-y8**2*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8**2*y9/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark)*y6)*y2**2
      dF(247) = s1+(((y9/2._ark-y8)*y7**2+(-y9**2/2._ark+y8**2)*y7+y8**2*y9/2._ark+y8*y9**2/ &
          2._ark)*y5+(sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/ &
          2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y3**2+(((y8+y9/2._ark)*y7**2+(y9**2/2._ark- &
          y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y5+(sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6)*y4**2
      dF(248) = ((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark-2._ark/3._ark*sqrt(3._ark)*y8)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y1**2+((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(y8*y9- &
          y7*y9)*y6*y5+((-2._ark/3._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/6._ark)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y2**2+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(- &
          y8*y9-y7*y9)*y6*y5+((sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7- &
          sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y3**2+((-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/ &
          2._ark)*y5**2+(y8*y9+y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y4**2
      dF(249) = (-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y1**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y2**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y3**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y4**2
      dF(250) = ((-y9/3._ark-y8/3._ark-y7/3._ark)*y5**3+(y9+y7+y8)*y6**2*y5)*y1**2+((y8/ &
          3._ark-y9/3._ark+y7/3._ark)*y5**3+(-y7+y9-y8)*y6**2*y5)*y2**2+((y9/3._ark+y8/3._ark-y7/ &
          3._ark)*y5**3+(-y8-y9+y7)*y6**2*y5)*y3**2+((y7/3._ark+y9/3._ark-y8/3._ark)*y5**3+(-y7- &
          y9+y8)*y6**2*y5)*y4**2
      dF(251) = ((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y1**2+((y8**2+ &
          y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y2**2+((y8**2+y7**2+y9**2)*y5**2+ &
          (y8**2+y7**2+y9**2)*y6**2)*y3**2+((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+ &
          y9**2)*y6**2)*y4**2
      dF(252) = ((-4._ark/9._ark*sqrt(3._ark)*y8-4._ark/9._ark*sqrt(3._ark)*y7+5._ark/ &
          9._ark*sqrt(3._ark)*y9)*y5**3+(y7-y8)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9+(y7- &
          y8)*y6**3)*y1**2+((4._ark/9._ark*sqrt(3._ark)*y7+5._ark/9._ark*sqrt(3._ark)*y9+4._ark/ &
          9._ark*sqrt(3._ark)*y8)*y5**3+(y8-y7)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9+(y8- &
          y7)*y6**3)*y2**2+((4._ark/9._ark*sqrt(3._ark)*y8-5._ark/9._ark*sqrt(3._ark)*y9-4._ark/ &
          9._ark*sqrt(3._ark)*y7)*y5**3+(y7+y8)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y9+(y7+ &
          y8)*y6**3)*y3**2+((-5._ark/9._ark*sqrt(3._ark)*y9+4._ark/9._ark*sqrt(3._ark)*y7-4._ark/ &
          9._ark*sqrt(3._ark)*y8)*y5**3+(-y8-y7)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y9+(-y8- &
          y7)*y6**3)*y4**2
      dF(253) = ((sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2+(- &
          sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3+(-sqrt(3._ark)*y6**2*y8/3._ark- &
          y5*y6*y8)*y4)*y1**2+((sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2**2+ &
          (-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3**2+(-sqrt(3._ark)*y6**2*y8/3._ark- &
          y5*y6*y8)*y4**2)*y1+((sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y3+(-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/3._ark)*y4)*y2**2+((sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y3**2+ &
          (-y5*y6*y7+sqrt(3._ark)*y6**2*y7/3._ark)*y4**2)*y2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
          sqrt(3._ark)*y5**2*y9/2._ark)*y4*y3**2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
          sqrt(3._ark)*y5**2*y9/2._ark)*y4**2*y3
      dF(254) = (2._ark/3._ark*sqrt(3._ark)*y2**2*y6**2+(y5*y6+sqrt(3._ark)*y5**2/2._ark+ &
          sqrt(3._ark)*y6**2/6._ark)*y3**2+(-y5*y6+sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4**2)*y1**2+((-y5*y6+sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y3**2+ &
          (y5*y6+sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y4**2)*y2**2+2._ark/ &
          3._ark*sqrt(3._ark)*y3**2*y4**2*y6**2
      dF(255) = (y2*y5*y7*y8+(sqrt(3._ark)*y6*y8*y9/2._ark-y5*y8*y9/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4)*y1**2+(y2**2*y5*y7*y8+ &
          (sqrt(3._ark)*y6*y8*y9/2._ark-y5*y8*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y7*y9/2._ark- &
          y5*y7*y9/2._ark)*y4**2)*y1+((sqrt(3._ark)*y6*y7*y9/2._ark+y5*y7*y9/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/2._ark)*y4)*y2**2+((sqrt(3._ark)*y6*y7*y9/2._ark+ &
          y5*y7*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/2._ark)*y4**2)*y2- &
          y3**2*y4*y5*y7*y8-y3*y4**2*y5*y7*y8
      dF(256) = ((y3*y7*y9+y4*y8*y9)*y2+y3*y4*y7*y8)*y1**2+((-y4*y7*y9- &
          y3*y8*y9)*y2**2+(-y3**2*y7*y8-y4**2*y7*y8)*y2-y3*y4**2*y8*y9-y3**2*y4*y7*y9)*y1+ &
          y2**2*y3*y4*y7*y8+(y3**2*y4*y8*y9+y3*y4**2*y7*y9)*y2
      dF(257) = (((y7**2+y9**2)*y3+(y9**2+y8**2)*y4)*y2+(y8**2+y7**2)*y4*y3)*y1**2+ &
          (((y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y2**2+((y8**2+y7**2)*y3**2+(y8**2+ &
          y7**2)*y4**2)*y2+(y7**2+y9**2)*y4*y3**2+(y9**2+y8**2)*y4**2*y3)*y1+(y8**2+ &
          y7**2)*y4*y3*y2**2+((y9**2+y8**2)*y4*y3**2+(y7**2+y9**2)*y4**2*y3)*y2
      s1 = (((sqrt(3._ark)*y5*y9/2._ark+(y7+y9/2._ark)*y6)*y3+(sqrt(3._ark)*y5*y9/2._ark+(-y9/ &
          2._ark-y8)*y6)*y4)*y2+((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/ &
          2._ark)*y6)*y4*y3)*y1**2+(((sqrt(3._ark)*y5*y9/2._ark+(y8-y9/2._ark)*y6)*y3+ &
          (sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/2._ark)*y6)*y4)*y2**2+(((-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark+y8/2._ark)*y6)*y3**2+((sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/2._ark)*y6)*y4**2)*y2+(-sqrt(3._ark)*y5*y9/2._ark+ &
          (-y9/2._ark+y7)*y6)*y4*y3**2+(-sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark-y8)*y6)*y4**2*y3)*y1
      dF(258) = s1+((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
          2._ark)*y6)*y4*y3*y2**2+((-sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y4*y3**2+(- &
          sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark-y7)*y6)*y4**2*y3)*y2
      s1 = ((((-y9/2._ark+y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3+((y8-y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4)*y2+((-y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3)*y1**2+((((-y9/2._ark-y8)*y5-sqrt(3._ark)*y6*y9/ &
          2._ark)*y3+((-y9/2._ark-y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**2+(((y8/2._ark-y7/ &
          2._ark)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6)*y3**2+((y7/2._ark-y8/ &
          2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y4**2)*y2+((y7+y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4*y3**2+((y8+y9/2._ark)*y5+sqrt(3._ark)*y6*y9/ &
          2._ark)*y4**2*y3)*y1
      dF(259) = s1+((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y4*y3*y2**2+(((y9/2._ark-y8)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4*y3**2+((-y7+ &
          y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y4**2*y3)*y2
      dF(260) = (((y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3+(-y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y4)*y2+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y4*y3)*y1**2+(((- &
          y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y2**2+((- &
          sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y3**2+(-sqrt(3._ark)*y6**2/6._ark+ &
          sqrt(3._ark)*y5**2/2._ark)*y4**2)*y2+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**2+(-y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y4**2*y3)*y1+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4*y3*y2**2+((-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**2+(y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y4**2*y3)*y2
      dF(261) = (y9+y7+y8)*y4*y3*y2*y1**2+((-y7+y9-y8)*y4*y3*y2**2+((-y8-y9+ &
          y7)*y4*y3**2+(-y7-y9+y8)*y4**2*y3)*y2)*y1
      dF(262) = (y4**2*y7*y9+y2**2*y7*y8+y3**2*y8*y9)*y1**2+(-y4**2*y8*y9- &
          y3**2*y7*y9)*y2**2-y3**2*y4**2*y7*y8
      dF(263) = (y2**2*y5*y9+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3**2+(-y5*y8/2._ark- &
          sqrt(3._ark)*y6*y8/2._ark)*y4**2)*y1**2+((y5*y8/2._ark+sqrt(3._ark)*y6*y8/2._ark)*y3**2+ &
          (y5*y7/2._ark-sqrt(3._ark)*y6*y7/2._ark)*y4**2)*y2**2-y3**2*y4**2*y5*y9
      dF(264) = ((y6**2+y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y1**2+ &
          ((y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2**2+(y6**2+y5**2)*y4**2*y3**2
      dF(265) = ((y3*y9+y4*y9)*y2**2+(y4**2*y8+y3**2*y7)*y2+y3*y4**2*y8+ &
          y3**2*y4*y7)*y1**2+((-y3**2*y8-y4**2*y7)*y2**2-y3**2*y4**2*y9)*y1+(-y3*y4**2*y7- &
          y3**2*y4*y8)*y2**2-y2*y3**2*y4**2*y9
      dF(266) = (((y7-y8)*y3+(y8-y7)*y4)*y2**2+((-y8+y9)*y3**2+(y9-y7)*y4**2)*y2+(y8- &
          y9)*y4*y3**2+(y7-y9)*y4**2*y3)*y1**2+(((y7+y9)*y3**2+(y9+y8)*y4**2)*y2**2+(y7+ &
          y8)*y4**2*y3**2)*y1+((-y9-y7)*y4*y3**2+(-y9-y8)*y4**2*y3)*y2**2+(-y8- &
          y7)*y4**2*y3**2*y2
      dF(267) = ((y4*y5+y3*y5)*y2**2+((-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y2+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2*y3)*y1**2+(((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/ &
          2._ark+sqrt(3._ark)*y6/2._ark)*y4**2)*y2**2+y3**2*y4**2*y5)*y1+((-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4*y3**2+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4**2*y3)*y2**2+ &
          y2*y3**2*y4**2*y5
      dF(268) = (y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2)*y1**2+((y3*y4**2+y3**2*y4)*y2**2+ &
          y2*y3**2*y4**2)*y1
      dF(269) = ((y3**2+y4**2)*y2**2+y3**2*y4**2)*y1**2+y2**2*y3**2*y4**2
      dF(270) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1**3+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2**3+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3**3+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**3
      dF(271) = (-sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(sqrt(3._ark)*y9/6._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y1**3+(-sqrt(3._ark)*y5**2*y9/2._ark+ &
          (y8-y7)*y6*y5+(sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/ &
          3._ark)*y6**2)*y2**3+(sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y8/3._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y9/6._ark)*y6**2)*y3**3+(sqrt(3._ark)*y5**2*y9/2._ark+(- &
          y8-y7)*y6*y5+(-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/ &
          3._ark)*y6**2)*y4**3
      dF(272) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1**3+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2**3+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3**3+((-y7-y9+y8)*y5**2+ &
          (-y7-y9+y8)*y6**2)*y4**3
      dF(273) = (y5**3-3._ark*y5*y6**2)*y1**3+(y5**3-3._ark*y5*y6**2)*y2**3+(y5**3- &
          3._ark*y5*y6**2)*y3**3+(y5**3-3._ark*y5*y6**2)*y4**3
      dF(274) = (y4**3+y3**3+y2**3)*y1**3+(y4**3+y3**3)*y2**3+y3**3*y4**3
      dF(275) = (2._ark/3._ark*sqrt(3._ark)*y2*y5*y9+(-sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3+(- &
          sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4)*y1**3+(2._ark/3._ark*sqrt(3._ark)*y2**3*y5*y9+(- &
          sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3**3+(-sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4**3)*y1+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4)*y2**3+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3**3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4**3)*y2- &
          2._ark/3._ark*sqrt(3._ark)*y3*y4**3*y5*y9-2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y5*y9
      dF(276) = ((y3*y8+y4*y7)*y2+y3*y4*y9)*y1**3+((-y3*y7-y4*y8)*y2**3+(-y3**3*y9- &
          y4**3*y9)*y2-y3**3*y4*y8-y3*y4**3*y7)*y1+y2**3*y3*y4*y9+(y3*y4**3*y8+ &
          y3**3*y4*y7)*y2
      dF(277) = (((y7+y9)*y3+(y9+y8)*y4)*y2+(y7+y8)*y4*y3)*y1**3+(((-y8+y9)*y3+(y9- &
          y7)*y4)*y2**3+((y7-y8)*y3**3+(y8-y7)*y4**3)*y2+(y7-y9)*y4*y3**3+(y8- &
          y9)*y4**3*y3)*y1+(-y8-y7)*y4*y3*y2**3+((-y9-y8)*y4*y3**3+(-y9-y7)*y4**3*y3)*y2
      dF(278) = (((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+sqrt(3._ark)*y6/ &
          2._ark)*y4)*y2+y3*y4*y5)*y1**3+(((-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4)*y2**3+(y3**3*y5+y4**3*y5)*y2+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4*y3**3+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4**3*y3)*y1+y2**3*y3*y4*y5+((-y5/ &
          2._ark+sqrt(3._ark)*y6/2._ark)*y4*y3**3+(-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y4**3*y3)*y2
      dF(279) = ((y7+y8)*y2**2+(y9+y8)*y3**2+(y7+y9)*y4**2)*y1**3+((-y8-y7)*y2**3+(- &
          y9-y8)*y3**3+(-y9-y7)*y4**3)*y1**2+((y9-y7)*y3**2+(-y8+y9)*y4**2)*y2**3+((y7- &
          y9)*y3**3+(y8-y9)*y4**3)*y2**2+(y7-y8)*y4**2*y3**3+(y8-y7)*y4**3*y3**2
      dF(280) = (y4**2*y8+y3**2*y7+y2**2*y9)*y1**3+(y3**3*y7+y2**3*y9+y4**3*y8)*y1**2+ &
          (-y3**2*y8-y4**2*y7)*y2**3+(-y4**3*y7-y3**3*y8)*y2**2-y3**3*y4**2*y9- &
          y3**2*y4**3*y9
      dF(281) = (2._ark/3._ark*sqrt(3._ark)*y2**2*y5+(-sqrt(3._ark)*y5/3._ark+y6)*y3**2+(- &
          sqrt(3._ark)*y5/3._ark-y6)*y4**2)*y1**3+(2._ark/3._ark*sqrt(3._ark)*y2**3*y5+(- &
          sqrt(3._ark)*y5/3._ark+y6)*y3**3+(-sqrt(3._ark)*y5/3._ark-y6)*y4**3)*y1**2+((- &
          sqrt(3._ark)*y5/3._ark-y6)*y3**2+(-sqrt(3._ark)*y5/3._ark+y6)*y4**2)*y2**3+((- &
          sqrt(3._ark)*y5/3._ark-y6)*y3**3+(-sqrt(3._ark)*y5/3._ark+y6)*y4**3)*y2**2+2._ark/ &
          3._ark*sqrt(3._ark)*y3**3*y4**2*y5+2._ark/3._ark*sqrt(3._ark)*y3**2*y4**3*y5
      dF(282) = ((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+y3**2*y4)*y1**3+((y4+ &
          y3)*y2**3+(y4**3+y3**3)*y2+y3**3*y4+y3*y4**3)*y1**2+((y3**2+y4**2)*y2**3+(y4**3+ &
          y3**3)*y2**2+y3**3*y4**2+y3**2*y4**3)*y1+(y3*y4**2+y3**2*y4)*y2**3+(y3**3*y4+ &
          y3*y4**3)*y2**2+(y3**3*y4**2+y3**2*y4**3)*y2
      dF(283) = y1**3*y2*y3*y4+(y2**3*y3*y4+(y3**3*y4+y3*y4**3)*y2)*y1
      dF(284) = (y8**2+y7**2+y9**2)*y1**4+(y8**2+y7**2+y9**2)*y2**4+(y8**2+y7**2+ &
          y9**2)*y3**4+(y8**2+y7**2+y9**2)*y4**4
      dF(285) = (y6**2+y5**2)*y1**4+(y6**2+y5**2)*y2**4+(y6**2+y5**2)*y3**4+(y6**2+ &
          y5**2)*y4**4
      dF(286) = (y3**2+y2**2+y4**2)*y1**4+(y4**4+y2**4+y3**4)*y1**2+(y3**2+ &
          y4**2)*y2**4+(y4**4+y3**4)*y2**2+y3**4*y4**2+y3**2*y4**4
      dF(287) = ((y4+y3)*y2+y3*y4)*y1**4+((y4+y3)*y2**4+(y4**4+y3**4)*y2+y3*y4**4+ &
          y3**4*y4)*y1+y2**4*y3*y4+(y3*y4**4+y3**4*y4)*y2
      dF(288) = (y9+y7+y8)*y1**5+(-y7+y9-y8)*y2**5+(-y8-y9+y7)*y3**5+(-y7-y9+y8)*y4**5
      dF(289) = y3**6+y4**6+y1**6+y2**6

     f = 0.0_ark

     do ipar = 1,289
       !
       !
       f = f + force(ipar)*dF(ipar)
!	print*,f,df(k)
       !
     enddo      

return
end 
