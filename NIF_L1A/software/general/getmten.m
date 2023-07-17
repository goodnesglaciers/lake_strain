function	[Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=getmten(strike,dip,rake,m0);
%	input strike, dip, rake in A+R convention, and seismic moment
%
% 	output is Mrr, Mtt, Mpp, Mrt, Mrp, Mtp

% convert to radians
      strike=strike*pi/180.0;
      dip=dip*pi/180.0;
      rake=rake*pi/180.0;
%
% original fortran from mtu
%      m(2)=-m0*((dsin(dip)*dcos(rake)*dsin(2d0*strike))+ (dsin(2d0*dip)*dsin(rake)*dsin(strike)*dsin(strike)))
%      m(6)=-m0*((dsin(dip)*dcos(rake)*dcos(2d0*strike))+ (0.5d0*dsin(2d0*dip)*dsin(rake)*dsin(2d0*strike)))
%      m(4)=-m0*((dcos(dip)*dcos(rake)*dcos(strike))+ (dcos(2d0*dip)*dsin(rake)*dsin(strike)))
%      m(3)=+m0*((dsin(dip)*dcos(rake)*dsin(2d0*strike))- (dsin(2d0*dip)*dsin(rake)*dcos(strike)*dcos(strike)))
%      m(5)=+m0*((dcos(dip)*dcos(rake)*dsin(strike))-(dcos(2d0*dip)*dsin(rake)*dcos(strike)))
%      m(1)=+m0*dsin(2d0*dip)*dsin(rake)

      m(2)=-m0*((sin(dip)*cos(rake)*sin(2*strike))+ (sin(2*dip)*sin(rake)*sin(strike)*sin(strike)));
      m(6)=-m0*((sin(dip)*cos(rake)*cos(2*strike))+ (0.5*sin(2*dip)*sin(rake)*sin(2*strike)));
      m(4)=-m0*((cos(dip)*cos(rake)*cos(strike))+ (cos(2*dip)*sin(rake)*sin(strike)));
      m(3)=+m0*((sin(dip)*cos(rake)*sin(2*strike))- (sin(2*dip)*sin(rake)*cos(strike)*cos(strike)));
      m(5)=+m0*((cos(dip)*cos(rake)*sin(strike))-(cos(2*dip)*sin(rake)*cos(strike)));
      m(1)=+m0*sin(2*dip)*sin(rake);

      Mrr=m(1); Mtt=m(2); Mpp=m(3); Mrt=m(4); Mrp=m(5); Mtp=m(6);

end
