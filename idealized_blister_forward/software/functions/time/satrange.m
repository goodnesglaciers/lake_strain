function [range,e1,e2,e3]=satrange(apcoords,t,To,xyz)
%SATRANGE    [range,e1,e2,e3]=satrange(apcoords,t,To,xyz)
%
%Returns, at the position specified in apcoords, the range to and
%unit vectors toward the satellites whose positions are given in xyz 
%at the times specified in t.  Input 'To' gives the satellite position times.

%Set constants, etc.

	c=299792458;
	omegaE=7292115.1467e-11;
    ns=length(xyz);
    ne=length(t);
    delta=inf;
    alpha=0;
    tau=0;
    range=0;

%Begin iteration

	while max(abs(delta(:)))>1.e-3
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         	
	%Get station coordinates in ECEF

        for i=ns:-1:1
            for j=ne:-1:1
                X(j,i)=lagrange([To,xyz{i}(:,1)],t(j)-tau,9);
                Y(j,i)=lagrange([To,xyz{i}(:,2)],t(j)-tau,9);
                Z(j,i)=lagrange([To,xyz{i}(:,3)],t(j)-tau,9);
            end
        end

	%Calculate range (alpha cancels the effect of earth rotation during travel time)

		Xdiff=(cos(alpha).*X+sin(alpha).*Y)-apcoords(1);
		Ydiff=(-sin(alpha).*X+cos(alpha).*Y)-apcoords(2);
		Zdiff=Z-apcoords(3);

   %Calculate corrections

		old_range=range;
		range=sqrt(Xdiff.^2+Ydiff.^2+Zdiff.^2);
		delta=range-old_range;
		tau=tau+delta/c;
        alpha=omegaE*tau;

	end

%Calculate unit vectors

	e1=Xdiff./range;
	e2=Ydiff./range;
	e3=Zdiff./range;
