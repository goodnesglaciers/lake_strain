        function [a,b,az] = ErrlpsGMT(cov)
%
% Compute two dimensional error ellipses given covariance matrix
% In form appropriate for plotting with GMT command psvelomeca
% Input:
%	cov	- covariance matrix(horizontal displacements only)
%		  use cov3d2cov2d for 3-dimensional covariance matrices
%		  [a,b,az] = ErrlpsGMT( cov3d2cov2d(cov)  );
% Output:
%	a 	- semi-major axis of error ellipse
%	b 	- semi-minor axis of error ellipse
%	b 	- orientation of semi-major axis 


	[a,b,az] = errlps( cov );
	a = a/2.448;
	b = b/2.448;
	az = 90*ones(size(az)) - az;




