function varargout=okada85(e,n,depth,strike,dip,L,W,rake,slip,U3,nu)
%OKADA85 Surface deformation due to a finite rectangular source.
%	[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(E,N,...
%		DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,NU)
%	computes Okada's 1985 solution for displacements, tilts and strains in a
%	geographic referential (East,North,Up). Coordinates (0,0) and depth 
%	correspond to the fault centroid.
%
%		E, N    : Coordinates of observation points (relative to fault centroid)
%		DEPTH   : Depth of the fault centroid (DEP > 0)
%		STRIKE  : Strike-angle from North (in degrees)
%		DIP     : Dip-angle (in degrees)
%		LENGTH  : Fault length in strike direction (LEN > 0)
%		WIDTH   : Fault width in dip direction (WIDTH > 0)
%		RAKE    : Slip-angle direction on fault plane (in degrees)
%		SLIP    : Dislocation in rake direction
%		OPEN    : Dislocation in tensile component
%		NU      : Poisson's ratio (optional, default is 0.25)
%
%	returns the following variables (same size as E and N):
%	uN,uE,uZ        : Displacements (Unit of SLIP and OPEN)
%	uZE,uZN         : Tilts (in radian)
%	uNN,uNE,uEN,uEE : Strains (Unit of SLIP and OPEN)/(Unit of N,E,..,WIDTH)
%
%	It is also possible to produce partial outputs, with following syntax:
%		[uE,uN,uZ] = OKADA85(...) for displacements only;
%		[uE,uN,uZ,uZE,uZN] = OKADA85(...) for displacements and tilts;
%		[uE,uN,uZ,uNN,uNE,uEN,uEE] = OKADA85(...) for displacements and strains;
%		[uZE,uZN] = OKADA85(...) for tilts only;
%		[uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...) for tilts and strains;
%		[uNN,uNE,uEN,uEE] = OKADA85(...) for strains only.
%
%	OKADA85(...) without output argument produces a 3-D figure with fault
%	geometry and dislocation. Example:
%		okada85(0,0,2,60,45,5,3,-45,1,1)
%		for a 5x3 fault, 60° strike, 45° dip, and unit dislocation in all
%		directions (reverse, senestral and open).
%
%	Equations are all vectorized excepted for argument DIP which has to be
%	a scalar; all other arguments can be scalar or matrix of the same size.
%
%	Formulas and notations similar to Okada [1985]. Convention for fault
%	parameters from Aki & Richards [1980] definition:
%		STRIKE = fault trace direction (0 to 360° relative to North), 
%			defined so that the fault dips to the right side of the trace.
%		DIP = angle of the fault (0 to 90°, relative to horizontal).
%		RAKE = direction the hanging wall moves during rupture, measured 
%			relative to the fault strike (-180 to 180°).
%	Examples:
%		DIP=90 & RAKE=0   : left lateral (senestral) strike slip
%		DIP=90 & RAKE=180 : right lateral (dextral) strike slip
%		DIP=45 & RAKE=90  : reverse fault
%		DIP=45 & RAKE=-90 : normal fault 
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%	   Institut de Physique du Globe de Paris
%	Created: 1997
%	Updated: 2010-11-29
%
%	References:
%	   Aki K., and P. G. Richards, Quantitative seismology, Freemann & Co,
%	      New York, 1980.
%	   Okada Y., Surface deformation due to shear and tensile faults in a
%	      half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.
%
%	Acknowledgments: Dmitry Nicolsky, University of Alaska

%	Development history:
%		[2010-11-29]: change coordinates and depth to fault centroid 
%		(instead of middle top edge).
%		[2010-09-24]: bugs correction in the syntax of I1, K2 and uyy_tf
%		functions, affecting some components. Detected by Dmitry Nicolsky.
%
%	Copyright (c) 1997-2010, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

error(nargchk(10,11,nargin))

if nargin < 11
	nu = 0.25;
end

if ~isnumeric(e) | ~isnumeric(n) | ~isnumeric(depth) ...
		| ~isnumeric(strike) | ~isnumeric(dip) | ~isnumeric(L) | ~isnumeric(W) ...
		| ~isnumeric(rake) | ~isnumeric(slip) | ~isnumeric(U3) | ~isnumeric(nu)
		error('All input arguments must be numeric.')
end

if numel(dip) ~= 1
	error('DIP argument must be scalar.')
end

% From degrees to radians
strike = strike*pi/180;
delta = dip*pi/180;
rake = rake*pi/180;

% Dislocation in the fault plane system
U1 = cos(rake).*slip;
U2 = sin(rake).*slip;

% Converts fault coordinates into Okada's reference system
depth = depth + sin(delta).*W/2;
e = e + cos(strike)*cos(delta)*W/2;
n = n + sin(strike)*cos(delta)*W/2;
x = cos(strike).*n + sin(strike).*e + L/2;
y = sin(strike).*n - cos(strike).*e + cos(delta).*W;

% Variable substitution (independent from xi and eta)
p = y.*cos(delta) + depth.*sin(delta);
q = y.*sin(delta) - depth.*cos(delta);

% Displacements
if any(nargout==[3, 5, 7, 9])
	ux = -U1/(2*pi) * chinnery(@ux_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		- U2/(2*pi) * chinnery(@ux_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		+ U3/(2*pi) * chinnery(@ux_tf,x,p,L,W,q,delta,nu); ... % tensile fault

	uy = -U1/(2*pi) * chinnery(@uy_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		- U2/(2*pi) * chinnery(@uy_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		+ U3/(2*pi) * chinnery(@uy_tf,x,p,L,W,q,delta,nu); ... % tensile fault

	uz = -U1/(2*pi) * chinnery(@uz_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		- U2/(2*pi) * chinnery(@uz_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		+ U3/(2*pi) * chinnery(@uz_tf,x,p,L,W,q,delta,nu); ... % tensile fault

	% Rotation from Okada's axes to geographic
	ue = sin(strike).*ux - cos(strike).*uy;
	un = cos(strike).*ux + sin(strike).*uy;
end

% Tilt
if any(nargout==[2, 5, 6, 9])
	uzx = -U1/(2*pi) * chinnery(@uzx_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		 - U2/(2*pi) * chinnery(@uzx_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		 + U3/(2*pi) * chinnery(@uzx_tf,x,p,L,W,q,delta,nu); ... % tensile fault

	uzy = -U1/(2*pi) * chinnery(@uzy_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		 - U2/(2*pi) * chinnery(@uzy_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		 + U3/(2*pi) * chinnery(@uzy_tf,x,p,L,W,q,delta,nu); ... % tensile fault

	% Rotation from Okada's axes to geographic
	uze = sin(strike).*uzx - cos(strike).*uzy;
	uzn = cos(strike).*uzx + sin(strike).*uzy;
end

% Strain
if any(nargout==[4, 6, 7, 9])
	uxx = -U1/(2*pi) * chinnery(@uxx_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		 - U2/(2*pi) * chinnery(@uxx_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		 + U3/(2*pi) * chinnery(@uxx_tf,x,p,L,W,q,delta,nu); ... % tensile fault
	uxy = -U1/(2*pi) * chinnery(@uxy_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		 - U2/(2*pi) * chinnery(@uxy_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		 + U3/(2*pi) * chinnery(@uxy_tf,x,p,L,W,q,delta,nu); ... % tensile fault
	uyx = -U1/(2*pi) * chinnery(@uyx_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		 - U2/(2*pi) * chinnery(@uyx_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		 + U3/(2*pi) * chinnery(@uyx_tf,x,p,L,W,q,delta,nu); ... % tensile fault
	uyy = -U1/(2*pi) * chinnery(@uyy_ss,x,p,L,W,q,delta,nu) ... % strike-slip
		 - U2/(2*pi) * chinnery(@uyy_ds,x,p,L,W,q,delta,nu) ... % delta-slip
		 + U3/(2*pi) * chinnery(@uyy_tf,x,p,L,W,q,delta,nu); ... % tensile fault

	% Rotation from Okada's axes to geographic
	unn = cos(strike)^2*uxx + sin(2*strike)*(uxy + uyx)/2 + sin(strike)^2*uyy;
	une = sin(2*strike)*(uxx - uyy)/2 + sin(strike)^2*uyx - cos(strike)^2*uxy;
	uen = sin(2*strike)*(uxx - uyy)/2 - cos(strike)^2*uyx + sin(strike)^2*uxy;
	uee = sin(strike)^2*uxx - sin(2*strike)*(uyx + uxy)/2 + cos(strike)^2*uyy;
end

% Assigns output arguments
switch nargout
	case 2
		varargout = {uze, uzn};
	case 3
		varargout = {ue, un, uz};
	case 4
		varargout = {unn, une, uen, uee};
	case 5
		varargout = {ue, un, uz, uze, uzn};
	case 6
		varargout = {uze, ezn, unn, une, uen, uee};
	case 7
		varargout = {ue, un, uz, unn, une, uen, uee};
	case 9
		varargout = {ue, un, uz, uze, uzn, unn, une, uen, uee};
	case 0
		% no output argument: plots geometry of the fault and dislocation
		figure
		plot(e,n,'.r','MarkerSize',.1)
		alpha = pi/2 - strike;
		x_fault = L/2*cos(alpha)*[-1,1,1,-1] + sin(alpha)*cos(delta)*W/2*[0,0,1,1];
		y_fault = L/2*sin(alpha)*[-1,1,1,-1] + cos(alpha)*cos(delta)*W/2*[0,0,-1,-1];
		z_fault = -depth + sin(delta)*W*[1,1,0,0];
		ddx = U1*cos(alpha) - U2*sin(alpha)*cos(delta) + U3*sin(alpha)*sin(delta);
		ddy = U1*sin(alpha) + U2*cos(alpha)*cos(delta) - U3*cos(alpha)*sin(delta);
		ddz = U2*sin(delta) + U3*cos(delta);
		patch(x_fault,y_fault,z_fault,.3*[1,1,1],'EdgeColor','k','LineWidth',2)
		patch(x_fault+ddx/2,y_fault+ddy/2,z_fault+ddz/2,.6*[1,1,1], ...
			'EdgeColor','k','LineWidth',1,'FaceAlpha',.5)
		patch(x_fault-ddx/2,y_fault-ddy/2,z_fault-ddz/2,.6*[1,1,1], ...
			'EdgeColor','k','LineWidth',1,'FaceAlpha',.5)
		xlabel('East'), ylabel('North'), zlabel('Vertical')
		view(3), grid on, axis equal, rotate3d
	otherwise
		disp('Unvalid number of output arguments.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes for I... and K... subfunctions:
%
%	1. original formulas use Lame's constant as mu/(mu+lambda) which
%	   depends only on the Poisson's ratio = 1 - 2*nu
%	2. tests for cos(delta) == 0 are made with "cos(delta) > eps" 
%	   because cos(90*pi/180) is not zero but = 6.1232e-17 (!)


% =================================================================
% Chinnery's notation [equation (24) p. 1143]

% -----------------------------------------------------------------
function u=chinnery(f,x,p,L,W,q,delta,nu)
u = feval(f,x,p,q,delta,nu) ...
	- feval(f,x,p-W,q,delta,nu) ...
	- feval(f,x-L,p,q,delta,nu) ...
	+ feval(f,x-L,p-W,q,delta,nu);


% =================================================================
% Displacement subfunctions

% strike-slip displacement subfunctions [equation (25) p. 1144]

% -----------------------------------------------------------------
function u=ux_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./(R.*(R + eta)) ...
	+ atan(xi.*eta./(q.*R)) ...
	+ I1(xi,eta,q,delta,nu,R).*sin(delta);

% -----------------------------------------------------------------
function u=uy_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta*cos(delta) + q*sin(delta)).*q./(R.*(R + eta)) ...
	+ q.*cos(delta)./(R + eta) ...
	+ I2(eta,q,delta,nu,R).*sin(delta);

% -----------------------------------------------------------------
function u=uz_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
u = (eta*sin(delta) - q*cos(delta)).*q./(R.*(R + eta)) ...
	+ q.*sin(delta)./(R + eta) ...
	+ I4(db,eta,q,delta,nu,R).*sin(delta);

% dip-slip displacement subfunctions [equation (26) p. 1144]
% -----------------------------------------------------------------
function u=ux_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q./R ...
	- I3(eta,q,delta,nu,R).*sin(delta).*cos(delta);

% -----------------------------------------------------------------
function u=uy_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta*cos(delta) + q*sin(delta)).*q./(R.*(R + xi)) ...
	+ cos(delta).*atan(xi.*eta./(q.*R)) ...
	- I1(xi,eta,q,delta,nu,R).*sin(delta).*cos(delta);

% -----------------------------------------------------------------
function u=uz_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
u = db.*q./(R.*(R + xi)) ...
	+ sin(delta).*atan(xi.*eta./(q.*R)) ...
	- I5(xi,eta,q,delta,nu,R,db).*sin(delta).*cos(delta);

% tensile fault displacement subfunctions [equation (27) p. 1144]
% -----------------------------------------------------------------
function u=ux_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q.^2 ./(R.*(R + eta)) ...
	- I3(eta,q,delta,nu,R).*sin(delta).^2;

% -----------------------------------------------------------------
function u=uy_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = -(eta*sin(delta) - q*cos(delta)).*q./(R.*(R + xi)) ...
	- sin(delta).*(xi.*q./(R.*(R + eta)) ...
	- atan(xi.*eta./(q.*R))) ...
	- I1(xi,eta,q,delta,nu,R).*sin(delta).^2;

% -----------------------------------------------------------------
function u=uz_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
u = (eta*cos(delta) + q*sin(delta)).*q./(R.*(R + xi)) ...
	+ cos(delta).*(xi.*q./(R.*(R + eta)) ...
	- atan(xi.*eta./(q.*R))) ...
	- I5(xi,eta,q,delta,nu,R,db).*sin(delta).^2;


% I... displacement subfunctions [equations (28) (29) p. 1144-1145]
% -----------------------------------------------------------------
function I=I1(xi,eta,q,delta,nu,R)
db = eta*sin(delta) - q*cos(delta);
if cos(delta) > eps
	I = (1 - 2*nu) * (-xi./(cos(delta).*(R + db))) - sin(delta)./cos(delta) .*I5(xi,eta,q,delta,nu,R,db);
else
	I = -(1 - 2*nu)/2 * xi.*q./(R + db).^2;
end

% -----------------------------------------------------------------
function I=I2(eta,q,delta,nu,R)
I = (1 - 2*nu) * (-log(R + eta)) - I3(eta,q,delta,nu,R);

% -----------------------------------------------------------------
function I=I3(eta,q,delta,nu,R)
yb = eta*cos(delta) + q*sin(delta);
db = eta*sin(delta) - q*cos(delta);
if cos(delta) > eps
	I = (1 - 2*nu) * (yb./(cos(delta)*(R + db)) - log(R + eta)) ...
		+ sin(delta)./cos(delta) * I4(db,eta,q,delta,nu,R);
else
	I = (1 - 2*nu)/2 * (eta./(R + db) + yb.*q./(R + db).^2 - log(R + eta));
end

% -----------------------------------------------------------------
function I=I4(db,eta,q,delta,nu,R)
if cos(delta) > eps
	I = (1 - 2*nu) * 1/cos(delta) * (log(R + db) - sin(delta)*log(R + eta));
else
	I = -(1 - 2*nu) * q./(R + db);
end

% -----------------------------------------------------------------
function I=I5(xi,eta,q,delta,nu,R,db)
X = sqrt(xi.^2 + q.^2);
if cos(delta) > eps
	I = (1 - 2*nu) * 2./cos(delta) ...
		.* atan((eta.*(X + q.*cos(delta)) + X.*(R + X).*sin(delta))./(xi.*(R + X).*cos(delta)));
else
	I = -(1 - 2*nu) * xi.*sin(delta)./(R + db);
end


% =================================================================
% Tilt subfunctions

% strike-slip tilt subfunctions [equation (37) p. 1147]

% -----------------------------------------------------------------
function u=uzx_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = -xi.*q.^2.*A(eta,R).*cos(delta) ...
	+ ((xi.*q)./R.^3 - K1(xi,eta,q,delta,nu,R))*sin(delta);

% -----------------------------------------------------------------
function u=uzy_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
yb = eta*cos(delta) + q*sin(delta);
u = (db.*q./R.^3).*cos(delta) ...
	+ (xi.^2.*q.*A(eta,R)*cos(delta) - sin(delta)./R + yb.*q./R.^3 - K2(xi,eta,q,delta,nu,R))*sin(delta);

% dip-slip tilt subfunctions [equation (38) p. 1147]

% -----------------------------------------------------------------
function u=uzx_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
u = db.*q./R.^3 ...
	+ q*sin(delta)./(R.*(R + eta)) ...
	+ K3(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta);

% -----------------------------------------------------------------
function u=uzy_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
yb = eta*cos(delta) + q*sin(delta);
u = yb.*db.*q.*A(xi,R) ...
	- (2*db./(R.*(R + xi)) + xi*sin(delta)./(R.*(R + eta)))*sin(delta) ...
	+ K1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta);

% tensile fault tilt subfunctions [equation (39) p. 1147]

% -----------------------------------------------------------------
function u=uzx_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q.^2./R.^3*sin(delta) ...
	- q.^3.*A(eta,R)*cos(delta) ...
	+ K3(xi,eta,q,delta,nu,R)*sin(delta)^2;

% -----------------------------------------------------------------
function u=uzy_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
yb = eta*cos(delta) + q*sin(delta);
u = (yb*sin(delta) + db*cos(delta)).*q.^2.*A(xi,R) ...
	+ xi.*q.^2.*A(eta,R)*sin(delta)*cos(delta) ...
	- (2*q./(R.*(R + xi)) - K1(xi,eta,q,delta,nu,R))*sin(delta)^2;

% -----------------------------------------------------------------
function a=A(x,R)
a = (2*R + x)./(R.^3.*(R + x).^2);

% K... tilt subfunctions [equations (40) (41) p. 1148]
% -----------------------------------------------------------------
function K=K1(xi,eta,q,delta,nu,R)
db = eta*sin(delta) - q*cos(delta);
if cos(delta) > eps
	K = (1 - 2*nu) * xi/cos(delta) .* (1./(R.*(R + db)) - sin(delta)./(R.*(R + eta)));
else
	K = (1 - 2*nu) * xi.*q./(R + db).^2;
end

% -----------------------------------------------------------------
function K=K2(xi,eta,q,delta,nu,R)
K = (1 - 2*nu) * (-sin(delta)./R + q*cos(delta)./(R.*(R + eta))) - K3(xi,eta,q,delta,nu,R);

% -----------------------------------------------------------------
function K=K3(xi,eta,q,delta,nu,R)
db = eta*sin(delta) - q*cos(delta);
yb = eta*cos(delta) + q*sin(delta);
if cos(delta) > eps
	K = (1 - 2*nu) * 1/cos(delta) .* (q./(R.*(R + eta)) - yb./(R.*(R + db)));
else
	K = (1 - 2*nu) * sin(delta)./(R + db) .* (xi.^2./(R.*(R + db)) - 1);
end


% =================================================================
% Strain subfunctions

% strike-slip strain subfunctions [equation (31) p. 1145]

% -----------------------------------------------------------------
function u=uxx_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.^2.*q.*A(eta,R) ...
	- J1(xi,eta,q,delta,nu,R)*sin(delta);

% -----------------------------------------------------------------
function u=uxy_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
u = xi.^3.*db./(R.^3.*(eta.^2 + q.^2)) ...
	- (xi.^3.*A(eta,R) + J2(xi,eta,q,delta,nu,R))*sin(delta);

% -----------------------------------------------------------------
function u=uyx_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./R.^3*cos(delta) ...
	+ (xi.*q.^2.*A(eta,R) - J2(xi,eta,q,delta,nu,R))*sin(delta);

% -----------------------------------------------------------------
function u=uyy_ss(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta*cos(delta) + q*sin(delta);
u = yb.*q./R.^3*cos(delta) ...
	+ (q.^3.*A(eta,R)*sin(delta) - 2*q*sin(delta)./(R.*(R + eta)) ...
		- (xi.^2 + eta.^2)./R.^3*cos(delta) - J4(xi,eta,q,delta,nu,R))*sin(delta);
	
% dip-slip strain subfunctions [equation (32) p. 1146]

% -----------------------------------------------------------------
function u=uxx_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./R.^3 ...
	+ J3(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta);

% -----------------------------------------------------------------
function u=uxy_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta*cos(delta) + q*sin(delta);
u = yb.*q./R.^3 ...
	- sin(delta)./R ...
	+ J1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta);

% -----------------------------------------------------------------
function u=uyx_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta*cos(delta) + q*sin(delta);
u = yb.*q./R.^3 ...
	+ q*cos(delta)./(R.*(R + eta)) ...
	+ J1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta);

% -----------------------------------------------------------------
function u=uyy_ds(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta*cos(delta) + q*sin(delta);
u = yb.^2.*q.*A(xi,R) ...
	- (2*yb./(R.*(R + xi)) + xi*cos(delta)./(R.*(R + eta)))*sin(delta) ...
	+ J2(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta);

% tensile fault strain subfunctions [equation (33) p. 1146]

% -----------------------------------------------------------------
function u=uxx_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q.^2.*A(eta,R) ...
	+ J3(xi,eta,q,delta,nu,R)*sin(delta)^2;

% -----------------------------------------------------------------
function u=uxy_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
u = -db.*q./R.^3 ...
	- xi.^2.*q.*A(eta,R)*sin(delta) ...
	+ J1(xi,eta,q,delta,nu,R)*sin(delta)^2;

% -----------------------------------------------------------------
function u=uyx_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q.^2./R.^3*cos(delta) ...
	+ q.^3.*A(eta,R)*sin(delta) ...
	+ J1(xi,eta,q,delta,nu,R)*sin(delta)^2;

% -----------------------------------------------------------------
function u=uyy_tf(xi,eta,q,delta,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(delta) - q*cos(delta);
yb = eta*cos(delta) + q*sin(delta);
u = (yb*cos(delta) - db*sin(delta)).*q.^2.*A(xi,R) ...
	- q*sin(2*delta)./(R.*(R + xi)) ...
	- (xi.*q.^2.*A(eta,R) - J2(xi,eta,q,delta,nu,R))*sin(delta)^2;


% J... tensile fault subfunctions [equations (34) (35) p. 1146-1147]
% -----------------------------------------------------------------
function J=J1(xi,eta,q,delta,nu,R)
db = eta*sin(delta) - q*cos(delta);
if cos(delta) > eps
	J = (1 - 2*nu) * 1/cos(delta) * (xi.^2./(R.*(R + db).^2) - 1./(R + db)) ...
		- sin(delta)/cos(delta)*K3(xi,eta,q,delta,nu,R);
else
	J = (1 - 2*nu)/2 * q./(R + db).^2 .* (2*xi.^2./(R.*(R + db)) - 1);
end

% -----------------------------------------------------------------
function J=J2(xi,eta,q,delta,nu,R)
db = eta*sin(delta) - q*cos(delta);
yb = eta*cos(delta) + q*sin(delta);
if cos(delta) > eps
	J = (1 - 2*nu) * 1/cos(delta) * xi.*yb./(R.*(R + db).^2) ...
		- sin(delta)/cos(delta)*K1(xi,eta,q,delta,nu,R);
else
	J = (1 - 2*nu)/2 * xi*sin(delta)./(R + db).^2 .* (2*q.^2./(R.*(R + db)) - 1);
end

% -----------------------------------------------------------------
function J=J3(xi,eta,q,delta,nu,R)
J = (1 - 2*nu) * -xi./(R.*(R + eta)) ...
	- J2(xi,eta,q,delta,nu,R);

% -----------------------------------------------------------------
function J=J4(xi,eta,q,delta,nu,R)
J = (1 - 2*nu) * (-cos(delta)./R - q*sin(delta)./(R.*(R + eta))) ...
	- J1(xi,eta,q,delta,nu,R);
