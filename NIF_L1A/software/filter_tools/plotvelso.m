	function [outmatrix] = plotvelso(E, N, E_Sig, N_Sig, NE_corr, t, t0, tf, llh)

%%modified 1/99 using Dave Schaff's multiplication routines

%% To get predicted average velocities

	Nsites = size(E,1);

% find the indicies of the closest matching times

	[Y, ind2] = min(  abs(t -tf) );
	[Y, ind1] = min(  abs(t -t0) );

% Load up a vector in e,n, e, n form
	x = zeros(4*Nsites,1); covx = zeros(4*Nsites);
	for i=1:Nsites
		x(4*(i-1)+1) =  E(i,ind1);
		x(4*(i-1)+2) =  N(i,ind1);
		x(4*(i-1)+3) =  E(i,ind2);
		x(4*(i-1)+4) =  N(i,ind2);

		covx(4*(i-1)+1,4*(i-1)+1) =  E_Sig(i,ind1)^2;
		covx(4*(i-1)+2,4*(i-1)+2) =  N_Sig(i,ind1)^2;
		covx(4*(i-1)+1,4*(i-1)+2) =  NE_corr(i,ind1)*sqrt( E_Sig(i,ind1)^2 + ...
							N_Sig(i,ind1)^2 );
		covx(4*(i-1)+2,4*(i-1)+1) =  covx(4*(i-1)+1,4*(i-1)+2);

		covx(4*(i-1)+3,4*(i-1)+3) =  E_Sig(i,ind2)^2;
		covx(4*(i-1)+4,4*(i-1)+4) =  N_Sig(i,ind2)^2;
		covx(4*(i-1)+3,4*(i-1)+4) =  NE_corr(i,ind2)*sqrt( E_Sig(i,ind2)^2 + ...
							N_Sig(i,ind2)^2 );
		covx(4*(i-1)+4,4*(i-1)+3) =  covx(4*(i-1)+3,4*(i-1)+4);
	end


% Form matrix to difference

	A = zeros(2*Nsites, 4*Nsites);
	s = [ -1 0 1 0; 0 -1 0 1]/(t(ind2) - t(ind1));

	for i = 1:Nsites
		A(2*(i-1)+1:2*(i-1)+2,4*(i-1)+1:4*(i-1)+4) = s;
	end

% Difference the coordinates to form average velocity;

%%	V = A*x;
%%	covV = A*covx*A';
%%	V = mmwans(A,x);
	V = oneloop(A,x);
%%	temp = mmwans(covx, A');
%%	covV = mmwans(A, temp);
	temp = oneloop(covx, A');
	covV = oneloop(A, temp);

% Extract components
	for i = 1:Nsites
		Ve(i) = V(2*(i-1)+1);
		Vn(i) = V(2*(i-1)+2);

		sigE(i) = sqrt(  covV(2*(i-1)+1,2*(i-1)+1)  );
		sigN(i) = sqrt(  covV(2*(i-1)+2,2*(i-1)+2)  );
		Corr(i) = covV(2*(i-1)+1,2*(i-1)+2)/ sqrt( covV(2*(i-1)+1,2*(i-1)+1) + ...
							covV(2*(i-1)+2,2*(i-1)+2)  );
	end

%% Need to plot in GMT format


% Plot predicted vectors
        outmatrix = [llh(2,1:Nsites)', llh(1,1:Nsites)', 1000 * Ve', ...
                 1000 * Vn', 1000*sigE', 1000*sigN', Corr'];
 

%% THis is bad, but done for expediancy
	for i =1:Nsites
		for j=3:7
			outmatrix(i,j) = subnan( outmatrix(i,j) );
		end
	end
%% Done
