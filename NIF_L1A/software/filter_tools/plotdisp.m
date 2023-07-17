        function [outmatrix] = plotdisp(Ep, Np, t, t0, tf, llh)
 
%% To get predicted displacements
 
        [Nsites, oner] = size(Ep);
 
% find the indicies of the closest matching times
 
        [Y, ind2] = min(  abs(t -tf) );
        [Y, ind1] = min(  abs(t -t0) );
 
% Difference the coordinates 
 
        Ve  = (Ep(:,ind2) - Ep(:, ind1) );
        Vn  = (Np(:,ind2) - Np(:, ind1) );
 
 
%% Need to plot in GMT format
 
 
% Plot predicted vectors
        null = zeros(size(Ve));
        outmatrix = [llh(2,1:Nsites)', llh(1,1:Nsites)', 1000 * Ve , ...
                 1000 * Vn, null, null, null];




%% Done
~
~
~
