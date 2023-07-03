
function PrintResids(r, d_gps, d_lev, G)

 
        disp(' ')
        disp('GOODNESS OF FIT: All Data')
        fprintf('Weighted Residual Sum of Squares (WRSS), %g\n', r'*r)
        fprintf('WRSS / (N-P),         %g\n', (r'*r)/(length(r)-rank(G)) )
        disp('__________________________________ ')
        disp(' ')
        disp(' ')

	n_gps = length(d_gps);
if n_gps > 0
        r_gps = r(1:n_gps);
	rank_gps = rank(G(1:n_gps,:));
        disp('GOODNESS OF FIT:  GPS')
        fprintf('GPS:   Weighted Residual Sum of Squares (WRSS), %g\n', r_gps'*r_gps)
        fprintf('GPS:   WRSS / (N-P),         %g\n', (r_gps'*r_gps)/(n_gps-rank_gps) )
        disp('__________________________________ ')
        disp(' ')
        disp(' ')
end

	n_lev = length(d_lev);
if n_lev > 0
        r_lev = r(n_gps+1:n_gps + n_lev);
	rank_lev = rank(G(n_gps+1:n_gps + n_lev,:));
        disp('GOODNESS OF FIT:  Leveling')
        fprintf('Lev:   Weighted Residual Sum of Squares (WRSS), %g\n', r_lev'*r_lev)
        fprintf('Lev:   WRSS / (N-P),         %g\n', (r_lev'*r_lev)/(n_lev-rank_lev) )
        disp(' ')
end
