	function [d_pseudo, Hs]  = smoothH(t, sqrtL, Nbasis, statedim)

%%
%  [d_pseudo, Hs]  = smoothH(t, sqrtL, Nbasis, statedim);
%
%  funtion to implement spatial smoothing as pseduo-observations
%  Note there is some question as to whether we should be smoothing
%  velocity or slip, or whether it in fact matters.
%  For now smooth velocity
%
% P.Segall 12/1/1997

	A = zeros(Nbasis, statedim);
	d_pseudo = zeros(Nbasis,1);


        for i = 1:Nbasis
                 A(i,3*(i-1)+1:3*(i-1)+3) = [1,0,1];
%%               A(i,3*(i-1)+1:3*(i-1)+3) = [t,1,0];
%% This would be to smooth the slip rather than velocity
        end

	Hs = sqrtL*A;
