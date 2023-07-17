  function Q =  makeQ_sparse(statedim, Nbasis, Nframe,Nsubf1, Nsubf2, Nsubf3, delt, alpha1, alpha2, alpha3, tau, sig_frame);

%%
%% Q = makeQ2(statedim, Nbasis, Nframe, 
%%			delt, alpha, tau, sig_frame)
%%
%%  subroutine to make "Process" Covariance for Kalman filter
%%  Q is (statedim x statedim)
%%
%%  Input:
%%	statedim  = dimension of state vector
%%	Nbasis	  = number of basis functions
%%	alpha     = scale parameter for integrated random walk
%%	tau       = scale parameter for random walk
%%	sig_frame = scale parameter for white noise 
%%	delt      = t_k - t_{k-1} at present epoch
%%	
	statedim;

	Q = zeros(statedim);
	Ncomps=statedim-Nbasis*2-Nframe;
    
	for i = 1:Nsubf1
		index = 2*(i-1)+1;
        Q(index,index) = alpha1^2*delt^3/3;
        Q(index,index+1) = alpha1^2*delt^2/2;
        Q(index+1,index) = Q(index,index+1);            
        Q(index+1,index+1) = alpha1^2*delt;     
    end

    for i = Nsubf1+1:Nsubf1+Nsubf2
		index = 2*(i-1)+1;
        Q(index,index) = alpha2^2*delt^3/3;
        Q(index,index+1) = alpha2^2*delt^2/2;
        Q(index+1,index) = Q(index,index+1);            
        Q(index+1,index+1) = alpha2^2*delt;     
        end
    
    for i = Nsubf1+Nsubf2+1:Nsubf1+Nsubf2+Nsubf3
		index = 2*(i-1)+1;
        Q(index,index) = alpha3^2*delt^3/3;
        Q(index,index+1) = alpha3^2*delt^2/2;
        Q(index+1,index) = Q(index,index+1);            
        Q(index+1,index+1) = alpha3^2*delt;     
	end    
    
    
    %keyboard
		
	Ncomps = (statedim-2*Nbasis-Nframe);
    %whos tau delt
	for i = 1:Ncomps
		index=2*Nbasis+i;
		Q(index, index) = tau^2*delt;
	end
		    
		
		
% process noise for frame error

	Ntrans = 3; Nrot = 3; Nscale = 1;

	for i = 1:Ntrans
		index=2*Nbasis+Ncomps+i;
		Q(index,index) = sig_frame(i)^2;
    end

	Q = sparse(Q);
