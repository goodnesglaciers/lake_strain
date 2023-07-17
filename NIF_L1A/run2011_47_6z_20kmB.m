%% RUNFILE FOR 2011 NIF
%  March 22 2018  run the 2011 NIF with only 4 stations - LAS
%  November 20 2020: run the 2011 NIF with all stations and a 20km by 20 km basal patch.
%  Use with makegeom_caseC_20km.m for 20km basal plane + L1A hydro-fracture
%  Runs NIF and saves output as: MLE_47_all_6z_2_20kmB.mat

path(path,'software/general')
path(path,'software/functions/dataIO')
path(path,'software/objects')
path(path,'software/functions/elmat')
path(path,'software/functions/time')
path(path,'software/filter_tools')
path(path,'software/matools')
path(path,'software/okada85')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load a .mat file that has tsb and tsb_enu in them which are time series
% objects from readsinex.m that are in GPS and ENU coordinates respectively
% preferably already downsampled in time to make this go quicker.

         load ts_displ_week_NLBS.mat % load time series object
         ts_displ = ts_displ_week_NLBS;
         
         tstart=167.1; % start time for inversion
         tend=172.0; % end time for inversion
         [a,I1]=min(abs(ts_displ.epochs-tstart));
         [a,I2]=min(abs(ts_displ.epochs-tend));
         
% Select epochs by subsetting 
% May 24 2018: Changed subsetting from '1' to '3' if downsampling timeseries object
         ts2=ts_displ(:,I1:1:I2); % change to 1 for MLE runs
         
         figure % number of observations per time
         for i=1:length(ts2.epochs)-1
             plot(i,length(ts2(:,i).d),'*')
             hold on;
         end
         xlabel('epoch #')
         ylabel('number of observations')
         

%%%%%%%% SELECT STATIONS %%%%%%%%%%%%%%%%%%%%%%%
tsb_enu=ts2;   %can't do anything here until you get lats+lons in here instead of hardwired.
%tsb_enu=tsb_enu(['NL01'; 'NL02'; 'NL03'; 'NL04';   ],:);
%tsb=tsb(['NL01'; 'NL02'; 'NL03'; 'NL04'; ],:);

load apcoords_lle.mat;
lons=apcoords_lle(2,:);
lats=apcoords_lle(1,:);
Nsites=length(lats);
plotyn=1;
        if(plotyn)
         figure
         plot(lons,lats,'ko','MarkerFaceColor','r','MarkerSize',4)
         hold on;
         for j=1:length(lats)
           text(lons(j),lats(j),tsb_enu.sites(j));
         end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PICK EPOCHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        J=find(tsb_enu.epochs>=2011.3 & tsb_enu.epochs<=2011.6);
%        J=J(1:2:end);       
%        ts_enu=tsb_enu(:,J); ts=tsb(:,J);  %subsetting

ts_enu=tsb_enu;   %already subsetted above

        % uncomment below if you think you have some bad solutions to throw
        % out a priori; fix threshold to suit
        %disp('finding strange epochs')
	    %for i=1:length(ts_enu.epochs);
    	%	[U,S,V]=svd(ts_enu(:,i).dcov);
     	%	mineig(i)=min(diag(S));
	    %end
        %KK=find(mineig>=1e-10);
        %ts_enu=ts_enu(:,KK);
        %disp('done')


%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP FAULT GEOMETRY%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
dosetup=1;
if(dosetup)
  origin=[68.72, -49.53];  %% SET COORDINATE SYSTEM ORIGIN IN MIDDLE OF ARRAY
  % we used to send down the time series object in GPS XYZ coordinates, but
  % this dataset doesn't have that.   
  [G1, G2, G3, G2b, G3b, Rc, Rb, patchesC, patchesB]=makegeom_caseC_20km(ts_enu,origin);
  %keep different Rs, gammas for crack and base
%
% make Kern  Nsta*3 by Nsubfaults  ordered Sta1ENU, Sta2ENU,...
justcrack=0;
if(justcrack)
 Nsubf=size(G3,2);
 KernO=zeros(Nsites*3,Nsubf);
 % Kernal for opening is only G3
 for i=1:Nsites
  ind=(i-1)*3+1;
  KernO(ind,1:Nsubf)=G3(i,1:Nsubf);   %E comps
  KernO(ind+1,1:Nsubf)=G3(Nsites+i,1:Nsubf);  %N comps
  KernO(ind+2,1:Nsubf)=G3(Nsites*2+i,1:Nsubf); %Z comps
 end
 
else
 Nsubfa=size(G3,2);
 Nsubfb=size(G3b,2);
 Nsubfc=size(G2b,2);
 Nsubf=Nsubfa+Nsubfb+Nsubfc;
 KernO=zeros(Nsites*3,Nsubf);
 % Kernal for all 3 is ordered  Crack Opening 1-Na; Basal openin 1-Nb; Basal thrust 1-Nc;
 for i=1:Nsites
  ind=(i-1)*3+1;
  KernO(ind,1:Nsubfa)=G3(i,1:Nsubfa);   %E comps
  KernO(ind+1,1:Nsubfa)=G3(Nsites+i,1:Nsubfa);  %N comps
  KernO(ind+2,1:Nsubfa)=G3(Nsites*2+i,1:Nsubfa); %Z comps
  % basal openin
  KernO(ind,Nsubfa+1:Nsubfa+Nsubfb)=G3b(i,1:Nsubfb);   %E comps
  KernO(ind+1,Nsubfa+1:Nsubfa+Nsubfb)=G3b(Nsites+i,1:Nsubfb);  %N comps
  KernO(ind+2,Nsubfa+1:Nsubfa+Nsubfb)=G3b(Nsites*2+i,1:Nsubfb); %Z comps
  % basal thrust
  KernO(ind,Nsubfa+Nsubfb+1:Nsubf)=G2b(i,1:Nsubfc);   %E comps
  KernO(ind+1,Nsubfa+Nsubfb+1:Nsubf)=G2b(Nsites+i,1:Nsubfc);  %N comps
  KernO(ind+2,Nsubfa+Nsubfb+1:Nsubf)=G2b(Nsites*2+i,1:Nsubfc); %Z comps  
 end
end
    
  % load in a priori velocities
  %bring in V0_enu, X0_enu, V0_cov_enu, ref_epoch, sites0_enu;
  %load apvelenu_every4;   
  %zz=reshape(X0_enu,3,29); 
  %X0_enu=zz; 
  %sitesV=char(sites0_enu); 
else
 load setup_v1.mat
end
disp('done setup')

Kern=KernO;  % Here we concatonate all of the source we want to use into a Nsites*3 by Nbasis matrix
%keyboard

%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsites=length(ts_enu.sites);  
Nbasis=size(Kern,2)
Nframe=3;  sig_f=[.01, .01, .05]; % Just translation for this array; Assumed White Noise in m^2;  Variances set here.
Nepochs=length(ts_enu.epochs);

% SECULAR VELOCITIES
zerovs=0;   clear V0_enu V0_cov_enu X0_enu sitesV ref_epoch;
sitesV=ts_enu.sites;
ref_epoch=ts_enu.epochs(1);
V0_cov_enu=zeros(3*Nsites,3*Nsites);
tstart=167.0997; tend=168.8309;
[a,emin]=min(abs(ts_enu.epochs-tstart));
[a,emax]=min(abs(ts_enu.epochs-tend));

if(zerovs)
    V0_enu=zeros(3*Nsites,1);  % Do we have a prior velocities that matter? 
    V0_cov_enu=zeros(3*Nsites,3*Nsites);
    X0_enu=zeros(3*Nsites,1);  % This is how Laura defines her time series? <-- Yes (LAS)
else
  % simple station by station fit?
  % emax=250; emin=1;
  for j=1:Nsites
    for k=1:3
%       if(length(ts_enu(j,emin:emax).d)>=10)
        test = (ts_enu(j,emin:emax).epochs);
        X=test-ref_epoch;
        Y=ts_enu(j,emin:emax).d(k:3:end);

         [P,S]=polyfit(X,Y',1);
         R=S.R; Rinv=inv(R); df=S.df; normr=S.normr;
         covP=(Rinv*Rinv')*normr^2/df;
         
%         figure
%         plot(X,Y,'*');hold on; plot(X,P(1)*X+P(2),'r');
         
         ind=(j-1)*3+k;
%         else
%          P(1)=V0_enu((j-2)*3+k);
%          P(2)=X0_enu((j-2)*3+k);
%          covP(1,1)=.01;   % hardwired big
%       end
%       ind=(j-1)*3+k;
        V0_enu(ind)=P(1);
        X0_enu(ind)=P(2);
        V0_cov_enu(ind,ind)=covP(1,1);
        %keyboard
     end
  end
    
end

%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% STATEVECTOR  %%%%%%%%%%
%%%%% 3*2*Nfaults + Nsites +XREF + 5 (hyper params).
% This is NIF (no positivity wanted) not eNIF so no hyper parameters in state vector
% Using vertical 3 benchmark terms per station
% No positivity = no dummy variable   = only 2 elements (RW+IRW) per basis
% function
% State Vector is ordered:  [Basis coefficients (IRWi, RWi...);  L;  Frame;]
% Do we have any offsets to add in?
statedim=2*Nbasis +3*Nsites(1) +Nframe;
%V0_enu=zeros(3*Nsites,1);  % Do we have a prior velocities that matter? 
%V0_cov_enu=zeros(3*Nsites,3*Nsites);
%X0_enu=zeros(3*Nsites,1);  % This is how Laura defines her time series?
%sitesV=ts_enu.sites;
%ref_epoch=ts_enu.epochs(1);
%
% PRIORS 
     clear x0 var0;
     x0(1:statedim)=0;
     var0=1e-8*eye(statedim);
     Ncomps=3*Nsites(1)
            
     % And constrain Slip to be zero at start but slip-rate to be
     % undetermined
     for j=1:Nbasis
         x0(j*2)=0;
         x0(j*2-1)=0;
         var0(j*2,j*2)=100^2;   % Rate? What units is this? m/day?
         var0(j*2-1,j*2-1)=.0001^2;   % Slip in units of m?
     end
                                                                                                  
     % Give the L values freedom in the first step
     for j=1:Ncomps
       x0(j+2*Nbasis)=0;
       var0(j+2*Nbasis,j+2*Nbasis)=.0001^2;   % in case the defenition of zero is wrong
     end
            
     % LOOSE CONSTRAINT FOR FRAME ERROR
      for j=2*Nbasis+Ncomps+1:statedim
         var0(j,j)=.2^2;
      end
 
     %% HYPER-PARAMETER PRIORS!!
%         tau=.1   % m/yr^1/2
%         sigma=2;  % never trust GPS covariances
%         alphaC=500;
%         alphaB=50;
%         gammaC=20;  %20 is good for 8 by 4
%         gammaB=20;
%         hypers=[tau,sigma,alphaC,alphaB,gammaC,gammaB];
%         
      %%RANGE HYPER-PARAMETER PRIORS!!
       % tau= 0.01 %logspace(-2,1,10); %linspace(0.01, 1.1, 2)   % m/yr^1/2 (random walk)
        %sigma= 2; %linspace(1, 4, 2);  % never trust GPS covariances (white noise)
       % alphaC=linspace(5,10005,5); %linspace(3000,10000, 3); % temporal smoothing 
       % gammaC=linspace(5,10005,5); %linspace(5, 35, 3);
       % alphaB=linspace(1,1001,5); %linspace(1000, 10000, 3);  %20 is good for 8 by 4
       % gammaB=linspace(1,1001,5); %linspace(5, 55, 3); % spatial smoothing
       % hypers=vertcat(tau,sigma,alphaC,alphaB,gammaC,gammaB);

        %%RANGE HYPER-PARAMETER PRIORS!!
        tau= 0.05 %logspace(-2,1,10); %linspace(0.01, 1.1, 2)   % m/yr^1/2 (random walk)
        sigma= 2; %linspace(1, 4, 2);  % never trust GPS covariances (white noise)
        alphaC=200; %linspace(3000,10000, 3); % temporal smoothing 
        gammaC=450; %linspace(5, 35, 3); % spatial smoothing
        alphaB=25; %linspace(1000, 10000, 3);  %20 is good for 8 by 4
        gammaB=50; %linspace(5, 55, 3); % spatial smoothing
       % hypers=vertcat(tau,sigma,alphaC,alphaB,gammaC,gammaB);

       theta1= tau./sigma;
       theta2C=alphaC./sigma;
       theta3C=gammaC./sigma;       
       theta2B=alphaB./sigma;
       theta3B=gammaB./sigma;
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LAPLACIAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1=eye(Nsites*3); R2=R1;  %just damping total moment
R1=Rc;
R2=Rb;
c=1; c2=1;  % what are these?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Ready To Start Filter with HYPER PARAMETERS MAXIMUM LIKELIHOOD')
%keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%%%%%%%%%% HYPER PARAMETERS MAXIMUM LIKELIHOOD  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---------------------------------------------------------------

%% netfilt_withvert
%%  input:      theta2   = vector of trial alpha/sigma;
%%              theta3   = vector of trial gamma/sigma;
smooth=1; 
% set up hyperparameters to loop through        
        ntheta2C=length(theta2C);
        ntheta3C=length(theta3C);
        ntheta2B=length(theta2B);
        ntheta3B=length(theta3B);
        ntheta1=length(theta1);
         
        val = zeros(ntheta2C, ntheta3C, ntheta2B, ntheta3B, ntheta1);
        sig2 = zeros(ntheta2C, ntheta3C, ntheta2B, ntheta3B, ntheta1);

% % %%%%%%%%%%%%%%%%%%%%% NOW RUN FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

%for iii=4 % individual run for each value in fuxi
%    for jjj=1
 %       for kkk=1
  %          for lll=1
[val, sig2, x_s, xcov] = ...
	  netfilt_withvert_MLEiii14(ts_enu, Kern, x0, var0, sitesV, X0_enu, V0_enu, ...
      V0_cov_enu, ref_epoch, origin, R1, R2, smooth, c, sig_f,...
      theta2C, theta3C, theta2B, theta3B, theta1)

save junk47_2_6z_20kmB i
save val_47_2_6z_20kmB val
save sig2_47_2_6z_20kmB sig2
disp('One Iteration Complete')
          %  end
 %       end
  %  end
%end
      
disp('Done Filter')

%%%%%%%%%%%%DETERMINE BEST ESTIMATE HYPERS from minimum values
%5save MLE_47.mat

 %%%%%%%%%% EXTRACT COMPONENTS OF STATE VECTOR %%%%%%%%%%%%%%%%%%%%
 %%%%% remember statedim=2*Nbasis +3*Nsites(1) +Nframe;
 
clear L Lf F slip_b slip_f rate_f rate_b Cslip alphavar gammavar

 for i=1:Nepochs;
L(:,i)=x_s{i,end}(Nbasis*2+1:Nbasis*2+Nsites*3);
 Lf(:,i)=x_s{i,i}(Nbasis*2+1:Nbasis*2+Nsites*3);
 slip_b(:,i)=x_s{i,end}(1:2:Nbasis*2);
slip_f(:,i)=x_s{i,i}(1:2:Nbasis*2);
rate_f(:,i)=x_s{i,i}(2:2:Nbasis*2);
rate_b(:,i)=x_s{i,end}(2:2:Nbasis*2);
  F(1,i)=x_s{i,end}(statedim-2);
  F(2,i)=x_s{i,end}(statedim-1);
 F(3,i)=x_s{i,end}(statedim);
 Ff(1,i)=x_s{i,i}(statedim-2);
  Ff(2,i)=x_s{i,i}(statedim-1);
 Ff(3,i)=x_s{i,i}(statedim);
 
 end
% % % % 
 L0=L(:,1);  % L0 DOESN"T MEAN ANYTHING FOR STATIONS THAT DON"T HAVE DATA ON FIRST EPOCH
% % % % 
% % % % % displacements due to transient fault motion
 transhat_f=KernO*slip_f;
 transhat_b=KernO*slip_b;
% % % % % displacements due to secular velocity
 t=ts_enu.epochs-ref_epoch;
Vhat=X0_enu.'*ones(1,Nepochs)+V0_enu'*t';

save MLE_47_all_6z_2_20kmB.mat
% % % % %
disp('Done saving output')       
