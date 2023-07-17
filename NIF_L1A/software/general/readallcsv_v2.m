        path(path,'~/matlab/objects')
        path(path,'~/matlab/functions')
        path(path,'~/matlab/filter_tools')
	path(path,'~jeff/matlab/functions/time');
	path(path,'~jeff/matlab/matools');
%a=rdcsv('p561.2009.167.001253.csv');
%b=rdcsv('p560.2009.167.001253.csv');


	clear files
        zz=strcat('*.csv');
        f1=dir(zz);
	nf=length(f1);
        %nf=100;
        for i=1:1:nf
        files(i,1:24)=char(char(f1(i).name));
        end

% Structure Variable  
        clear S
% GPSTS arrays
        clear d dcov sites epochs epochindex apcoords masterlist
%	d is ordered station by station thankfully [e11, n11, u11, e12, n12, ...];
%       NOT TRUE?? do opposite ordering
%	dcov is sparse matrix [nsta*nepochs by nsta*nepochs] ordered?
%	masterlist is a nsta by 4 char array of station names
%	epochs is a double array of epoch times in years.
%	siteindex has the same length as d and is a numerical index of the
%		mastersitelist for each data point e.g. [1, 1, 1, 2, 2, 2, ...1, 1, 1, ]
%	epochindex is a cell array that is Nepochs long. each cell is a double array that has 
%      		integers for each data from that epoch e.g. [1,1,1,1,1,1,1,1] for epoch 1, very dumb
%	apcoords is a 3 by nsta array
%
%
%

        i=1;
	fname=files(i,1:24);
        S=rdcsv(fname);
	H=waitbar(0,'READING CSV FILES')
        masterlist(i,1:4)=files(i,1:4);
        %%lons(i)=[S.lon]; lats(i)=[S.lat]; elevs(i)=[S.ele];
 	epochs=S(i).time;       

epochs=[];
d=[];

%% concatinating structures did not work!
%% have to make d and siteindex variables as you go. 
 	for i=1:nf
          h=waitbar(i./nf);
          fname=files(i,1:24);
	  a=rdcsv(fname);
	  S=[S,a];
          masterlist(i,1:4)=files(i,1:4);
          %%lons(i)=S.lon; lats(i)=S.lat; elevs(i)=S.ele;
          epochs=[epochs; S(i).time];
          epochs=unique(epochs);
          epochs=sort(epochs);
        end
        close(H);
        whos S
	%keyboard

     Nepochs=length(epochs);
     Nsites=nf;
     Nsitesj=zeros(Nepochs,1);

%    make data
% THIS ORDERING DID NOT WORK
%     d=[]; siteindex=[];
%     for i=1:Nsites
%      Nepochsi=length(S(i).time);
%      for j=1:Nepochsi
%        thisj=find(S(i).time(j)==epochs);
%        Nsitesj(thisj)=Nsitesj(thisj)+1;
%        d=[d,S(i).e(j), S(i).n(j), S(i).u(j)];
%        siteindex=[siteindex, i, i, i];
%      end
%     end

	%H=waitbar(0,'Making Data Vector')
%% opposite ordering?
     d=[]; siteindex=[]; Nsitesj=zeros(Nepochs,1);

%%% do this in peices and then concatenate is faster?
Nis=1000;
di=round(Nepochs/Nis);
Is=1:di:Nepochs;
Is=[Is, Nepochs];

for ii=1:Nis
  ii
   i1=Is(ii);
   i2=Is(ii+1)-1;
   % this part has to be small
      ddum=[]; sitedum=[]; 
      ind=1;
   for i=i1:1:i2
     
      for j=1:Nsites
       ip=find(S(j).time==epochs(i));
       if(length(ip)==1)
        Nsitesj(i)=Nsitesj(i)+1;
        ddum(ind)=S(j).e(ip); sitedum(ind)=j;
        ddum(ind+1)=S(j).n(ip); sitedum(ind+1)=j;
        ddum(ind+2)=S(j).u(ip); sitedum(ind+2)=j;
        ind=ind+3;
       end
      end
     end
     %close(h);
% concatenate
  d=[d,ddum];
  siteindex=[siteindex, sitedum];
end


%  make covariance matrix
   dcov=sparse(length(d),length(d));

% make apcoords
% don't have a function to go backwards from llh to xyz?
% code actually gets coords from [S.lon];
    apcoords=zeros(3,Nsites);

% make epochindex
% epochindex is a cell array that is Nepochs long. each cell is a double array that has 
% integers for each data from that epoch e.g. [1,1,1,1,1,1,1,1] for epoch 1, very dumb


    for i=1:Nepochs
      epochindex{i}=i*ones(3*Nsitesj(i),1);
    end

    ts=gpsts(d,dcov,cellstr(masterlist),epochs,siteindex,cat(1,epochindex{:}),apcoords);

    save thists ts
