function [ts,ts_enu]=ReadGPSTimeSeries3(files,include,exclude,origin)
%ReadGPSTimeSeries     [ts,ts2]=ReadGPSTimeSeries3(files,include,exclude,origin)
%
%Reads GPS data from stacov or sinex files.  Outputs a gps
%time series object with actual data(ts) and a dummy one with
%the a priori constraints (ts2).
% this version outputs the enu

%Parse input arguments

    switch nargin

        case 2
            include='';
            exclude='';

        case 3
            exclude='';

        case 4

        otherwise
            help ReadGPSTimeSeries
            return
   end

    include=char(include);
    exclude=char(exclude);

%Read the files 

    if size(files,1)==1
        filelist=dir(files);
        if isempty(filelist)
            error(['No match: ',files])
        end
        index=find(files=='/' | files=='\');
        pathhead=files(1:index(end));
        files=[repmat(pathhead,size(filelist,1),1),char({filelist.name}')];
    else
        files=char(files);
    end

    nf=size(files,1);
    wb=waitbar(0,'Loading GPS Data','Name','ReadGPSTimeSeries'); 

    c=0;

    for i=1:nf
        disp(i)
        disp(nf)

        [D,DCOV,EPOCH,SITES,TYPES]=read_sinex(files(i,:));
        EPOCH

        if isempty(D)
            continue
        end

        I=strmatch('STA',TYPES);
        D=D(I);
        DCOV=DCOV(I,I);
        EPOCH=EPOCH(I,:);
        SITES=upper(SITES(I,:));
        SITES(SITES==0)='_';
        SITES=SITES(1:3:end,:);
        % whos D DCOV
	[denu,dcovenu]=xyz2enu(D,DCOV,origin);
	%disp('got here') %dcovenu allready fuckedup
	%keyboard	

        if ~isempty(include)
            [SITES,NULL,b]=intersect(include,SITES,'rows');
        elseif ~isempty(exclude)
            [SITES,b]=setdiff(SITES,exclude,'rows');
        else
            [SITES,b]=sortrows(SITES);
        end

        if ~isempty(b)
            c=c+1;
            K=keepstations(length(D)/3,b);
            d{c}=K*D;
            dcov{c}=K*DCOV*K';    
            epoch(c,1)=cal2decyr(sinex2cal(EPOCH(1,:)));
            sites{c}=SITES;
            denusv{c}=K*denu;
            denucovsv{c}=K*dcovenu*K';
        end

    waitbar(i/nf,wb);
  
    end

    delete(wb)

%Return an empty object if nothing was read

    if c==0
        ts=gpsts;
        return
    end

%Sort by epoch

    [epochs,I]=sort(epoch); 
    d=d(I);
    dcov=dcov(I);
    sites=sites(I);
    denusv=denusv(I);
    denucovsv=denucovsv(I);
	%denucovsv is fucked up allready
%Parse the data

    sites=cat(1,sites{:});
    masterlist=unique(sites,'rows');

    for i=1:length(epochs)
        epochindex{i}=repmat(i,length(d{i}),1);
    end

    d=cat(1,d{:});
    dcov=cell2blkdiag(dcov);

    denusv=cat(1,denusv{:});
    denucovsv=cell2blkdiag(denucovsv);

    I=ones(size(sites,1),1);

    for i=size(masterlist,1):-1:1
        site=masterlist(i,:);
        b=(find(~sum(site(I,:)~=sites,2))-1)*3+1;
        b=[b,b+1,b+2]';
        siteindex(b,1)=i;
        AC=d(b(:));
        apcoords(:,i)=mean(reshape(AC,3,length(AC)/3),2);
    end
    dum=[];
    for j=1:length(apcoords);
      dum=[dum, apcoords(1,j), apcoords(2,j), apcoords(3,j)];
    end
    [apcoordsenu]=xyz2enu(dum',[],origin);
    apenu2=[];
    for j=1:length(apcoords);
       apenu2(:,j)=apcoordsenu(j*3-2:j*3);
    end

%Create the time series object

    %disp('making GPSTS')
    %keyboard
    ts=gpsts(d,dcov,cellstr(masterlist),epochs,siteindex,cat(1,epochindex{:}),apcoords);
    ts_enu=gpsts(denusv,denucovsv,cellstr(masterlist),epochs,siteindex,cat(1,epochindex{:}),apenu2);
