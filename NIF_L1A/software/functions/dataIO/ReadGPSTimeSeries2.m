function [ts,ts2]=ReadGPSTimeSeries2(files,include,exclude,stacov)
%ReadGPSTimeSeries     [ts,ts2]=ReadGPSTimeSeries2(files,include,exclude)
%
%Reads GPS data from stacov or sinex files.  Outputs a gps
%time series object with actual data(ts) and a dummy one with
%the a priori constraints (ts2).

%Parse input arguments

    switch nargin

        case 1
            include='';
            exclude='';

        case 2
            exclude='';

        case 3

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
        
        [D,DCOV,EPOCH,SITES,TYPES,CON]=read_sinex(files(i,:));
        
        if isempty(D)
            continue
        end

        I=strmatch('STA',TYPES);
        D=D(I);
        DCOV=DCOV(I,I);
	CON=CON(I,I);
        EPOCH=EPOCH(I,:);
        SITES=upper(SITES(I,:));
        SITES(SITES==0)='_';
        SITES=SITES(1:3:end,:);


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
	    apcon{c}=K*CON*K'; %sparse?  
            epoch(c,1)=cal2decyr(sinex2cal(EPOCH(1,:)));
            sites{c}=SITES;
        end
	whos dcov apcon

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
    whos I dcov apcon
    dcov=dcov(I);
    apcon=apcon(I);
    whos I dcov apcon
    sites=sites(I);

%Parse the data

    sites=cat(1,sites{:});
    masterlist=unique(sites,'rows');

    for i=1:length(epochs)
        epochindex{i}=repmat(i,length(d{i}),1);
    end
	%keyboard
    d=cat(1,d{:});
    dcov=cell2blkdiag(dcov);
    apcon=cell2blkdiag(apcon);  %PROBLEM HERE!
    
    whos d dcov apcon

    I=ones(size(sites,1),1);

    for i=size(masterlist,1):-1:1
        site=masterlist(i,:);
        b=(find(~sum(site(I,:)~=sites,2))-1)*3+1;
        b=[b,b+1,b+2]';
        siteindex(b,1)=i;
        AC=d(b(:));
        apcoords(:,i)=mean(reshape(AC,3,length(AC)/3),2);
    end

%Create the time series object

    ts=gpsts(d,dcov,cellstr(masterlist),epochs,siteindex,cat(1,epochindex{:}),apcoords);
    whos d apcon dcov epochs masterlist siteindex epochindex apcoords
    ts2=gpsts(d,apcon,cellstr(masterlist),epochs,siteindex,cat(1,epochindex{:}),apcoords);
