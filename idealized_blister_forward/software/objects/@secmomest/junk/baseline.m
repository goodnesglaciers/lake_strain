%baseline junk

    
    if ~isempty(options.baseline)

       %Must filter requested indices for availability of base station

			baseindex=[];
			
			switch class(options.baseline)
			
                case {'cell','char'}
                    baseindex=strmatch(options.baseline,ts.sites);
			
                case 'double'
                    baseindex=intersect(options.baseline,1:size(ts.sites,1));
			
			end
			
			if isempty(baseindex)
                error('Could not parse base station index.')
			end

            %Force base station index to be among requested sites

                siteindex=union(baseindex,siteindex);

            %Weed out any epochs requested that don't contain the base station

                baseepochs=ts.epochindex(ismember(ts.siteindex,baseindex));
                epochindex=intersect(epochindex,baseepochs(1:3:end));

            %Weed out any epochs requested that have ONLY the base station (to do)

    end


%BEGIN MakeS --------------------------------------------------

function S=MakeS(ts,siteindex,epochindex,masterindex,baseindex)

    n=sum(masterindex);
    j=1:n;
    bi=find(ts.siteindex(masterindex)==baseindex);
    j(bi)=[];
    m=length(j);
    i=[1:m 1:m];
    j(2*m)=0;
    nobs=(diff([0;find(diff(ts.epochindex(masterindex))>0);n]))/3-1;
    bi=bi(1:3:end);
    c=m;
    for k=1:length(nobs)
        j(c+1:c+nobs(k)*3)=repmat(bi(k):bi(k)+2,1,nobs(k));
        c=c+nobs(k)*3;
    end
    S=sparse(i,j,[ones(m,1),-ones(m,1)]);

%END MakeS ----------------------------------------------------