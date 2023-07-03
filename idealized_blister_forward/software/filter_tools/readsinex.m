function [d,dcov,epochs,sites,masterlist]=readsinex(files,SITES,refsite)
%readsinex     [d,dcov,epochs,sites,masterlist]=readsinex(files,SITES,refsite)
%
%Reads through multiple sinex files. Optional SITES lists sites to keep.
%Optional 'refsite' causes differencing with the specified reference site.

if nargin < 2
	SITES=[];
end
if nargin < 3
      refsite=[];
elseif ~isempty(SITES)
   SITES=union(SITES,refsite,'rows');
end

sites=[];
d=[];
dcov=sparse(0,0);
epochs=[];

%Loop through files

	for i=1:size(files,1)
		
		if iscell(files)
			FILENAME=files{i};
		else
			FILENAME=files(i,:);
		end

	%Check to see if file exists

		if size(dir(FILENAME),1)==0
			disp(['Could not find ',FILENAME,'.'])
		else
		
			[x,xcov,obstime,sitelist]=read_sinex(FILENAME);
			obstime=obstime(1);
			if ~isempty(SITES)
            [S,a,b]=intersect(SITES,sitelist,'rows');
            K=keepstations(size(x,2),b);
            if ~isempty(b)
               x=x(:,b);
               xcov=K*xcov*K';
            end
      	else
            S=sitelist;
      	end
			if ~isempty(S)
            if ~isempty(refsite)
               J=strmatch(refsite,S);
               if ~isempty(J)
                  s=SubtractionMatrix(J,size(x,2));
                  x=s*x(:);
                  xcov=s*xcov*s';
                  S(J,:)=[];
               else
                  disp(['Specified reference site was not found in ',FILENAME,'.'])
                  x=x(:);
               end
            else
               x=x(:);
            end
            d=[d;x(:)];
            j=size(dcov,1);
            k=size(xcov,1);
            dcov(j+1:j+k,j+1:j+k)=xcov;
            sites=[sites;S];
				epochs=[epochs;repmat(obstime,k,1)];
			else
			   disp(['None of the stations specified were found in ',FILENAME,'.'])
         end
		end
	end
   masterlist=unique(sites,'rows');

