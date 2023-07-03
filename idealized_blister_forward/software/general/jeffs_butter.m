function [seis_filtered] = jeffs_butter(seis,dt,f1,f2,npole);
% f1 and f2 are in Hz
% dt is in seconds


nop = npole-2*(npole/2);
nepp = npole/2;
[mpts] = length(seis);
[npwr2,mpts2] = pwr2(seis(:,1));
nf = mpts2/2 + 1;

% demean and taper
seis=seis-mean(seis);
tp=taper(length(seis),.15);
data1=seis.*tp;

% set up frequency vector
df = 1.0/(mpts2*dt);
fmax = 1.0/(2.0*dt);
f = 0:df:fmax;

c0 = complex(0.,0.);
c1 = complex(1.,0.);


butter=zeros(length(f),1);

for kk=1:length(butter)

%%% make butter worth response

      w = f(kk);
      wch = f2;
      wcl = f1;
      np=0;
      if(nop>0)
        np = np+1;
        sph(np) = c1;
      end
      if (nepp>0) 
        for ii=1:nepp
          ak = 2.*sin((2.*ii-1)*pi/(2.*npole));
          ar = 0.5*ak*wch;
          ai = 0.5*wch*sqrt(4.-ak*ak);
          np = np+1;
          sph(np) = complex(-ar,-ai);
          np = np+1;
          sph(np) = complex(-ar,+ai);
        end 
      end
      np = 0;
      if (nop > 0) 
        np = np+1;
        spl(np) = c1;
      end
      if (nepp> 0)
        for ii = 1:npole/2
          ak = 2.*sin((2.*ii-1.)*pi/(2.*npole));
          ar = 0.5*ak*wcl;
          ai = 0.5*wcl*sqrt(4.-ak*ak);
          np = np+1;
          spl(np) = complex(-ar,-ai);
          np = np+1;
          spl(np) = complex(-ar,+ai);
        end
      end
      cjw = complex(0.,-w);
      cph = c1;
      cpl = c1;
      for j = 1:npole
        cph = cph*sph(j)/(sph(j)+cjw);
        cpl = cpl*cjw/(spl(j)+cjw);
      end
   butter(kk) = cph*cpl;

end

% filter seismograms
seisf = fft(seis,mpts2);

seisf(1:nf) = seisf(1:nf) .*butter;



% inverse fft
seis_filtered(nf+1:mpts2,:) = conj(flipud((seisf(2:nf-1,:))));
seis_filtered = real(ifft(seis_filtered,mpts2));
seis_filtered(mpts+1:mpts2,:) = [];
