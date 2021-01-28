function  [Hm0,Tp,imax,km,Ef,overlap,df,Hsws]=HsTpkm_from_spectrum_windsea(E,freq,dir);
g=9.81;
nk=size(freq,1);
nth=size(dir,1);
dth=2*pi/double(nth);
df=zeros(nk,1);
df(2:nk-1)=0.5*(freq(3:nk)-freq(1:nk-2));
df(1)=freq(2)-freq(1);
df(nk)=freq(nk)-freq(nk-1);
Ef=sum(E,2)*dth;
overlap=dth.*sum(E(:,1:nth/2).*E(:,1+nth/2:nth),2)./(Ef.^2);
   

[Emax,imax]=max(Ef);
fp=sum(Ef(max(imax-1,1):min(imax+1,nk)).*freq(max(imax-1,1):min(imax+1,nk))) ...
    ./sum(Ef(max(imax-1,1):min(imax+1,nk)));
Etot=(Ef'*df(:));
km=(Ef'*df(:)*freq(:))./Etot.*(2*pi/sqrt(g));
Hm0=4.*sqrt(Etot);
Tp=1./fp;
I=find(freq > 0.068);
Etot=(Ef(I)'*df(I));
Hsws=4.*sqrt(Etot);
