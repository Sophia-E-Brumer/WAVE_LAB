grav=9.81; % Yes we will use SI units… this is the apparent gravity
nth=24;    % Number of directions for our discrete spectrum
nk=36;     % Number of wavenumbers : in practice we work with 
           %        spectral density in wavenumber space
XFR=1.1;   % ratio of frequencies f(i+1)/f(i) 

d2r=pi./180;

% Now we define our vector of discrete frequencies 
freq=0.034.*exp(linspace(0,nk-1,nk).*log(XFR))';
sig=2.*pi.*freq;
% and corresponding wavenumbers 
k=sig.^2/grav; 

% and a few more arrays to manipulate these
df=freq.*(XFR-1/XFR)./2;
Cg=grav./(2.*sig);      % group speed assuming deep water
dk=df.*2.*pi./Cg;
dk2d=repmat(dk,1,nth);
sig2d=repmat(sig,1,nth);
k2d=repmat(k,1,nth);
c2d=sig2d./k2d;
dkdf2d=repmat((dk./df),1,nth);

% now we define the discrete directions 
dth=360./nth;		% direction increment in degrees
dir=linspace(0,(nth-1).*dth,nth);
dthr=dth.*d2r;          % direction increment in radians



% We will now also discretize the time ... 
% You can increase tmax from 12 to 96 hours ... 

tmax=48*3600;		% time at end of integration (in seconds)
dt=100;nt=tmax/dt+1;	% time step and number of steps

E=zeros(nk,nth);
allE=zeros(nk,nth,nt); % stores the spectrum at all time steps
allEf=zeros(nk,nt);    % stores the spectrum at all time steps
Hs=zeros(nt,1);

% defining our winf forcing
% WARNING: when you change the wind in time, you have to move this piece 
% inside of the loop on time steps
U10=10;
thetaw=180;
Cd=0.0012;                     % drag coefficient
ustar=sqrt(Cd).*U10;           % NB: the wind stress in rhoa*ustar^2 = rhoa*Cd*U10^2
fPM=0.123.*grav./(28.*ustar);  % this is the Pierson-Moskowitz frequency
coswind=cos((dir-thetaw).*d2r);
cos2d=repmat(cos((dir-thetaw).*d2r),nk,1);

% some parameters
rhoa=1.29;  % air density in kg/m^3
rhow=1026;  % water density in kg / m^3

coef1=80.*(Cd.*rhoa./rhow).^2./grav.^2;



coefSnyder=0.25.*rhoa./rhow;
coefKomen=-2.36E-5;
s2PM=3.02E-4;

% If the wind stays constant, then the linear wave growth is constant
% and can be pre-computed outside of the time loop
for ik=1:nk
  G=exp(-(freq(ik)./fPM).^-4); % filter for linear input
  %G=1; % use G=1 to test the impact of direct input to low frequencies
  Sinl(ik,:)=coef1.*sig(ik)./k2d(ik).*U10.^4.*max(0,coswind).^4.*G;
end


time=linspace(0,dt*(nt-1),nt);

% Loop on time steps. 
for it=1:nt-1  
  dt=(time(it+1)-time(it));
  allEf(:,it)=sum(E(:,:).*dkdf2d(:,:),2).*dthr;
  
  Etot=sum(allEf(:,it).*df(:)); 
  Eftot=sum(allEf(:,it).*df(:)./freq); 
   

  fmean=min(freq(end),(Etot+1E-10)./(Eftot+3E-10));
  kmean=(2.*pi.*fmean).^2./grav;
  s2mean=Etot.*kmean.^2./grav.^2;
  
  Hs(it)=4.*sqrt(Etot);
 
  Sin_Snyder=coefSnyder.*max(0,28.*ustar./c2d.*cos2d-1);
  Sds_Komen=coefKomen.*fmean.*(s2mean./s2PM).^2.*k2d./kmean;
  
  % Only linear growth ... 
  %E(:,:)=E(:,:)+Sinl(:,:)*dt;
   
   diag= Sin_Snyder+Sds_Komen;
   E(:,:)=E(:,:)+E(:,:).*diag.*dt+Sinl(:,:)*dt;


end
% values for last time step
allEf(:,it+1)=sum(E(:,:).*dkdf2d(:,:),2).*dthr;
Etot=sum(allEf(:,it+1).*df(:)); 
Hs(it+1)=4.*sqrt(Etot);



figure(1)
clf
plot(time./3600,Hs,'LineWidth',3);
set(gca,'FontSize',16,'LineWidth',1);
title('Significant wave height (m)')
xlabel('time (hours)')
grid on


figure(2)
clf
plot(freq,allEf(:,1:7200/dt:end)')
set(gca,'FontSize',16,'LineWidth',1);
title('Evolution of frequency spectrum')
xlabel('frequency (Hz)')
ylabel('E(f) (m^2/Hz)')
hold on
plot(freq,allEf(:,end),'k-','LineWidth',2)

grid on

figure(3)
clf
pcolor(freq,dir,10.*log10(E(:,:).*dkdf2d)');
title('10 log10(E(f,theta))');
xlabel('frequency (Hz)');
ylabel('direction (deg, to)');
colorbar
caxis([-40 0])

figure(4)
clf
sinf=sum(Sin_Snyder(:,:).*E(:,:).*dkdf2d.*dthr,2);
sdsf=sum(Sds_Komen(:,:).*E(:,:).*dkdf2d.*dthr,2);
plot(freq,sinf,'r-',freq,sdsf,'b-',freq,sinf+sdsf,'g-','LineWidth',2);
set(gca,'FontSize',16,'LineWidth',1);
legend('Sin','Sds','Stot');
xlabel('frequency (Hz)');
ylabel('Source term (m^2/Hz/s)');
grid on
 

figure(5)
clf
pcolor(freq,dir,10.*log10(Sin_Snyder)');
title('wind input growth rate: 10 log10(Sin/E)');
xlabel('frequency (Hz)');
ylabel('direction (deg, to)');
colorbar
caxis([-50 -20])



