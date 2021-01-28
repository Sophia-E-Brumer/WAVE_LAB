clear all;
close all;
grav=9.81; % Yes we will use SI units… this is the apparent gravity
d2r=pi./180;

% defines folders and names of the different runs

tag1='WAM4.5';path1='ST3-Default';
tag2='T471';path2='STY-T471';
tag3='T475';path3='ST4-T475';
tag4='T700 (Romero 2019)';path4='ST4-T700';  %tag2='T600h';path2='ST4-T600h';
tag5='T700NL2'; path5='STYNL2';
%tag2='T705';path2='ST4-T705';
%tag5='T702';path5='ST4-T702';
tag6='T704';path6='ST4-T704';


nfi=5;nfip=5; 			    % number of different runs for plots
syms1=['ro','ks','b<','mx','m+'];   % symbols for plots
syms=reshape(syms1,2,nfi);

name1='/ww3.P1_2000_src.nc';

% some interesting frequencies ... 

f9=0.2;
f10=0.5;
f11=0.8;
f12=0.95;

% Creates a few figures ... 

figure(1)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
title('Wave height growth')
grid on

figure(2)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
grid on

figure(3)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
grid on

figure(4)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
grid on

figure(5)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
set(gca,'Yscale','log','Xscale','log');
grid on

figure(6)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
grid on

figure(7)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
grid on

figure(8)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
grid on

figure(13)
clf
hold on
set(gca,'FontSize',16,'LineWidth',2);
grid on


for ifig=6:15
figure(ifig)
clf
hold on
grid on
set(gca,'FontSize',16)
end

% loop on the model runs 
for ifi=1:nfi
   eval(['path=[ path' num2str(ifi) '];'])
   file=[ path name1];
   % Now reads spectrum and source terms from the current model run
   [lat,lon,freq,dir,df,time,Efth,depth,U10,UC,Cdir,unit1]=readWWNC_SPECv2(file,'efth');
   [lat,lon,freq,dir,df,time,Sin,depth,U10,UC,Cdir,unit1]=readWWNC_SPECv2(file,'sin');
   [lat,lon,freq,dir,df,time,sds,depth,U10,UC,Cdir,unit1]=readWWNC_SPECv2(file,'sds');
   [lat,lon,freq,dir,df,time,snl,depth,U10,UC,Cdir,unit1]=readWWNC_SPECv2(file,'snl');
   dates=time;
   k=((2.*pi.*freq).^2./grav); % deep water only
   % k=dispNewtonTH(freq,100); % this would be some known water depth ... 
  
   taille=size(Efth);
   nk=taille(2);
   nth=taille(1);
   nt=taille(4);

   dth=2*pi./nth;
   cos2=repmat(cos(dir'.*pi./180).^2,nk,1);
   sin2=repmat(sin(dir'.*pi./180).^2,nk,1);


   % Reads other parameters from the model run (wind forcing, wind stress ...) 
   filename=[path '/ww3.200001.nc'];
   [lat,lon,time,ust,varu,unitu]=read_WWNC_var(filename,'uust');
   ust2=squeeze(ust);
   [lat,lon,~,hs,varu,unitu]=read_WWNC_var(filename,'hs');
   hs2=squeeze(hs);
   [lat,lon,time,taw,varu,unitu]=read_WWNC_var(filename,'utaw');
   tw2=squeeze(taw);
   varmss1='mssu';  % NB: these mss (as of today when no tail is added) can also be computed from the spectrum
   varmss2='mssc';
   [lat,lon,time,mssx,varu,unitu]=read_WWNC_var(filename,varmss1);
   [lat,lon,time,mssy,varu,unitu]=read_WWNC_var(filename,varmss2);
   mss2=squeeze(mssx+mssy);
   [lat,lon,time,wndu,varu,unitu]=read_WWNC_var(filename,'uwnd');
   [lat,lon,time,wndv,varu,unitu]=read_WWNC_var(filename,'vwnd');
   wnd=sqrt(wndu.^2+wndv.^2);
   wnd2=squeeze(wnd);
   ndateso=taille(4);
   E=zeros(nk,nth);

   %
   % Extracts one spectrum
   %
   Hs=zeros(ndateso,1);
   Tp=zeros(ndateso,1);
   for it=1:ndateso;
      E(:,:)=Efth(:,:,1,it)';
      dth=2*pi./nth;
      [Hs(it),Tp(it),~,~,~,~,~]=HsTp_from_spectrum_windsea(E,freq,dir);
      if (it == 49)  % Chooses the 49th time step
        Efth2=E;
        [~,~,imax,km,Ef,overlap,df]=HsTp_from_spectrum_windsea(E,freq,dir);

      end
   end
    
   times=(dates-dates(1))*24;
   timenc=(time-time(1))*24;

   figure(1)
   plot(timenc,hs2,syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(2)
   plot(timenc,mss2,syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(3)
   plot(timenc,ust2.^2./(wnd2.^2),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(4)
   plot(timenc,squeeze(mssy)./squeeze(mssx),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(5)
   plot((freq),(Ef),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(6)
   plot((freq),(Ef).*(freq.^5),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(7)
   plot((freq),(Ef).*(freq.^4),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(8)
   plot(times,wnd2(1)./(9.81.*Tp./(2*pi)),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
   figure(13)
   plot(freq,overlap,syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)


figure(9)
I=find(freq > f9);
plot(dir,Efth2(I(1),:),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
axis([0 360 0 1])
set(gca,'XTick',linspace(0,360,25))
figure(10)
I=find(freq > f10);
plot(dir,Efth2(I(1),:),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
axis([0 360 0 0.009])
set(gca,'XTick',linspace(0,360,25))
figure(11)
I=find(freq > f11);
plot(dir,Efth2(I(1),:),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
axis([0 360 0 1.3E-3])
set(gca,'XTick',linspace(0,360,25))
figure(12)
I=find(freq > f12);
plot(dir,Efth2(I(1),:),syms(ifi*2-1:ifi*2),'LineWidth',2,'MarkerSize',8)
axis([0 360 0 8E-4])
set(gca,'XTick',linspace(0,360,25))


end


figure(1)
xlabel('Duration (hours)')
ylabel('Hs (m)')
legend(path1,path2,path3,path4,path5,path6)
legend(tag1,tag2,tag3,tag4,tag5,tag6)
axis([0 48 0 2.8])
set(gca,'XTick',linspace(0,48,9));


figure(2)
xlabel('Duration (hours)')
ylabel('mss_{f< 1 Hz}')
legend(path1,path2,path3,path4,path5,path6)
legend(tag1,tag2,tag3,tag4,tag5,tag6)
axis([0 48 0 0.025])
set(gca,'XTick',linspace(0,48,9));

figure(3)
xlabel('Duration (hours)')
ylabel('C_D')
legend(path1,path2,path3,path4,path5,path6)
legend(tag1,tag2,tag3,tag4,tag5,tag6)

figure(4)
xlabel('Duration (hours)')
ylabel('mss2/mss1')
legend(path1,path2,path3,path4,path5,path6)
legend(tag1,tag2,tag3,tag4,tag5,tag6)
axis([0 48 0 1])
set(gca,'XTick',linspace(0,48,9));

figure(5)
xlabel('frequency (Hz)')
ylabel('E(f) (m^2/Hz)')
legend(path1,path2,path3,path4,path5,path6)
loglog((freq),(freq.^-4)./1000,'k-','LineWidth',1)
loglog((freq),(freq.^-5)./1000,'k-','LineWidth',1)
axis([0.04 2 0.001 10])

figure(6)
xlabel('frequency (Hz)')
ylabel('E(f) \times f^5 (m^2 )')
legend(tag1,tag2,tag3,tag4,tag5,tag6)

figure(7)
xlabel('frequency (Hz)')
ylabel('E(f) \times f^4 (m^2 )')
legend(tag1,tag2,tag3,tag4,tag5,tag6)


figure(8)
xlabel('Duration (hours)')
ylabel('Age')
legend(path1,path2,path3,path4,path5,path6)


figure(13)
xlabel('frequency (Hz)')
ylabel('Opposition integral I(f)')
legend(tag1,tag2,tag3,tag4,tag5,tag6)


