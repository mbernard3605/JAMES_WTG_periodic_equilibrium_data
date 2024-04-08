%Script to plot the simple metrics from wtg simulation:
%Currently plotting Vertical Velocity, Vertical velocity scaled by
%precipitation, and calculates the gms, II, SF, surface Evap
%current model runs are supposed to be 75 days, we are averaging the last
%20 days of this for statistics. In case of periodic this may or may not
%match the period of the oscillation 
function [outstats]=simple_wtg_plot_tropolimit(directory,wrfins,plotDir,cl,plt)

if isempty(wrfins)
    outstats=0;
    return;
end


if cl==1
    close all;
end
load('ERA5_EOFS.mat');
weights_old=weights;



[strt,nd]=regexp(plotDir,'t/\w*/D');
expName=plotDir(strt+2:nd-2);


R=287.04; kap=2/7; p0=1000e2; cp=R/kap;  G=9.81;      sclht = R*256./G;
Lv = 2.5e6;
Rd = 287.06;
Cp = 1005; 
Lv = 2.5e6;
Lf = 3.33e5;
lvlSize = 60;
g = 9.80665;
Tr = 273.15;
Pr = 1e5;

started=0;

a1=zeros(121,15,length(wrfins));a2=zeros(121,6,length(wrfins));
%set(gca,'ylim',[10000 100000],'ydir','reverse');

stats = zeros(length(wrfins),12);
temperature = {};
moisture = {};
omega = {};
gms = {};
gms_scm = {};
gms_crit = {};
gms_crit_scm = {};
time={};
lh={};
sh={};



for i = 1:length(wrfins)
   if regexp(wrfins{i},'\w*base\w*')
      qvtemp=nc_varget([directory wrfins{i}],'QV_LARGESCALE');
      ztmp = ncread([directory wrfins{i}],'PH')+ncread([directory wrfins{i}],'PHB')/9.81;
      ztmp=(mean(ztmp(2:end,:)+ztmp(1:end-1,:),2))/2;
      zls=nc_varget([directory wrfins{i}],'Z_FORCE');
      qv0=interp1(zls(2,:),qvtemp(2,:),ztmp,[],'extrap');

      break;
      
   else
      qvtemp=nc_varget([directory wrfins{i}],'QV_LARGESCALE');
      if size(qvtemp,1)>size(qvtemp,2)
      ztmp = ncread([directory wrfins{i}],'PH')+ncread([directory wrfins{i}],'PHB')/9.81;
      ztmp=(mean(ztmp(2:end,:)+ztmp(1:end-1,:),2))/2;
      zls=nc_varget([directory wrfins{i}],'Z_FORCE');
      qv0=interp1(zls(2,:),qvtemp(2,:),ztmp,[],'extrap');
      else
      ztmp = ncread([directory wrfins{i}],'PH')+ncread([directory wrfins{i}],'PHB')/9.81;
      ztmp=(mean(ztmp(:,2:end)+ztmp(:,1:end-1),1))/2;
      zls=nc_varget([directory wrfins{i}],'Z_FORCE');
      qv0=interp1(zls(:,2),qvtemp(:,2),ztmp,[],'extrap')';
      end
   end
    
    
end



for fileIndex=1:length(wrfins)
    
wrf=[directory wrfins{fileIndex}];
    [strt,nd]=regexp(wrf,'E_\w*_avg');
z = ncread(wrf,'PH')+ncread(wrf,'PHB')/9.81;
if size(z,1)>size(z,2)
    

zfull = (z(:,2:end)+z(:,1:end-1))/2;
zbig = mean(z);

rain=ncread(wrf,'RAINNC');
rr=diff(rain)*4;

time = (1:size(z,1))/4;

z_ongrid = (mean(z(2:end,:)+z(1:end-1,:),2))/2;

colors = jet(length(time));
[x,y] = meshgrid(time,z_ongrid);
[xshort,yshort] = meshgrid(time,z_ongrid(1:end-1));
P = ncread(wrf,'P')+ncread(wrf,'PB');

TH = ncread(wrf,'T')+300;
t = TH.*(P./Pr).^(2/7);


HFX = ncread(wrf,'HFX');
LH = ncread(wrf,'LH');
OMEGA_WTG = ncread(wrf,'OMEGA_WTG')';
qf = ncread(wrf,'QICE')+ncread(wrf,'QSNOW')+ncread(wrf,'QGRAUP');
ql = ncread(wrf,'QRAIN')+ncread(wrf,'QCLOUD');
qv = ncread(wrf,'QVAPOR');

if plt
figure(1),hold on;
omega{fileIndex}=OMEGA_WTG;
p(:,fileIndex)=mean(P(end-120:end,:));
z_save{fileIndex}=z;
h_omega(fileIndex)=plot(mean(OMEGA_WTG(:,end-120:end),2),mean(P(end-120:end,:),1),'linewidth',2);
set(gca,'ylim',[10000 100000],'ydir','reverse');
[sname1,ename1]=regexp(wrf,'wrfout_\w*_DE');
name = [wrf(sname+7:ename-3) wrf(strt+1:nd-4)];
spaces=regexp(name,'_');
name(spaces)=' ';
names{fileIndex} = name;

figure(2),hold on;
moisture{fileIndex}=qv;
sat_moisture{fileIndex}=qvs;
h_moisture(fileIndex) = plot((mean(qv(end-120:end,:)-qv0'))*1000,mean(P(end-12:end,:)),'linewidth',2);
set(gca,'ylim',[10000 100000],'ydir','reverse');

figure(3),hold on;
temperature{fileIndex}=t;
h_temperature(fileIndex) = plot((mean(t(end-120:end,:)-t(1,:))),mean(P(end-12:end,:)),'linewidth',2);
set(gca,'ylim',[10000 100000],'ydir','reverse');
end

S = Cp*t+zfull*9.81;
T = Cp*t;
q = Lv*qv;
L = Lv*ql;
F = Lf*qf;
H = S+q-F;

[~,~,Hstar,Hp,Tp] = zero_plume_buoyancy_frozen(t,qv,qf,zfull,P,4000);

dp = diff(P,1,2);

wtg_th=ncread(wrf,'SCMTHTEN');
wtg_t=wtg_th.*(P./Pr).^(2/7);
wtg_qv=ncread(wrf,'SCMQVTEN');

dH_p=diff(H,1,2);
dS_p=diff(S,1,2);
dT_p=diff(T,1,2);
dQ_p=diff(q,1,2);
dF_p=diff(F,1,2);
dL_p=diff(L,1,2);

dHdp = dH_p./dp;
dSdp = dS_p./dp;
dTdp = dT_p./dp;
dQdp = dQ_p./dp;
dFdp = dF_p./dp;
dLdp = dL_p./dp;

dHdt=diff(H,1,1)./(6*3600);
dSdt=diff(S,1,1)./(6*3600);
dTdt=diff(T,1,1)./(6*3600);
dQdt=diff(q,1,1)./(6*3600);
dFdt=diff(F,1,1)./(6*3600);
dLdt=diff(L,1,1)./(6*3600);


QH_p = zeros(size(dHdp,1),1);QS_p=QH_p;QT_p=QH_p;QQ_p=QH_p;QL_p=QH_p;QF_p=QH_p;
QH_tWTG=QH_p;QH_t=QH_p(2:end);QH_plocal=QH_p;QS_t=QH_t;QS_t=QH_t;QQ_t=QH_t;QF_t=QH_t;QL_t=QH_t;
QH_qvWTG=QH_p;QQVDIF=QH_p;QTDIF=QH_p;
pmid=(P(:,2:end)+P(:,1:end-1))/2;
omegamid=(OMEGA_WTG(2:end,:)+OMEGA_WTG(1:end-1,:))/2;

ptop=find(pmid(2,:)<=10000);

for i = 1:size(QH_p,1)-1
    QS_p(i)=trapz(pmid(i,1:ptop(end)),-dSdp(i,1:ptop(end)).*omegamid(1:ptop(end),i)'/9.81);
    QT_p(i)=trapz(pmid(i,1:ptop(end)),-dTdp(i,1:ptop(end)).*omegamid(1:ptop(end),i)'/9.81);
    QH_p(i)=trapz(pmid(i,1:ptop(end)),-dHdp(i,1:ptop(end)).*omegamid(1:ptop(end),i)'/9.81);
    QQ_p(i)=trapz(pmid(i,1:ptop(end)),-dQdp(i,1:ptop(end)).*omegamid(1:ptop(end),i)'/9.81);
    QF_p(i)=trapz(pmid(i,1:ptop(end)),-dFdp(i,1:ptop(end)).*omegamid(1:ptop(end),i)'/9.81);
    QL_p(i)=trapz(pmid(i,1:ptop(end)),-dLdp(i,1:ptop(end)).*omegamid(1:ptop(end),i)'/9.81);
    QH_tWTG(i)=trapz(P(i,1:ptop(end)),-wtg_t(i,1:ptop(end))/9.81)*Cp;
    QH_qvWTG(i)=trapz(P(i,1:ptop(end)),-wtg_qv(i,1:ptop(end))/9.81)*Lv;
    QH_t(i)=trapz(P(i,1:ptop(end)),-dHdt(i,1:ptop(end))/9.81);
    QS_t(i)=trapz(P(i,1:ptop(end)),-dSdt(i,1:ptop(end))/9.81);
    QT_t(i)=trapz(P(i,1:ptop(end)),-dTdt(i,1:ptop(end))/9.81);
    QQ_t(i)=trapz(P(i,1:ptop(end)),-dQdt(i,1:ptop(end))/9.81);
    QF_t(i)=trapz(P(i,1:ptop(end)),-dFdt(i,1:ptop(end))/9.81);
    QL_t(i)=trapz(P(i,1:ptop(end)),-dLdt(i,1:ptop(end))/9.81);
   % QH_plocal(i)=trapz((P(i,2:end)+P(i,1:end-1))/2,-wdhdp(i,1:end)/9.81);
end


else
   
        

zfull = (z(2:end,:)+z(1:end-1,:))/2;
zbig = mean(z,2);


time = (1:size(z,2))/4;

z_ongrid = (mean(z(2:end,:)+z(1:end-1,:),2))/2;

colors = jet(length(time));
[x,y] = meshgrid(time,z_ongrid);
[xshort,yshort] = meshgrid(time,z_ongrid(1:end-1));

rain=ncread(wrf,'RAINNC');
rr=diff(rain)*4;

P = ncread(wrf,'P')+ncread(wrf,'PB');

TH = ncread(wrf,'T')+300;
t = TH.*(P./Pr).^(2/7);
tref = ncread(wrf,'TH_LARGESCALE');
zref = ncread(wrf,'Z_FORCE');

OMEGA_WTG = ncread(wrf,'OMEGA_WTG')';






[pcs,lds_ERA_interp]=Calc_TH_angle(OMEGA_WTG,P');
eof1=-(lds_ERA_interp(1:end-1,1)+lds_ERA_interp(2:end,1))/2;
eof2=(lds_ERA_interp(1:end-1,2)+lds_ERA_interp(2:end,2))/2;
wres=lds_ERA_interp(:,3:end)*pcs(3:end,:);
wresmid=(wres(2:end,:)+wres(1:end-1,:))/2;
o1=-pcs(1,:);
o2=pcs(2,:);
o3=pcs(3,:);
w1=eof1*o1;
w2=eof2*o2;
angle_ERA=atan2d(o2,o1);
angle_mean_ERA=atan2d(mean(o2(end-120:end)),mean(o1(end-120:end)));

HFX = ncread(wrf,'HFX');
LH = ncread(wrf,'LH');
qv = ncread(wrf,'QVAPOR');
qf = ncread(wrf,'QICE')+ncread(wrf,'QSNOW')+ncread(wrf,'QGRAUP');
ql = ncread(wrf,'QRAIN')+ncread(wrf,'QCLOUD');

S = Cp*t+zfull*9.81;
T = Cp*t;%+zfull*9.81;
q = Lv*qv;
L = Lv*ql;
F = Lf*qf;
H = S+q-F;
pmean=mean(P(:,end-120:end),2);
qvs=calc_qvstar(t,pmean);

[~,~,Hstar,Hp,Tp] = zero_plume_buoyancy_frozen(nanmean(t(:,end-121:end),2),nanmean(qv(:,end-121:end),2),nanmean(qf(:,end-121:end),2),nanmean(zfull(:,end-121:end),2),nanmean(P(:,end-121:end),2),4000);

SF=trapz(pmean,mean(qv(:,end-120:end),2))./trapz(pmean,mean(qvs(:,end-120:end),2));


[sname1,ename1]=regexp(wrf,'wrfout_\w*_DE');
name = [wrf(sname1+7:ename1-3) wrf(strt+1:nd-4)];
spaces=regexp(name,'_');
name(spaces)=' ';
names{fileIndex} = name(11:end);
temperature{fileIndex}=t;

p(:,fileIndex)=mean(P(:,end-120:end),2);
omega{fileIndex}=OMEGA_WTG;
z_save{fileIndex}=z;
moisture{fileIndex} = qv;
sat_moisture{fileIndex}=qvs;


dp = diff(P,1,1);


wtg_th=ncread(wrf,'SCMTHTEN');
wtg_t=wtg_th.*(P./Pr).^(2/7);
wtg_qv=ncread(wrf,'SCMQVTEN');

dH_p=diff(H,1,1);
dS_p=diff(S,1,1);
dT_p=diff(T,1,1);
dTH_p=diff(TH*cp,1,1);
dQ_p=diff(q,1,1);
dF_p=diff(F,1,1);
dL_p=diff(L,1,1);

dHdp = dH_p./dp;
dSdp = dS_p./dp;
dTdp = dT_p./dp;
dTHdp = dTH_p./dp;
dQdp = dQ_p./dp;
dFdp = dF_p./dp;
dLdp = dL_p./dp;

dHdt=diff(H,1,2)./(6*3600);
dSdt=diff(S,1,2)./(6*3600);
dTdt=diff(T,1,2)./(6*3600);
dQdt=diff(q,1,2)./(6*3600);
dFdt=diff(F,1,2)./(6*3600);
dLdt=diff(L,1,2)./(6*3600);





QH_p = zeros(size(dHdp,2),1);QS_p=QH_p;QT_p=QH_p;QQ_p=QH_p;QL_p=QH_p;QF_p=QH_p;
QH_tWTG=QH_p;QTH_p=QH_p;QH_t=QH_p(2:end);QH_plocal=QH_p;QS_t=QH_t;QT_t=QH_t;QQ_t=QH_t;QF_t=QH_t;QL_t=QH_t;
QH_qvWTG=QH_p;QQVDIF=QH_p;QTDIF=QH_p;
QH1_p=QH_p;QH2_p=QH_p ; QS1_p=QH_p ;QS2_p=QH_p ;
pmid=(P(2:end,:)+P(1:end-1,:))/2;
omegamid=(OMEGA_WTG(:,2:end)+OMEGA_WTG(:,1:end-1))/2;

ptop=find(pmid(:,2)>=10000);
ptop=length(pmid(:,2));

for i = 1:size(QH_p,1)
    QS_p(i)=trapz(pmid(1:ptop(end),i),-dSdp(1:ptop(end),i).*omegamid(i,1:ptop(end))'/9.81);
    QS1_p(i)=trapz(pmid(1:ptop(end),i),-dSdp(1:ptop(end),i).*w1(1:ptop(end),i)/9.81);
    QS2_p(i)=trapz(pmid(1:ptop(end),i),-dSdp(1:ptop(end),i).*w2(1:ptop(end),i)/9.81);
    QSres_p(i)=trapz(pmid(1:ptop(end),i),-dSdp(1:ptop(end),i).*wresmid(1:ptop(end),i)/9.81);
    QT_p(i)=trapz(pmid(1:ptop(end),i),-dTdp(1:ptop(end),i).*omegamid(i,1:ptop(end))'/9.81);
    QTH_p(i)=trapz(pmid(1:ptop(end),i),-dTHdp(1:ptop(end),i).*omegamid(i,1:ptop(end))'/9.81);
    QH_p(i)=trapz(pmid(1:ptop(end),i),-dHdp(1:ptop(end),i).*omegamid(i,1:ptop(end))'/9.81);
    QH1_p(i)=trapz(pmid(1:ptop(end),i),-dHdp(1:ptop(end),i).*w1(1:ptop(end),i)/9.81);
    QH2_p(i)=trapz(pmid(1:ptop(end),i),-dHdp(1:ptop(end),i).*w2(1:ptop(end),i)/9.81);
    QHres_p(i)=trapz(pmid(1:ptop(end),i),-dHdp(1:ptop(end),i).*wresmid(1:ptop(end),i)/9.81);
    QQ_p(i)=trapz(pmid(1:ptop(end),i),-dQdp(1:ptop(end),i).*omegamid(i,1:ptop(end))'/9.81);
    QF_p(i)=trapz(pmid(1:ptop(end),i),-dFdp(1:ptop(end),i).*omegamid(i,1:ptop(end))'/9.81);
    QL_p(i)=trapz(pmid(1:ptop(end),i),-dLdp(1:ptop(end),i).*omegamid(i,1:ptop(end))'/9.81);
    QH_tWTG(i)=trapz(P(1:ptop(end),i),wtg_t(1:ptop(end),i)/9.81)*Cp;
    QH_qvWTG(i)=trapz(P(1:ptop(end),i),wtg_qv(1:ptop(end),i)/9.81)*Lv;
end

if plt
figure(1),hold on;
h_omega(fileIndex) = plot(mean(dH_p(:,end-120:end).*w1(:,end-120:end),2),mean(pmid(:,end-12:end),2),'linewidth',2);
set(gca,'ylim',[10000 100000],'ydir','reverse','fontsize',16);
pause(.1);
figure(2),hold on;
h_temperature(fileIndex) = plot(mean(H(:,end-120:end),2),mean(P(:,end-12:end),2),'linewidth',2);
set(gca,'ylim',[10000 100000],'ydir','reverse','fontsize',16);
pause(.1);
figure(3),hold on;
h_temperature(fileIndex) = plot(mean(dH_p(:,end-120:end),2),mean(pmid(:,end-12:end),2),'linewidth',2);
set(gca,'ylim',[10000 100000],'ydir','reverse','fontsize',16);
pause(.1);
end
    
        for i = 1:size(QH_t,1)
            QH_t(i)=trapz(P(1:ptop(end),i),-dHdt(1:ptop(end),i)/9.81);
            QS_t(i)=trapz(P(1:ptop(end),i),-dSdt(1:ptop(end),i)/9.81);
            QT_t(i)=trapz(P(1:ptop(end),i),-dTdt(1:ptop(end),i)/9.81);
            QQ_t(i)=trapz(P(1:ptop(end),i),-dQdt(1:ptop(end),i)/9.81);
            QF_t(i)=trapz(P(1:ptop(end),i),-dFdt(1:ptop(end),i)/9.81);
            QL_t(i)=trapz(P(1:ptop(end),i),-dLdt(1:ptop(end),i)/9.81);
          % QH_plocal(i)=trapz((P(i,2:end)+P(i,1:end-1))/2,-wdhdp(i,1:end)/9.81);
        end

    
    
    end

    qvtend=ncread(wrf,'QV_LARGESCALE_TEND')*Lv;


    
            zforce=ncread(wrf,'Z_FORCE')';
            qraten=ncread(wrf,'RTHRATEN');

            if size(zforce,1)>size(zforce,2)   
                try
            zt=find(zforce(2,:)>=17000);
            zt=zt(1)-1;
                catch
                    
            zt=find(zforce(2,:)>=15000);
            zt=zt(1)-1;
                end
            else
            zt=find(zforce(:,2)>=17000);
            zt=zt(1)-1;               
            end
    
    
    if find(qraten~=0)
       
        prs=nc_varget(wrf,'P')+nc_varget(wrf,'PB');
        pii = (prs'/1e5).^(kap); 
        mu=nc_varget(wrf,'MU')'+nc_varget(wrf,'MUB')';
        zsize = nc_varsize(wrf,'ZNW');
        znw = nc_varget(wrf,'ZNW',[0 0],[1 zsize(2)])';
        znu = nc_varget(wrf,'ZNU',[0 0],[1 zsize(2)-1]);
        nt=numel(mu);
        for iz=1:(zsize(2)-1)
            qraten(iz,:) = qraten(iz,:)./mu;
            
        end
        qra_int=zeros(1,nt);
        for i=1:numel(znu)
         qra_int = qra_int + qraten(i,:).*pii(i,:).*(-znw(i+1)+znw(i)).*mu/G;
        end

       qrad_heating=qraten;
        QRAD=qra_int*Cp;

%clear QRAD;

   else
        QRAD=zeros(size(QH_tWTG));
        SW=zeros(length(time),1);
        LW=SW;
        qrad=ncread(wrf,'QRADIATION_SCM');
        qrad_heating=zeros(size(qrad));
        if size(qrad,1) < size(qrad,2)
            zforce=zforce';
            pforce=zeros(size(zforce));
            qraten=zeros(size(P));
            for index = 1:length(time)
                
                qraten(:,index)=interp1(zforce(:,index),qrad(:,index),zfull(:,index),[],'extrap');
                pforce(:,index)=interp1(zfull(:,index),P(:,index),zforce(:,index),[],'extrap');
                QRAD(index)=trapz(pforce(:,index),-qrad(:,index).*(pforce(:,index)./Pr).^(2/7)/9.81)*Cp;
                qrad_heating(:,index)=(qrad(:,index).*(pforce(:,index)./Pr).^(2/7));
            end
        else
            qrad=qrad';
            zfull=zfull';

            pforce=zeros(size(zforce));
            P=P';
            qraten=zeros(size(P));
            for index = 1:length(time)
                qraten(:,index)=interp1(zforce(:,index),qrad(:,index),zfull(:,index),[],'extrap');
                pforce(:,index)=interp1(zfull(:,index),P(:,index),zforce(:,index),[],'extrap');
                QRAD(index)=trapz(pforce(:,index),-qrad(:,index).*(pforce(:,index)./Pr).^(2/7)/9.81)*Cp;
                qrad_heating(:,index)=(qrad(:,index).*(pforce(:,index)./Pr).^(2/7));
            end
    
        end
    end

if size(QRAD,1) ~= size(QH_tWTG,1)
QRAD=QRAD';
end

        if size(zforce,1) < size(zforce,2)
            
            size(zforce);
            size(qvtend);
            size(zfull);
            
            pforce=zeros(size(zforce));
            qvten=zeros(size(P));
            for index = 1:length(time)
                
                qvten(:,index)=interp1(zforce(:,index),qvtend(:,index),zfull(:,index),[],'extrap');
                pforce(:,index)=interp1(zfull(:,index),P(:,index),zforce(:,index),[],'extrap');
                Qhadv(index)=trapz(pforce(:,index),-qvtend(:,index)/9.81);
        
            end
        else
            zforce=zforce';
            pforce=zeros(size(zforce));

            qvten=zeros(size(P));
            for index = 1:length(time)
                qvten(:,index)=interp1(zforce(:,index),qvtend(:,index),zfull(:,index),[],'extrap');
                pforce(:,index)=interp1(zfull(:,index),P(:,index),zforce(:,index),[],'extrap');
                Qhadv(index)=trapz(pforce(:,index),-qvtend(:,index)/9.81);
        
            end
    
        end

        dt=diff(mean(TH(:,end-119:end),2));
        dtdp=dt./mean(dp(:,end-119:end),2);
        dtdp(1:10)=mean(dtdp(11:15));
        try
        qrad_heating_interp=interp1(pforce(:,2),qrad_heating(:,2),nanmean(P(:,end-120:end),2),[],'extrap');
        omega_rad=midData(qrad_heating_interp)./dtdp;
        omega_rad(mean(pmid(:,end-120:end),2)<10000)=0;
        omega_rad=interp1(mean(pmid(:,end-120:end),2),omega_rad,mean(P(:,end-120:end),2),[],'extrap');
        catch
            omega_rad=zeros(size(P));
        end 
        

p_save{fileIndex}=pforce;
disp(name);
%disp(mean(-QH_tWTG(end-120:end))-(mean(rr(end-120:end))*29+mean(QRAD(end-120:end))+mean(HFX(end-120:end))));
%disp(mean(QTH_p(end-120:end))-(mean(rr(end-120:end))*29+mean(QRAD(end-120:end))+mean(HFX(end-120:end))));
%disp(mean(QS_p(end-120:end))-(mean(rr(end-120:end))*29+mean(QRAD(end-120:end))+mean(HFX(end-120:end))));
correct=-QS_p(end-120:end)\QH_tWTG(end-120:end);
disp(correct);
QS_p=QS_p*(correct);
QS1_p=QS1_p*(correct);
QS2_p=QS2_p*(correct);
QSres_p=QSres_p*correct;
hadv{fileIndex}=QH_qvWTG+QQ_p;
if exist('qrad_heating','var')
rad1{fileIndex}=qrad_heating;
end
rad2{fileIndex}=QRAD;
lh{fileIndex}=LH;
sh{fileIndex}=HFX;
rainrate{fileIndex}=rr*2.5e6/86400;
advectDSE{fileIndex}=QS_p;
advectMSE{fileIndex}=QS_p+QQ_p;
advectDSE1{fileIndex}=QS1_p;
advectMSE1{fileIndex}=QH1_p+QS1_p;
advectDSE2{fileIndex}=QS2_p;
advectMSE2{fileIndex}=QH2_p+QS2_p;
advectMSEres{fileIndex}=QHres_p+QSres_p;
advectDSEres{fileIndex}=QSres_p;
angle{fileIndex}=angle_ERA;
O1{fileIndex}=o1;
O2{fileIndex}=o2;
O3{fileIndex}=o3;
h{fileIndex}=H;
hsat{fileIndex}=Hstar;
hp{fileIndex}=Hp;
tp{fileIndex}=Tp;
angle_mean(fileIndex)=angle_mean_ERA;
advectDSE_SCM{fileIndex}=QH_tWTG;
advectMSE_SCM{fileIndex}=QH_qvWTG+QH_tWTG;
sf(fileIndex)=SF;
omega_rad_save{fileIndex}=omega_rad;


%a1(:,:,fileIndex)=[QS_p(end-120:end) QQ_p(end-120:end) QH_p(end-120:end)...
%    QH_tWTG(end-120:end)+QH_qvWTG(end-120:end) LH(end-120:end)+HFX(end-120:end) QRAD(end-120:end)...
%    (QH_tWTG(end-120:end)+QH_qvWTG(end-120:end)+LH(end-120:end)+HFX(end-120:end)+QRAD(end-120:end)-QH_t(end-120:end))...
%    (-QH_p(end-120:end)+LH(end-120:end)+HFX(end-120:end)+QRAD(end-120:end)-QH_t(end-120:end))...
%    QH_t(end-120:end) QS_t(end-120:end) QL_t(end-120:end) QF_t(end-120:end) QQ_t(end-120:end) QH_t(end-120:end) QH_t(end-120:end)];
%a1mean=squeeze(mean(a1));

 a2(:,:,fileIndex)=[QH_tWTG(end-120:end)+QH_qvWTG(end-120:end) LH(end-120:end)+HFX(end-120:end) QRAD(end-120:end)...
    (QH_tWTG(end-120:end)+QH_qvWTG(end-120:end)+LH(end-120:end)+HFX(end-120:end)+QRAD(end-120:end)-QH_t(end-120:end))...
    QH_t(end-120:end) rr(end-120:end)*Lv/86400];
%a2mean=squeeze(mean(a2));
if mean2(QQ_p)==0
   QQ_p=-QH_qvWTG; 
end
omega_int=intOmega(wrf);
GMS_SCM=1+QH_qvWTG./QH_tWTG;
GMS_SCM_avg=nanmean(GMS_SCM(end-80:end));
GMS_vert=(QH_tWTG-QQ_p)./QH_tWTG;
GMS_crit_SCM = -(LH+HFX+QRAD)./QH_tWTG;
dryingEfficiency_SCM = GMS_SCM-GMS_crit_SCM;
GMS_crit_SCM_avg=nanmean(GMS_crit_SCM(end-80:end));
dryingEfficiency_SCM_end_avg=nanmean(dryingEfficiency_SCM(end-80:end));



gms{fileIndex}=GMS_vert;
gms_scm{fileIndex}=GMS_SCM;
gms_crit_scm{fileIndex}=GMS_crit_SCM;

Rain=mean(rr(end-120:end));
try
omega_500_std=std(OMEGA_WTG(end-120:end,27));
catch
omega_500_std=std(OMEGA_WTG(27,end-120:end)); 
end
omega_std=std(omega_int(end-120:end));
rain_std=std(rr(end-120:end));
stats(fileIndex,:)=[Rain GMS_SCM_avg nanmean(GMS_vert(end-120:end)) GMS_crit_SCM_avg ...
    dryingEfficiency_SCM_end_avg nanmean(-QH_tWTG(end-120:end)+QQ_p(end-120:end))...
    mean(-QH_tWTG(end-120:end)) nanmean(QH_qvWTG(end-120:end)+QQ_p(end-120:end)) mean(QRAD(end-120:end)) mean(LH(end-120:end)) mean(HFX(end-120:end)) angle_mean_ERA];
    
end
if plt
marker={'o' '^' 's' 'd' 'p' 'h'};

legend(h_omega,names,'Interpreter','none');
set(gca,'ylim',[10000 100000],'ydir','reverse');
legend(names,'Interpreter','none');
legend(h_temperature,names,'Interpreter','none');


figure('units','normalized','outerposition',[0 0 1 1]), hold on;
for varIndex = 1:size(a2,2)-1
    scatter(1:size(a2,3),mean(a2(end-120:end,varIndex,:)),200,'filled',marker{varIndex})
    
end
[H,icons,plots,legend_text] = legend({'Advection','Surface','Radiation','residual','dhdt'},'fontsize',16,'Interpreter','none');
set(gca,'xlim',[0 size(a2,3)+1],'xtick',[1:length(names)],'xticklabel',names,'fontsize',16)
title(['MSE budget components for ' expName]);
for i =length(icons)/2+1:length(icons)
icons(i).Children.MarkerSize=15;
end
grid on;
xtickangle(-35);
%print(gcf,'-djpeg','-r300',[plotDir 'MSE_budget_chart.jpg']);


%figure(1),title([expName 'Vertical Velocity'])
%print(gcf,'-djpeg','-r300',[plotDir 'Vertical_Velocity_compare.jpg']);
%figure(2),title([expName 'Moisture anomaly'])
%print(gcf,'-djpeg','-r300',[plotDir 'Moisture_anomaly_compare.jpg']);
%figure(3),title([expName 'Temperature anomaly'])
%print(gcf,'-djpeg','-r300',[plotDir 'Temperature_anomaly_compare.jpg']);
end



outstats.stats=stats;
outstats.omega=omega;
outstats.qv=moisture;
outstats.qvs=sat_moisture;
outstats.th=temperature;    
outstats.names=names;
outstats.p = p;
outstats.gms = gms;
outstats.gms_scm = gms_scm;
outstats.gms_crit = gms_crit;
outstats.gms_crit_scm = gms_crit_scm;
outstats.time = time;
outstats.lh = lh;
outstats.sh = sh;
if exist('rad1','var')
outstats.qrad = rad1;
end
outstats.QRAD = rad2;
outstats.rain = rainrate;
outstats.advectDSE = advectDSE;
outstats.advectMSE = advectMSE;
outstats.advectDSE1 = advectDSE1;
outstats.advectMSE1 = advectMSE1;
outstats.advectDSE2 = advectDSE2;
outstats.advectMSE2 = advectMSE2;
outstats.advectDSEres = advectDSEres;
outstats.advectMSEres = advectMSEres;
outstats.advectDSE_SCM = advectDSE_SCM;
outstats.advectMSE_SCM = advectMSE_SCM;
outstats.h=h;
outstats.hstar=hsat;
outstats.hp=hp;
outstats.tp=tp;
outstats.z = z_save;
outstats.Qhadv = hadv;
outstats.O1 = O1;
outstats.O2 = O2;
outstats.O3 = O3;
outstats.meanangle = angle_mean;
outstats.sf=sf;
outstats.pforce=p_save;
outstats.omega_rad=omega_rad_save;

end
%figure,








