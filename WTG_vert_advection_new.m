%Script that calculates the vertical advection from WTG simulations
clear variables;
close all;

load('anglecolormap5.mat');
angles=linspace(-180,180,60);


%Load Constants
Rd = 287.06;
Cp = 1005; 
Lv = 2.5e6;
Lf = 3.33e5;
lvlSize = 60;
g = 9.80665;
Tr = 273.15;
Pr = 1e5;

%Local File structur
dir = '';
fname = 'wrfout_east_85RH_avg';

wrf=[dir fname];


wrftime=[];
    timstr = nc_varget(wrf,'Times');
    wrftime_pre = wrftime; 
    wrftime = datenum(timstr);


R=287.04; kap=2/7; p0=1000e2; cp=R/kap;  G=9.81;      sclht = R*256./G;
Lv = 2.5e6;

    tfactor = wrftime(2)-wrftime(1); % tfactor
 tfactor 

    ztime = datenum(nc_varget(wrf,'Times'));
    tmsize = nc_varsize(wrf,'Times');
    rsize = nc_varsize(wrf,'RAINNC');
    tsize = nc_varsize(wrf,'T');
    phsize = nc_varsize(wrf,'PH');
    zsize = nc_varsize(wrf,'ZNW');
    
    znw = nc_varget(wrf,'ZNW',[0 0],[1 zsize(2)]);
    znu = nc_varget(wrf,'ZNU',[0 0],[1 zsize(2)-1]);
    znfac=(znw(1:end-1)-znu)./(-diff(znw));
    
    nt=tmsize(1);
t=(1:nt)/4;
    
    %th_wrf_at_zforce = zeros(size(zforce));
    
    ght_all = zeros(tsize(1:2));
    ghts_all = zeros(tsize(1),tsize(2)+1);
    
    for itc=0:min(nt,numel(ztime))-1
        it=itc+1;
        if(ztime(it)~=wrftime(it))
            continue;
        end
        
        %th1 = nc_varget(wrf1,'T',[it 0 0 0],[1 tsize(2:end) ]); th1=th1+300;
       if(length(phsize)==4)
        ph=nc_varget(wrf,'PH',[itc 0 0 0],[1 phsize(2:end) ]) + nc_varget(wrf,'PHB',[itc 0 0 0],[1 phsize(2:end) ]);
       else
        ph=nc_varget(wrf,'PH',[itc 0 ],[1 phsize(2:end) ]) + nc_varget(wrf,'PHB',[itc 0 ],[1 phsize(2:end) ]);
       end
        ghts = exp(-ph/(G*sclht));

        if( (length(phsize)==4 && phsize(3)*phsize(4)==1) | length(phsize)==2 ); ghts=ghts'; end

        ght=zeros(1,tsize(2:end));
        for k=1:length(znfac)
            ght(k)=-sclht*log(ghts(k)*znfac(k)+ghts(k+1)*(1-znfac(k)));
        end
        zwrf_avg = ght;
        %th_avg = mean(mean(mean(th1,4),3),2);
        
        %qraten=nc_varget(wrf1,'RTHRATEN');
        %[tmk_intp]=interp1(zwrf_avg',th_avg, zforce(it,:),'linear','extrap');
        
        %th_wrf_at_zforce(it,:) = tmk_intp;
        ght_all(it,:) = zwrf_avg';
        ghts_all(it,:)=ph/G;
        %plot(tmk_intp,zforce(it,:),'.'); hold on; plot(th_avg,zwrf_avg','r.');
        %clear ght ghts ph
        
    end
    zmean=mean(ght_all);
    P=nc_varget(wrf,'P')'+nc_varget(wrf,'PB')';
    pii = (P/1e5).^(kap);
OMEGA_WTG = ncread(wrf,'OMEGA_WTG');


load('ERA5_EOFS.mat');
weights_old=weights;
lds_ERA=lds(:,1:2);

    presmean=mean(P,2);
  
dp=diff(presmean);
weightsShapeFull=[dp; dp(end)]./sum([dp; dp(end)]);
%Eliminate negative weights
%weightsShapeFull(weightsShapeFull<=0)=-weightsShapeFull(weightsShapeFull<=0);  
    
omega=OMEGA_WTG';
omega2=(omega).*(weightsShapeFull'.^0.5);
    
    height_index=find(level>=presmean(end));
lds_ERA=lds_ERA(height_index,:);
lds_ERA(:,1)=-lds_ERA(:,1);
omega2_interp=interp1(presmean,omega2',level(height_index),[],'extrap');
pcs_ERA=lds_ERA'*omega2_interp;
lds_ERA_interp=interp1(level(height_index),lds_ERA,presmean,[],'extrap');
lds_ERA_plot=interp1(level(height_index),lds_ERA./(weights.^0.5),presmean,[],'extrap');
cmap=sqrt((presmean(49)-presmean(1))./trapz(presmean(1:49)',lds_ERA_plot(1:49,:).^2));
lds_scaled=lds_ERA_plot.*cmap;
pcs_scaled=pcs_ERA./cmap';
pc_test=pcs_scaled;
%lds_plot=lds_scaled.*(weightShapeFull.^0.5);


   o1_ERA=pcs_scaled(1,:);
   o2_ERA=pcs_scaled(2,:);
ang_ERA=atan2d(o2_ERA,o1_ERA);
ang_ERA_colors=interp1(angles,testmap,ang_ERA,[],'extrap');

      o1_ERA_mid=midData(pcs_scaled(1,:)')';
   o2_ERA_mid=midData(pcs_scaled(2,:)')';
ang_ERA_mid=atan2d(o2_ERA_mid,o1_ERA_mid);
ang_ERA_mid_colors=interp1(angles,testmap,ang_ERA_mid,[],'extrap'); 
    
    


SCMTHTEN = ncread(wrf,'SCMTHTEN').*pii*Cp; %
SCMQVTEN = ncread(wrf,'SCMQVTEN')*Lv;
z_force = nc_varget(wrf,'Z_FORCE');
qvhadv=nc_varget(wrf,'QV_LARGESCALE_TEND');
rr=diff(nc_varget(wrf,'RAINNC'));
rr=[0; rr]*4*Lv/86400;
qvhadv_int=zeros(length(t),size(ght_all,2));
for i = 1:length(t)
qvhadv_int(i,:) = interp1(z_force(2,:),qvhadv(i,:)',ght_all(i,:),[],'extrap');
end
qv = ncread(wrf,'QVAPOR');
TH = ncread(wrf,'T')+300;
T = TH.*(P./Pr).^(2/7);
qvs=calc_qvstar(T,presmean);
rho = P./T./R;
rho_mid=midData(rho);
www = -OMEGA_WTG./rho/G;
www_mid = midData(www);
omega_mid = midData(OMEGA_WTG);
S = Cp*T+ght_all'*g;
H = S + qv*Lv;
sdiff = diff(S,1,1);
qvdiff=diff(qv,1,1);
dp = diff(P,1,1);
dz = diff(ght_all,1,2);
dSlargedp = sdiff./dp;
dqvlargedp = qvdiff./dp;
dSlargedz = sdiff./dz';
dqvlargedz = qvdiff./dz';

Pmid=midData(P);

for i = 1:size(dSlargedp,2)
%total advection calculated by the model
QStot(i)=trapz(ght_all(i,:),rho(:,i).*SCMTHTEN(:,i));
Qqvtot(i)=trapz(ght_all(i,:),rho(:,i).*SCMQVTEN(:,i));
Qqvhorz(i) = trapz(ght_all(i,:),rho(:,i)'.*qvhadv_int(i,:))*Lv;
%vertical advection from the WTG vertical velocity    

QS(i) = trapz(Pmid(:,i),(-dSlargedp(:,i).*midData(OMEGA_WTG(:,i))/G));
Qqv(i) = trapz(Pmid(:,i),(-dqvlargedp(:,i).*midData(OMEGA_WTG(:,i))/G))*Lv;

QS1(i) = trapz(Pmid(:,i),(-dSlargedp(:,i).*midData(lds_scaled(:,1))/G));
Qqv1(i) = trapz(Pmid(:,i),(-dqvlargedp(:,i).*midData(lds_scaled(:,1))/G))*Lv;


QS2(i) = trapz(Pmid(:,i),(-dSlargedp(:,i).*midData(lds_scaled(:,2))/G));
Qqv2(i)=trapz(Pmid(:,i),(-dqvlargedp(:,i).*midData(lds_scaled(:,2))/G))*Lv;

QSpc(i,:) = trapz(Pmid(:,i),(-dSlargedp(:,i).*midData(lds_scaled)/G));
Qqvpc(i,:) = trapz(Pmid(:,i),(-dqvlargedp(:,i).*midData(lds_scaled)/G))*Lv;


end


QSpc_mean = trapz(midData(presmean),(-mean(dSlargedp,2).*midData(lds_scaled)/G));
Qqvpc_mean = trapz(midData(presmean),(-mean(dqvlargedp,2).*midData(lds_scaled)/G))*Lv;
QHpc_mean=QSpc_mean+Qqvpc_mean;

QH=QS+Qqv;
QHpc=QSpc+Qqvpc;
QH1=QS1+Qqv1;
QH2=QS2+Qqv2;
QHtot=QStot+Qqvtot;
QHhorz=Qqvhorz;QHhorz(1)=QHhorz(2);
QHvert=QHtot-QHhorz;
QHvert_res=QH-QHpc_mean(1).*pc_test(1,:)-QHpc_mean(2).*pc_test(2,:);
QHvert_res1=(QH1.*pc_test(1,:)+QH2.*pc_test(2,:))-QHpc_mean(1).*pc_test(1,:)-QHpc_mean(2).*pc_test(2,:);
Qqvvert=Qqvtot-Qqvhorz;

%     qralw=nc_varget(wrf,'RTHRATLW');
%     qrasw=nc_varget(wrf,'RTHRATSW');
    qraten=nc_varget(wrf,'RTHRATEN');
    
            ist=1;

    

    lh=nc_varget(wrf,'ACLHF')*4/86400;     
    Qlh=[0 ; diff(lh)];
    hfx=nc_varget(wrf,'ACHFX')*4/86400;     
    Qsh=[0; diff(hfx)];
    
    
mu=nc_varget(wrf,'MU')+nc_varget(wrf,'MUB');
    %whos mu
    %mun=reshape(mu,[481 1 31 31]);
    %mun=reshape(mu,[481 1 31 31]);
    msize = size(mu);
    mun = reshape(mu, [msize(1) 1 msize(2:end)]);
    for iz=1:(zsize(2)-1)
        qraten(:,iz) = qraten(:,iz)./mun;
    end
   

    ny=size(mu,2);
    nx=size(mu,3);
    nt = size(mu,1);
    qlw_int=zeros(nt,ny,nx);    qra_int=zeros(nt,ny,nx);    qsw_int=zeros(nt,ny,nx); qt_int=zeros(nt,ny,nx); qqv_int=zeros(nt,ny,nx);
    for i=1:numel(znu)
         qra_int = qra_int + squeeze(qraten(:,i)).*squeeze(pii(i,:))'.*(-znw(i+1)+znw(i)).*squeeze(mu)/G;
%          qlw_int(:,:,:) = qlw_int(:,:,:) + squeeze(qralw(:,i,:,:)).*squeeze(pii(:,i,:,:)).*(-znw(i+1)+znw(i)).*squeeze(mu(:,:,:))/G;
%          qsw_int(:,:,:) = qsw_int(:,:,:) + squeeze(qrasw(:,i,:,:)).*squeeze(pii(:,i,:,:)).*(-znw(i+1)+znw(i)).*squeeze(mu(:,:,:))/G;
      
    end
    middle=find(presmean<=60000,1);
    SF=trapz(mean(ght_all),qv.*rho)./trapz(mean(ght_all),qvs.*rho);
		   MDC=(trapz(mean(ght_all(:,1:middle-1)),qv(1:middle-1,:).*rho(1:middle-1,:))-trapz(mean(ght_all(:,middle:end)),qv(middle:end,:).*rho(middle:end,:)))./trapz(mean(ght_all),qvs.*rho);

		   MDCnew=trapz(mean(ght_all(:,1:middle-1)),qv(1:middle-1,:).*rho(1:middle-1,:)./qvs(1:middle-1,:))./trapz(mean(ght_all(:,1:middle-1)),ones(1,size(ght_all(:,1:middle-1),2))'.*rho(1:middle-1,:))-trapz(mean(ght_all(:,middle:end)),qv(middle:end,:).*rho(middle:end,:)./qvs(middle:end,:))./trapz(mean(ght_all(:,middle:end)),ones(1,size(ght_all(:,middle:end),2))'.*rho(middle:end,:));
    

		   MDCalt=(trapz(mean(ght_all(:,1:middle-1)),qv(1:middle-1,:).*rho(1:middle-1,:))./trapz(mean(ght_all(:,1:middle-1)),qvs(1:middle-1,:).*rho(1:middle-1,:))-trapz(mean(ght_all(:,middle:end)),qv(middle:end,:).*rho(middle:end,:))./trapz(mean(ght_all(:,middle:end)),qvs(middle:end,:).*rho(middle:end,:)));
		  
    B=trapz(mean(ght_all(:,1:middle-1)),H(1:middle-1,:).*rho(1:middle-1,:));
    T=trapz(mean(ght_all(:,middle:end)),H(middle:end,:).*rho(middle:end,:));

%%
H_prime=detrend(H','constant')';
dH = diff(H,1,2);
    
H_Profs(1,:) = [o1_ERA]'\(H_prime');
dH_Profs(1,:) = [o1_ERA_mid]'\(dH)';    
H_Profs(2,:) = [o2_ERA]'\(H_prime');
dH_Profs(2,:) = [o2_ERA_mid]'\(dH)';

test_prof1 = zeros(3,length(presmean));
test_prof2 = zeros(3,length(presmean));
resid1 = zeros(301,length(presmean));
resid2 = zeros(300,length(presmean));

stats1 = zeros(length(presmean),4);
stats2 = zeros(length(presmean),4);
for i = 1:length(presmean)
    [temp1,~,temp2,~,temp3]=regress(H_prime(i,:)',[o1_ERA;o2_ERA;ones(size(o1_ERA))]');
    test_prof1(:,i)=temp1;
    resid1(:,i)=temp2;
    stats1(i,:)=temp3;
    [temp1,~,temp2,~,temp3]=regress(dH(i,:)',[o1_ERA_mid;o2_ERA_mid;ones(size(o1_ERA_mid))]');
    test_prof2(:,i)=temp1;
    resid2(:,i)=temp2;
    stats2(i,:)=temp3;
    
end

%%
%figure1
figure
			  s1=subplot(3,1,1);
			  contourf((1:301)/4,presmean/100,OMEGA_WTG);
			  colorbar;xlim([0 301/4]);
			  set(gca,'ylim',[100 1000],'ydir','reverse');
			  sgtitle('Model Motion');
			  s2=subplot(3,1,2);
			  plot((1:301)/4,rr,'linewidth',2);
			  title('Precipitation [mm/day]');
			  s3=subplot(3,1,3);xlim([0 301/4]);
			  plot((1:301)/4,SF,'linewidth',2);
			  title('Column Relative Humidity');
			  xlabel('time (days)');xlim([0 301/4]);
			  set(s1,'units','normalized');
			  pos1=get(s1,'Position');
			  set(s2,'units','normalized');
			  set(s2,'Position',[s2.Position(1), s2.Position(2), pos1(3), s2.Position(4)]);
			  set(s3,'units','normalized');
			  set(s3,'Position',[s3.Position(1), s3.Position(2), pos1(3), s3.Position(4)]);
			  pause;

%%
%figure 2

figure,
subplot(1,2,1),plot(H_Profs(1:2,:),presmean/100,'linewidth',3);xline(0,'linewidth',2);
set(gca,'ylim',[100 1000],'ydir','reverse','xticklabel',{},'fontsize',16);legend('O1','O2','fontsize',16);
title('h profiles','fontsize',16);grid on;
ylabel('Pressure (hPa)');
subplot(1,2,2),plot(dH_Profs(1:2,:),presmean/100,'linewidth',3);xline(0,'linewidth',2);
set(gca,'ylim',[100 1000],'ydir','reverse','xticklabel',{},'yticklabel',{},'fontsize',16);
title('dh/dt profiles','fontsize',16);grid on;
pause;
%saveas(gcf,'H_profiles_projection_ones_wmean.png');


    qv_int=trapz(mean(ght_all),qv.*rho);
 h_int=trapz(mean(ght_all),H.*rho)*4/86400; 
 s_int=trapz(mean(ght_all),S.*rho)*4/86400;
 dSdt=[0 diff(s_int)];
 dHdt=[0 diff(h_int)];
    QRtot = qra_int(ist:end)*cp;
%     QRlw = mean(mean(qlw_int(ist:end,:,:),3),2)*cp;
%     QRsw = mean(mean(qsw_int(ist:end,:,:),3),2)*cp;
    rrp=detrend(rr,'constant');
    
    pp1 = [ones(size(rr)) detrend(QHvert','constant') detrend(QHhorz','constant') detrend(Qlh,'constant') detrend(Qsh,'constant') detrend(QRtot,'constant')];
    pp = [ones(size(rr)) detrend(QHvert','constant') detrend(Qlh,'constant') detrend(Qsh,'constant') detrend(QRtot,'constant')];

    [r,~,resid_r,~,stats_r]=regress(QRtot,[ones(length(o1_ERA),1) o1_ERA' o2_ERA']);
    [P_coefs,~,resid_P,~,stats_P]=regress(rr-r(1),[o1_ERA' o2_ERA']);


    
    meanTerms=[mean(QHvert) mean(QHhorz) mean(Qlh) mean(Qsh) mean(QRtot)];
    
    propogate=zeros(1,15);%QH_vert,QHvert_resid,QHvert_o1 qbar,QHvert_o2 qbar,QHvert_o1 qp,QHvert_o2 qp,QHvert_o2,QHhorz,Qlh,Qsh,Qrad,Qrad_o1,Qrad_o2,Qrad_mean,Qrad_resid,dh/dt
    propogate(1)=mean(dHdt.*QHvert)./mean(dHdt.^2);    
    propogate(2)=mean(dHdt.*QHvert_res)./mean(dHdt.^2);
    propogate(3)=mean(-dHdt.*(QH1.*pc_test(1,:)))./mean(dHdt.^2);
    propogate(4)=mean(-dHdt.*(QH2.*pc_test(2,:)))./mean(dHdt.^2);
    propogate(5)=mean(-dHdt.*(QHpc_mean(1).*pc_test(1,:)))./mean(dHdt.^2);
    propogate(6)=mean(-dHdt.*(QHpc_mean(2).*pc_test(2,:)))./mean(dHdt.^2);
    propogate(7)=mean(dHdt.*QHhorz)./mean(dHdt.^2);
    propogate(8)=mean(dHdt.*Qlh')./mean(dHdt.^2);
    propogate(9)=mean(dHdt.*Qsh')./mean(dHdt.^2);
    propogate(10)=mean(dHdt.*QRtot')./mean(dHdt.^2);
    propogate(11)=mean(dHdt.*(pc_test(1,:).*r(2)))./mean(dHdt.^2);
    propogate(12)=mean(dHdt.*(pc_test(2,:).*r(3)))./mean(dHdt.^2);
    propogate(13)=mean(dHdt.*(ones(length(dHdt),1).*r(1))')./mean(dHdt.^2);
    propogate(14)=mean(dHdt.*resid_r')./mean(dHdt.^2);
    propogate(15)=mean(dHdt.*dHdt)./mean(dHdt.^2);
   
    maintanence=zeros(1,15);%QH_vert,QHvert_resid,QHvert_o1 qbar,QHvert_o2 qbar,QHvert_o1 qp,QHvert_o2 qp,QHvert_o2,QHhorz,Qlh,Qsh,Qrad,Qrad_o1,Qrad_o2,Qrad_mean,Qrad_resid,dh/dt
    maintanence(1)=mean(detrend(h_int,'constant').*QHvert)./mean(detrend(h_int,'constant').^2);    
    maintanence(2)=mean(detrend(h_int,'constant').*QHvert_res)./mean(detrend(h_int,'constant').^2);
    maintanence(3)=mean(-detrend(h_int,'constant').*(QH1.*pc_test(1,:)))./mean(detrend(h_int,'constant').^2);
    maintanence(4)=mean(-detrend(h_int,'constant').*(QH2.*pc_test(2,:)))./mean(detrend(h_int,'constant').^2);
    maintanence(5)=mean(-detrend(h_int,'constant').*(QHpc_mean(1).*pc_test(1,:)))./mean(detrend(h_int,'constant').^2);
    maintanence(6)=mean(-detrend(h_int,'constant').*(QHpc_mean(2).*pc_test(2,:)))./mean(detrend(h_int,'constant').^2);
    maintanence(7)=mean(detrend(h_int,'constant').*QHhorz)./mean(detrend(h_int,'constant').^2);
    maintanence(8)=mean(detrend(h_int,'constant').*Qlh')./mean(detrend(h_int,'constant').^2);
    maintanence(9)=mean(detrend(h_int,'constant').*Qsh')./mean(detrend(h_int,'constant').^2);
    maintanence(10)=mean(detrend(h_int,'constant').*QRtot')./mean(detrend(h_int,'constant').^2);
    maintanence(11)=mean(detrend(h_int,'constant').*(pc_test(1,:).*r(2)))./mean(detrend(h_int,'constant').^2);
    maintanence(12)=mean(detrend(h_int,'constant').*(pc_test(2,:).*r(3)))./mean(detrend(h_int,'constant').^2);
    maintanence(13)=mean(detrend(h_int,'constant').*(ones(length(detrend(h_int,'constant')),1).*r(1))')./mean(detrend(h_int,'constant').^2);
    maintanence(14)=mean(detrend(h_int,'constant').*resid_r')./mean(detrend(h_int,'constant').^2);
    maintanence(15)=mean(detrend(h_int,'constant').*dHdt)./mean(detrend(h_int,'constant').^2);
   
%%


%%
%Figure 2

figure('units','normalized','outerposition',[.035 .48 .90 .6]),
sgtitle('Moisture Mode Criteria','fontsize',28,'fontsize',28);
subplot(1,3,1)
scatter(rrp*86400/Lv,detrend(qv_int,'constant'),[],ang_ERA_colors,'filled');
ylim([-4 4]);
xlim([-20 20]);
xlabel('P'' (mm/day)','fontsize',28);
ylabel('<q_v>'' (mm)','fontsize',28);
axis square;
subplot(1,3,2)
scatter(-detrend(QStot,'constant'),detrend(QRtot+rr+Qsh,'constant'),[],ang_ERA_colors,'filled');
ylim([-600 600]);
xlim([-1000 1000]);
xlabel('<\omega \partial S/\partial p>'' (W/m^2)','fontsize',28);
ylabel('Q_{R}''+Q_{sh}''+Q_{P}'' (W/m^2)','fontsize',28);
axis square;
subplot(1,3,3)
scatter(detrend(qv_int,'constant'),detrend(h_int,'constant')*86400/Lv/4,[],ang_ERA_colors,'filled');
ylim([-4 4]);
xlim([-4 4]);
xlabel('<q_v>'' (mm)','fontsize',28);
ylabel('<h>'' (mm)','fontsize',28);
axis square;
pause;
%saveas(gcf,'Moisture_mode_evidence_ERA_EOFs.png');

%%
%figure 7
 figure('units','normalized','outerposition',[.035 .48 .90 .75]),
 subplot(2,1,1),
 plot((1:301)/4,o1_ERA,(1:301)/4,o2_ERA,'linewidth',3);
 xlim([1 301/4]);set(gca,'fontsize',28);
 title('Vertical Motion PCs','fontsize',28);
 legend('O1','O2','fontsize',28);
 subplot(2,1,2),
 [ax1,h1,h2]=plotyy((1:301)/4,SF,(1:301)/4,MDCnew);
 set(ax1,'fontsize',28,'xlim',[0 301/4]);
 set(h1,'linewidth',3);set(h2,'linewidth',3);
 title('Moisture variability','fontsize',28);
 legend('CRH','MDC','fontsize',28);
 xlabel('Time (days)','fontsize',28);
 
 saveas(gcf,'two_mode_compare_new_MDC_RH.png');
 
%%
%figure 4
figure('units','normalized','outerposition',[0.2424    0.2759    0.2576    0.7029]),
subplot(1,2,1),
plot(lds_ERA_plot,presmean/100,'linewidth',3);
set(gca,'ylim',[100 1000],'ydir','reverse','fontsize',24);
legend('Mode 1','Mode 2');
load('toy_model_MSE.mat');
subplot(1,2,2),
plot(s,presmean,h,presmean,'linewidth',3);
legend('DSE','MSE');
set(gca,'ylim',[100 1000],'ydir','reverse','fontsize',24,'yticklabel','');
saveas(gcf,'toy_plot_diagram.png');

a=-detrend(QStot,'constant');
b=-detrend(QHtot,'constant');
c=detrend(QRtot+Qlh+Qsh,'constant');
%figure,hold on;


    colors = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];

    
    
    %%
%    close all;
%    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
%    bar([propogate(15) propogate(1) propogate(8) propogate(9) propogate(10)])%dhdt vert horz lh qra%d
%    set(gca,'TickLabelInterpreter','latex','fontsize',20)
%    xticklabels({}) 
%    title('1','fontsize',14)
%    saveas(gcf,'propogation_projection_no_label.png');
%    
%    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
%    bar([maintanence(15) maintanence(1) maintanence(8) maintanence(9) maintanence(10)])%dhdt vert h%orz lh qrad
%    set(gca,'TickLabelInterpreter','latex','fontsize',20)
%    xticklabels({'$-o_i \langle \Omega_i \frac{\partial H}{\partial p} \rangle$','$$-o_i \langle \O%mega_i \frac{\partial \overline H}{\partial p} \rangle$$','','',''});
%    title('2','fontsize',14)
%    saveas(gcf,'maintanence_budget_projection_no_label.png');
%    
%        figure('units','normalized','outerposition',[.35 .48 .36 .41]),
%    bar([propogate(1) propogate(3)+propogate(4) propogate(5)+propogate(6) propogate(10) propogate(1%1)+propogate(12)+propogate(13)])%dhdt vert horz lh qrad
%    set(gca,'TickLabelInterpreter','latex','fontsize',20)

%    xticklabels({'$-o_i \langle \Omega_i \frac{\partial H}{\partial p} \rangle$','','$$-o_i \langle \Omega_i \frac{\partial \overline H}{\partial p} \rangle$$','','$o_i R_i$'});
%    title('3','fontsize',14)
%    saveas(gcf,'propogation_projection_labels_2mode.png');
%    
%    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
%    bar([maintanence(1) maintanence(3)+maintanence(4) maintanence(5)+maintanence(6) maintanence(10) maintanence(11)+maintanence(12)+maintanence(13)])%dhdt vert horz lh qrad
%    set(gca,'TickLabelInterpreter','latex','fontsize',20)
%    xticklabels({}) 
%    title('4','fontsize',14)
%    saveas(gcf,'maintanence_budget_projection_no_label_2mode.png');
    
%        figure('units','normalized','outerposition',[.35 .48 .36 .41]),
%    bar([propogate(5)+propogate(11)+propogate(6)+propogate(12) propogate(5)+propogate(11) propogate%(6)+propogate(12)])%dhdt vert horz lh qrad
%    set(gca,'TickLabelInterpreter','latex','fontsize',20)
%    xticklabels({}) 
%    title('5','fontsize',14)
%    saveas(gcf,'propogation_projection_no_label_2mode_compare.png');
    
%    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
%    bar([maintanence(5)+maintanence(11)+maintanence(6)+maintanence(12) maintanence(5)+maintanence(1%1) maintanence(6)+maintanence(12)])%dhdt vert horz lh qrad
%    set(gca,'TickLabelInterpreter','latex','fontsize',20)
%    xticklabels({}) 
%    title('6','fontsize',14)
%    saveas(gcf,'maintanence_projection_no_label_2mode_compare.png');
    
%    return;
    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([propogate(15) propogate(1) propogate(8) propogate(9) propogate(10)])%dhdt vert horz lh qrad
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'$\frac{\partial \langle H \rangle}{\partial t}$','$-\langle \omega \frac{\partial H}{\partial p} \rangle$','$Q_{LH}$','$Q_{SH}$','$Q_R$'})    
    title('Projection of MSE budget terms onto MSE tendency','fontsize',14)
    %print(gcf,'-djpeg','../periodic_simulation/MSE_budget_projection.jpg');
    
    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([maintanence(15) maintanence(1) maintanence(8) maintanence(9) maintanence(10)])%dhdt vert horz lh qrad
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'$\frac{\partial \langle H \rangle}{\partial t}$','$-\langle \omega \frac{\partial H}{\partial p} \rangle$','$Q_{LH}$','$Q_{SH}$','$Q_R$'})    
    title('Projection of MSE budget terms onto MSE anomaly','fontsize',14)
    
    
    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([propogate(1) propogate(3)+propogate(5) propogate(4)+propogate(6) propogate(2)])%vert o1qb o2qb o1qp o2qp resid
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'total','o1','o2','residual'})    
    title('Projection of MSE budget terms onto total tendency','fontsize',14)
    %saveas(gcf,'../periodic_simulation/MSE_budget_modes_projection.png');
    
    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([maintanence(1) maintanence(3)+maintanence(5) maintanence(4)+maintanence(6) maintanence(2)])%vert o1 o2 resid
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'total','o1','o2','residual'})    
    title('Projection of MSE budget terms on column MSE anomaly','fontsize',14)
    %saveas(gcf,'../periodic_simulation/maintanence_budget_modes_projection.png');
    
    
        figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([propogate(10) propogate(11) propogate(12) propogate(14)])%total o1 o2 r0 resid
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'total','R1','R2','residual'})    
    title('Projection of MSE budget terms on total tendency','fontsize',14)
    %saveas(gcf,'../periodic_simulation/MSE_budget_radiation_projection.png');
    
    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([maintanence(10) maintanence(11) maintanence(12) maintanence(14)])%total o1 o2 r0 resid
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'total','R1','R2','residual'})     
    title('Projection of MSE budget terms on column MSE anomaly','fontsize',14)
    %saveas(gcf,'../periodic_simulation/maintanence_budget_radiation_projection.png');
   
    
    
        figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([propogate(10)+propogate(1) propogate(11)+propogate(3)+propogate(5) propogate(12)+propogate(4)+propogate(6) propogate(14)+propogate(2)])%vert o1 o2 resid
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'total','o1','o2','residual'})    
    title('Projection of MSE budget terms on total tendency','fontsize',14)
    %saveas(gcf,'../periodic_simulation/MSE_budget_combined_projection.png');
    
    figure('units','normalized','outerposition',[.35 .48 .36 .41]),
    bar([maintanence(10)+maintanence(1) maintanence(11)+maintanence(3)+maintanence(5) maintanence(12)+maintanence(4)+maintanence(6) maintanence(14)+maintanence(2)])%vert o1 o2 resid
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    xticklabels({'total','o1','o2','residual'})     
    title('Projection of MSE budget terms on column MSE anomaly','fontsize',14)
    %saveas(gcf,'../periodic_simulation/maintanence_budget_combined_projection.png');

return;

    qpos=(qbot+qtop)/sqrt(2);
    qneg=-(qbot-qtop)/sqrt(2);
    

figure('units','normalized','outerposition',[.0635 .310 .65 .575]);
subplot(2,1,1),
plot(t,detrend(qpos,'constant'),t,detrend(qneg,'constant'),'linewidth',3.5)
title('Combined and Mixed moisture index');
xlim([0 t(end)]);
set(gca,'fontsize',20);
xlabel('days');
legend('B+T','B-T');
subplot(2,1,2),
plot(t,pcs_scaled(1,:),t,pcs_scaled(2,:),'linewidth',3.5)
xlim([0 t(end)]);
legend('O1','O2');
set(gca,'fontsize',20);
   % saveas(gcf,'../periodic_simulation/mixed_index.png');
    
figure('units','normalized','outerposition',[.0635 .310 .65 .575]);
subplot(2,1,1),
plot(t,detrend(qbot,'constant'),t,detrend(qtop,'constant'),'linewidth',3.5)
title('Top Moisture Index vs Bottom Moisture Index');
xlim([0 t(end)]);
set(gca,'fontsize',20);
xlabel('days');
legend('P<600 hPa','P>600 hPa');
subplot(2,1,2),
plot(t,pcs_scaled(1,:),t,pcs_scaled(2,:),'linewidth',3.5)
xlim([0 t(end)]);
legend('O1','O2');
set(gca,'fontsize',20);
   % saveas(gcf,'../periodic_simulation/top_bottom_index.png');
  
    
