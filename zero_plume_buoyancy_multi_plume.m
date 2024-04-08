function [s,henv,henv_sat,hp1,tp1] = zero_plume_buoyancy_multi_plume(T,qv,z,p,epsilon_base)

qvstar=calc_qvstar(T,p);

epsilon = 1./((epsilon_base)').*ones(length(4000),length(z)); 
%epsilon = .75./z';

hp1=zeros(size(epsilon,1),length(z));qp1=hp1;tp1=hp1;hp2=zeros(length(z),1);tp2=hp2;
if size(tp1,1)~=size(T,1)
tp1=tp1';hp1=hp1';
tp2=tp2';hp2=hp2';
end

%epsilon2 = (1/1000)*ones(size(z));


Lv = 2.5e6;
cp = 1005;
g=9.81;

s=cp*T+z*g;
henv=(cp*T+Lv*qv+z*g);
henv_sat=(cp*T+Lv*qvstar+z*g);


zt=find(henv(1)>=henv_sat);

if isempty(zt)

   zt=find(z<=1000);
   zt=zt(end);

end

zt=zt(1);
if z(zt)>=4000
   zt=find(z<=950);
   zt=zt(end);
end


hp1(1:zt-1)=repmat(henv_sat(1:zt-1),1,size(epsilon,1));
hp1(zt)=repmat(henv(1),1,size(epsilon,1));
[tp1(zt)]=plume_temperature(hp1(zt)-henv(zt),T(zt),p(zt));
qp1(zt)=calc_qvstar(tp1(zt),p(zt));



%hp(1:z500(end)-1)=henv(1:z500(end)-1);
for index = zt:length(z)-1
hp1(index+1) = hp1(index) - epsilon((index)).*(hp1(index)-repmat(henv(index),size(epsilon,1),1)).*repmat((z(index+1)-z(index)),size(epsilon,1),1);

[tp1(index+1)]=plume_temperature(hp1(index+1)-henv_sat(index+1),T(index+1),p(index+1));

qp1(index+1)=calc_qvstar(tp1(index+1),p(index+1));
end


hp1(1:zt)=henv(1);
tp1(1:zt)=T(1:zt);
% figure,
% plot(henv_sat,p,henv,p,'linewidth',3);
% hold on;
% plot(hp1,p,'k:','linewidth',3);
% xline(henv(1),'linewidth',3);
% set(gca,'ydir','reverse','ylim',[min(p) max(p)],'fontsize',16);
% legend('saturate MSE','MSE','plume MSE');
% xlabel('energy (J/kg)');
% ylabel('pressure (Pa)');
%print(gcf,'-dpng','plume_mse_comparison.jpg');
end