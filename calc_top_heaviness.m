function [th_angle,pcs] = calc_top_heaviness(omega,pres)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%load the era5 eofs and the colormap
load('ERA5_EOFS.mat');
lds_ERA=lds;

presmean=mean(pres,1);
dp=diff(pres,1,2);
weightsShapeFull=[dp dp(:,end)]./sum([dp dp(:,end)]);
weightsShapeFull(weightsShapeFull<=0)=-weightsShapeFull(weightsShapeFull<=0);
omega2=omega.*(weightsShapeFull.^0.5);

%interpolate and project the vertical motion onto the EOFs to get the PCs
height_index=find(level>=presmean(end));
lds_ERA=lds_ERA(height_index,:);
lds_ERA(:,1)=-lds_ERA(:,1);
lds_ERA(:,2)=lds_ERA(:,2);
omega2_interp=interp1(presmean',omega2',level(height_index),[],'extrap');
pcs_ERA=lds_ERA'*omega2_interp;
lds_ERA_interp=interp1(level(height_index),lds_ERA,presmean,[],'extrap');
lds_ERA_weight=lds_ERA_interp./(mean(weightsShapeFull)'.^0.5);
a=sqrt(trapz(presmean,lds_ERA_weight.^2)./(presmean(end)-presmean(1)));
eofs=lds_ERA_weight./a;
pcs=pcs_ERA.*a';



th_angle=atan2d(pcs(2,:),pcs(1,:));


end

