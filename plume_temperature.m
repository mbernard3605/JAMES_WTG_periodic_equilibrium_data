function [Tp] = plume_temperature(Hpdiff,Tenv,P)
%lookup temperature and saturation specific humidity from MSE minus
%geopotential and pressure
cp=1005.4;
Lv=2.5e6;
qvstar=calc_qvstar(Tenv,P);
A=cp+Lv.*qvstar.*(17.67*243.5)./(Tenv-273.15+243.5).^2;
Tp=Tenv+Hpdiff./A;



end

