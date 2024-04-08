function qvstar = calc_qvstar(t,p)
%T in kelvin
%p in pa
if size(p,2)~=size(t,2)
   p=repmat(p,1,size(t,2)); 
end

es=611.2*exp(17.67*(t-273.15)./(t-29.65));
rs=.622*es./(p-es);
qvstar=rs./(1+rs);



end

