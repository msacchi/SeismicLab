function out = NMO_operator(in,Param,flag);
% NMO_operator: performs forward and adjoing NMO operator
% 
% h:  offset
% dt: sampling interval
% v:  nmo-velocity
%
% In:  Data after NMO if flag==1;
%      Data before NMO if flag ==-1 
%
% out: Data before NMO if flag==1
%      Data after NMO if flag==-1 
%

   h = Param.h;
  dt = Param.dt;
vnmo = Param.vnmo; 

     v2 = vnmo.^2;
     h2 = h.^2;
[nt,nx] = size(in);

if flag==1; m = in; d = zeros(nt,nx); 

  for k = 1:nx
   for it0 = 1:nt
     t0 = (it0-1)*dt;
     t = sqrt(t0^2 + h2(k)/v2(it0));
     its = t/dt+1;
     it1 = floor(its);
     it2 = it1+1;
     a = its-it1;
    if it2<=nt; 
     d(it1,k) = d(it1,k) + (1-a)*m(it0,k); 
     d(it2,k) = d(it2,k) +     a*m(it0,k); 
   end
  end
 end

out = d;

else 

  d = in; ma = zeros(nt,nx); 

 for k = 1:nx
  for it0 = 1:nt
     t0 = (it0-1)*dt;
     t = sqrt(t0^2 + h2(k)/v2(it0));
     its = t/dt+1;
     it1 = floor(its);
     it2 = it1+1;
     a = its-it1;
    if it2<=nt; ma(it0,k) = (1.-a)*d(it1,k)+a*d(it2,k); end
   end
 end

out = ma;

end

