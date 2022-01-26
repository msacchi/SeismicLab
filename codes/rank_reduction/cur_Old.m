function [C,U,R]=cur(M,p,c);
%
%   CUR is a matrix decomposition where A = C*U*R.       
%   - Columns of C are randomly picked with the probability function based
%   on the L2 norm of each column. The bigger the norm, the higher the 
%   chance of being picked. Zeroed columns are not likely to be picked. 
%   - Columns of R are picked the same way as columns of C
%   - U = pseudoinv(C)*d0*pseudoinv(R)
%
%
%   d0     -> input data
%   nx, ny -> matrix dimensions
%   nc,nr  -> dimensions of wanted C, U and R, where C(nx,nc), U(nc,nr) 
%             and R(nr,ny). Usually we're using nc = nr, but there's
%             nothing that forbides you from using nc =/= nr
%   c      -> weighting factor applied on the probability function. On the
%             function it is default 0.2. c usually is in [0,1] interval, 
%             with smaller values giving us better results 
%             
%
%
%
%

 c = 0.2;
 C = zeros(p,p);
 R = zeros(p,p);
 [nx,ny] = size(M); 

 ex  = sum (abs(M).^2, 2);
 ey  = sum (abs(M).^2, 1);
 ex = sqrt(ex);
 ey = sqrt(ey);

 exmax = max(ex);
 eymax = max(ey);

 count_R = 0;

 ix = 1;

while ix < nx+1
   b = exmax*rand;
   b2 = c*ex(ix);
   if (b < b2)
      ex(ix) = -1 ;
      count_R = count_R + 1;
   end
   if count_R == nr
      break;
   end
   ix = ix+1;
   if (count_R < nr && ix == nx)
      ix = 1;
   end
end       
 
for ix = 1:nx
   if ex(ix) < 0
      R = [R; squeeze(M(ix,:))];
   end
end

count_C = 0;    
iy = 1;    

while iy < ny+1
   b = eymax*rand;
   if (b < c*ey(iy))
      ey(iy) = -1;
      count_C = count_C + 1;
   end
   if count_C == nc
      break;
   end
   iy = iy+1;
   if (count_C < nc && iy == ny)
      iy = 1;
   end
end 
    
for iy = 1:ny
   if ey(iy) < 0
      C = [C; M(:,iy)'];
   end
end

C = C';
C_inv = pinv(C);
R_inv = pinv(R);
U = C_inv*d0*R_inv;

return;
