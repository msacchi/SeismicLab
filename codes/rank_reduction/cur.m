function [DD]=cur(M,p);

  [n1,n2]=size(M);

  e  =  sum (abs(M).^2, 2);
  normalization = sum(e);
  Px = e/normalization ;

  e  = sum (abs(M).^2, 1);
  normalization = sum(e);
  Py = e'/normalization ;

  i1 = datasample(1:n1,p,'Replace', false, 'Weights',Px);
  i2 = datasample(1:n2,p,'Replace', false, 'Weights',Py);


 C = M(:,i2);
 R = M(i1,:);
 C_inv = pinv(C);
 R_inv = pinv(R);
 U = C_inv*M*R_inv;

 DD=C*U*R;

