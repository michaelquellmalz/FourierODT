function [xx1,xx2,xx3] = rotation_axis(x1,x2,x3,n1,n2,n3,a)
  s = size(x1.*n1);
  x1 = x1(:);
  x2 = x2(:);
  x3 = x3(:);
  n1 = n1(:);
  n2 = n2(:);
  n3 = n3(:);
  a = a(:);
  xx1 = (n1.*n1.*(1-cos(a))    +cos(a)).*x1...
       +(n1.*n2.*(1-cos(a))-n3.*sin(a)).*x2...
       +(n1.*n3.*(1-cos(a))+n2.*sin(a)).*x3;
  xx2 = (n2.*n1.*(1-cos(a))+n3.*sin(a)).*x1...
       +(n2.*n2.*(1-cos(a))+    cos(a)).*x2...
       +(n2.*n3.*(1-cos(a))-n1.*sin(a)).*x3;
  xx3 = (n3.*n1.*(1-cos(a))-n2.*sin(a)).*x1...
       +(n3.*n2.*(1-cos(a))+n1.*sin(a)).*x2...
       +(n3.*n3.*(1-cos(a))+    cos(a)).*x3;
  xx1 = reshape(xx1,s);
  xx2 = reshape(xx2,s);
  xx3 = reshape(xx3,s);
end