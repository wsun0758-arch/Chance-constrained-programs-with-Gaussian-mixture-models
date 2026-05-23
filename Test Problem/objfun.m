function y= objfun(x)
global B b;
n=size(B,1);
y=sum(b*x(1:n));




end
