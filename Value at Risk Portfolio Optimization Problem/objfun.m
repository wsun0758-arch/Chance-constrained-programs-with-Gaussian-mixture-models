function y= objfun(x)
global BestModel B;
n=size(B,1);
try
y=-sum(BestModel.PComponents*BestModel.mu*x(1:n));
catch
    y
end



end
