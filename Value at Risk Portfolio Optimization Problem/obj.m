function [c,ceq,gradc,gradceq] =obj(x)
global pLevel BestModel a B f0 eps;
pi=BestModel.PComponents;
f=a+B'*x;

y=0;
 for k=1:length(pi)
     X=-(f0+f'*BestModel.mu(k,:)')/sqrt(f'*BestModel.Sigma(:,:,k)*f);
    y=y+pi(k)*normcdf(X);
    y
 end  

c=pLevel-y;
ceq = [];
if nargout > 2
    try
    gradc = -gradient(x);
    gradceq = [];
    catch 
        gradient(x)
    end
end