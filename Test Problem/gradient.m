function y=gradient(x)
global a B BestModel f0;
f=a+B'*x;
pi=BestModel.PComponents  ;
y=0;
try
    for k=1:length(pi)
         X=-(f0+f'*BestModel.mu(k,:)')/sqrt(f'*BestModel.Sigma(:,:,k)*f);
        y0=normpdf(X)*(f'*BestModel.mu(k,:)'*B*BestModel.Sigma(:,:,k)*f-B*BestModel.mu(k,:)'*f'*BestModel.Sigma(:,:,k)*f);
        y0=y0/((f'*BestModel.Sigma(:,:,k)*f)^1.5);
         y=y+pi(k)*y0;
    end
catch
    y
end
