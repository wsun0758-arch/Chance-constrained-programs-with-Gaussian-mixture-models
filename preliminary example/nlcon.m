function cx = nlcon(x)
global BestModel a B lwz;
num=BestModel.NComponents;  
  n=size(B,1);
  f=a+B'*x(1:n);
  for k=1:num
      try
      cx(k)=-(norminv(full(lwz(k)))*sqrt(f'*BestModel.Sigma(:,:,k)*f)+f'*BestModel.mu(k,:)');
      catch
          lwz
      end
  end
end

