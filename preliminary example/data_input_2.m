% SeDuMi code for marginal risk constrained portfolio optimization: data
% input

function [At,c,K,lwz0,upz0]=data_input_2(model,a,B,alpha,ub);  

n=size(B,1);
m=size(B,2);
num=model.NComponents;

  % ---------- constraints c-A^T*x >= 0) ----------
  c = sparse([ub*ones(n,1);zeros(n,1);-1+alpha; zeros(num,1);ones(num,1)]);  

     At = sparse([eye(n), zeros(n,num)]);
      At = sparse([At;-eye(n), zeros(n,num)]); % x > lwx ~ -x < -lwx 
  At = sparse([At; ...
               zeros(1,n), -model.PComponents]); % pi^T*x >= 1-alpha
           
  At = sparse([At; ...
               zeros(num,n), -eye(num)]); % y > lwy ~ -y < -lwy 
    
  At = sparse([At; ...
               zeros(num,n), eye(num)]); % y < upy   
 
   L=zeros(m,m,num);
  for k=1:num
      L(:,:,k) = chol(model.Sigma(:,:,k),'lower');
      c = sparse([c;0;L(:,:,k)*a]);
  end
  %\lambda_k x^T u_ku_k' x + l(z_k) < \rho_k 
  Atp=zeros(num*(1+m),n+num); 
  for k=1:num
     Atp((k-1)*(m+1)+2:(k)*(m+1),1:n) = -(B*L(:,:,k))';      
  end
  
  At= sparse([At; Atp]);         
  
  clear('Atp');
      
  %-------- inequality constraint -----------
  K.l=[2*n+1+2*num];

  %-------- quadratic constraint -----------
  K.q=[];
  for k=1:num
      K.q=[K.q,m+1];
  end
lwz0=zeros(1,num);
upz0=ones(1,num); 