global BestModel a B f0 pLevel eps;
clc;
m=size(BestModel.mu,2);
d=size(sample,2);
a=-zeros(m,1);
B=-eye(m);%%%d*k
f0=-2.5/100;%-v
%f0=-3.8/100;%%
interv=41;
eps=1.0e-7;
alpha=[0.01,0.05,0.1];
opts = optimset('Algorithm','sqp','GradConstr','on','MaxFunEvals',10000);
opts.ConstraintTolerance=eps;
%solution_opt=zeros(d,interv,3);
%objective_value=zeros(interv,3);
%spend_time=zeros(interv,3);
%x0=zeros(d,interv,3);
%FVAL=zeros(interv,3);
%EXITFLAG=zeros(interv,3);
for i=1:interv
    for j=1:length(alpha)
        pLevel=1-alpha(j);
        [x0(:,i,j),FVAL(i,j),EXITFLAG(i,j)]=fmincon(@objfun,ones(d,1)/2,[],[],ones(1,d),1,zeros(d,1),ones(d,1),@obj,opts);
        f=a+B'*x0(:,i,j);
        tempx0=x0(:,i,j);
        for k=1:BestModel.NComponents  
            X=-(f0+f'*BestModel.mu(k,:)')/sqrt(f'*BestModel.Sigma(:,:,k)*f); 
            tempx0=[tempx0;normpdf(X)];
       end
        [solution_opt(:,i,j),objective_value(i,j),flag,spend_time(i,j)]=branch_bound(BestModel,a,B,f0,alpha(j),tempx0)
    end  
      f0=f0-0.1/100;
end
figure
 for i=1:4
     subplot(2,2,i);
     plot([feasibleiter(i)+1:iter_time(i)],current_objective_value{1,i});
      xlabel('Iterations');
     ylabel('Objective Function Value');
 end
