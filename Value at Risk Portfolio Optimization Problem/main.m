global BestModel a B f0 pLevel eps;
clc;
m=size(BestModel.mu,2);%Gaussian Mixture Model fitted using EM algorithm based on real data
d=size(sample,2);%dimension of decision vector
%%%%%%%%%%%%%%%%
%Pr{f0+f(x)^T \xi<=0}>1-alpha
%f(x)=(f1,...,fk)^T
%f(x)=a+B^T*x 
a=-zeros(m,1);
B=-eye(m);%%%B is d*k matrix
f0=-2.5/100;%f0=-v %%v refers to Right Hand Side of VaR constraint(14)
%f0=-3.8/100;%%
interv=41;%%the nummber of various Right Hand Side of VaR constraint(14)
eps=1.0e-7;
alpha=[0.01,0.05,0.1];%%the nummber of various alpha value in VaR constraint(14)
%%%%%%%%
%Parameter setting for fmincon function
%Using Sequential Quadratic Algorithm for which users themselves input 
%gradient function of Nonlinear constraint(deterministic counterpart of VaR constraint(14))
opts = optimset('Algorithm','sqp','GradConstr','on','MaxFunEvals',10000);
opts.ConstraintTolerance=eps;%the Error involved in the evaluation of 
%probability function of CCP(14)
%%%%%%%%
solution_opt=zeros(d,interv,3);
objective_value=zeros(interv,3);
spend_time=zeros(interv,3);
x0=zeros(d,interv,3);
FVAL=zeros(interv,3);
EXITFLAG=zeros(interv,3);
for i=1:interv
    for j=1:length(alpha)
        pLevel=1-alpha(j);
        %%intialize the Branch and Bound algorithm by the local optima
        %%obtained using fmincon function
        [x0(:,i,j),FVAL(i,j),EXITFLAG(i,j)]=fmincon(@objfun,ones(d,1)/2,[],[],ones(1,d),1,zeros(d,1),ones(d,1),@obj,opts);
        f=a+B'*x0(:,i,j);
        tempx0=x0(:,i,j);%local optima obtained by fmincon function
        for k=1:BestModel.NComponents  
             X=-(f0+f'*BestModel.mu(k,:)')/sqrt(f'*BestModel.Sigma(:,:,k)*f); 
             tempx0=[tempx0;normpdf(X)];%% initial auxilary variables y
        end
        %%lanuch branch_bound combined with Baron and Sedumi
       [solution_opt(:,i,j),objective_value(i,j),flag,spend_time(i,j)]=branch_bound(BestModel,a,B,f0,alpha(j),tempx0)
    end  
      f0=f0-0.1/100;
end

