load matlab2
global BestModel a B b f0 pLevel eps;
clc;
alpha=[0.1,0.2,0.3,0.4];
maxnum=length(alpha);
objective_value=zeros(replication,numd,maxnum);
spend_time=zeros(replication,numd,maxnum);
FVAL=zeros(replication,numd,maxnum);
EXITFLAG=zeros(replication,numd,maxnum);
f0=-20;
ub=5;
eps=1.0e-7;
opts = optimset('Algorithm','sqp','GradConstr','on','MaxFunEvals',10000);
opts.ConstraintTolerance=eps;


for k=1:maxnum
    pLevel=1-alpha(k);
     for j=1:numd
        for i=1:replication
            BestModel=model{i,j,k};
            b=objpara{i,j};
            m=size(BestModel.mu,2);
            d=length(b);
            a=zeros(m,1);
            B=eye(m);%%%d*k
            [x0,FVAL(i,j,k),EXITFLAG(i,j,k)]=fmincon(@objfun,ub*ones(d,1)/2,[],[],ones(1,d),1,zeros(d,1),ones(d,1),@obj,opts);
            f=a+B'*x0;
            tempx0=x0;
            for num=1:BestModel.NComponents  
                 X=-(f0+f'*BestModel.mu(num,:)')/sqrt(f'*BestModel.Sigma(:,:,num)*f); 
                 tempx0=[tempx0;normpdf(X)];
            end
           [solution_opt,objective_value(i,j,k),flag,spend_time(i,j,k)]=branch_bound(BestModel,a,B,b,f0,alpha(k),tempx0,ub)
        end
    end
end
time=[];time0=[];
for k=1:maxnum
     for j=1:numd
        for i=1:replication
            if EXITFLAG(i,j,k)==-2
                FVAL(i,j,k)=0;
                time=[time, spend_time(i,j,k)];
            else
                time0=[time0,spend_time(i,j,k)];
            end
        end
     end
end
