clear;
clc;

replication=5;
numd=3;
num=2;
model=cell(replication,numd,maxnum-1);
objpara=cell(replication,numd);
for k=1:4
    for d=10:10:10*numd
        for i=1:replication
            mu=unifrnd(0,5,num,d);
            sigma=zeros(d,d,num);
            for j=1:num
                temp=unifrnd(0,1,d,d);
                sigma(:,:,j)=temp*temp';
            end
            p=rand(1,num);
            p=p./sum(p);
            model{i,d/10,k}=gmdistribution(mu,sigma,p);
            objpara{i,d/10}=unifrnd(-10,0,1,d);
        end
    end
end
