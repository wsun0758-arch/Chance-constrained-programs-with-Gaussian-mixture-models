clear;
clc;
load matlab2
global BestModel
X=sample;
maxk=10;
repetition=500;
options = statset('MaxIter',10000); % Increase number of EM iterations
RegularizationValue = 1e-4;
index=zeros(1,repetition);
  gmm=cell(1,maxk-1);
  AIC=zeros(1,maxk-1);

 for k=2:maxk
     tempgmm=cell(1,repetition);
     tempmle=zeros(1,repetition);
      for j=1:repetition
         tempgmm{j}=fitgmdist(X,k,'Options',options,'RegularizationValue',RegularizationValue);
         tempmle(j)=tempgmm{j}.AIC;
      end
      [~,I]=min(tempmle);
      AIC(k-1)=tempgmm{unique(I)}.AIC;
      gmm{k-1}=tempgmm{I};
  end
  [~,I] = min(AIC);
I+1
BestModel=gmm{I};
singleGaussian=fitgmdist(sample,1,'Options',options,'RegularizationValue',RegularizationValue);