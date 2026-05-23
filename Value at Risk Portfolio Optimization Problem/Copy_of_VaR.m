figure
for j=1:length(alpha)
     subplot(1,3,j);
     Var=[];
     v=[];
     for i=1:size(objective_value,1) 
         v=[v,mean(sample*solution_opt(:,i,j))];
         Var=[Var,quantile(-sample*solution_opt(:,i,j),1-alpha(j))];
     end
     plot(Var,v,'o');hold on;
     xlabel('VaR_{1-\alpha}');
     ylabel('Mean Return Rate');
end


for j=1:length(alpha)
     subplot(1,3,j);
     Var=[];
     v=[];
     for i=1:size(objective_value,1)
           v=[v,mean(sample*solution_opt0(:,i,j))];
           Var=[Var,quantile(-sample*solution_opt0(:,i,j),1-alpha(j))];
     end
     plot(Var,v,'x');hold on;
     xlabel('VaR_{1-\alpha}');
     ylabel('Mean Return Rate');
end
hold off;