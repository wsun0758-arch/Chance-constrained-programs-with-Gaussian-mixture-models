v=[2.5/100:0.1/100: 6.5/100];
figure
for j=1:length(alpha)
    subplot(1,3,j);
    obj=[];
     for i=1:size(objective_value,1) 
            obj=[obj,objective_value(i,j)];
     end
     
     temp0=plot([2.5/100:0.1/100: 6.5/100],obj,'o');hold on;
     
     
     obj=[];
     for i=1:size(objective_value,1) 
             if EXITFLAG(i,j)~=-2
                obj=[obj,-FVAL(i,j)];
               if j==length(alpha) && i==size(objective_value,1) 
                 temp=plot(v(i),-FVAL(i,j),'rx');
                  legend([temp0,temp],'Branch and Bound Algorithm','Nonlinear Optimization Algorithm');
               else
                 plot(v(i),-FVAL(i,j),'rx');
               end
             end
     end
  
     xlabel('v');
     ylabel('Expected Return Rate');
     
end
  
    hold off;



figure
for j=1:length(alpha)
     subplot(1,3,j);
     plot([2.5/100:0.1/100: 6.5/100],spend_time(:,j)); 
     xlabel('v');
     ylabel('CPU time (seconds)');
end
     
