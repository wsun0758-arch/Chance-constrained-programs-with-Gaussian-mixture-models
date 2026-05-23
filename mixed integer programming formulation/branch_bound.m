function [solution_opt,objective_value,flag,spend_time]= branch_bound(model,a,B,f0,alpha,dz,ub)
global lwz eps;
eps0=1.0e-5;

n=size(B,1);
%m=size(B,2);
num=model.NComponents;
%b=[(model.PComponents*model.mu)';zeros(num,1)];
%Objective
fun = @(x) objfun(x);
%Nonlinear Constraints
NonlinearConstraints =@nlcon;
%Nonlinear Constraints Bounds
cl = f0*ones(num,1);
cu = ones(num,1)*Inf;
%Bounds on decision variables
lb = [zeros(n,1);zeros(num,1)];
ub = [Inf*ones(n,1);ones(num,1)];
%Bound on Linear Inequality
rl=[ub;1-alpha];
ru=[1;Inf];
%linear matrix
A=[[ones(1,n),zeros(1,num)];[zeros(1,n),model.PComponents]];
ctype=[repmat('I',dz),repmat('C',n+num-dz)];

lwz0=zeros(1,num);
upz0=ones(1,num); 

options = baronset('EpsA',1e-3,'EpsR',1e-3);
initialtime=cputime;

%initialization

v_opt = 1.0e+10;         % current optimal value 
v_low = -1.0e+10;            
x_opt = zeros(n+num,1);  % curent optimal solution 

Omega_delta_lwz = lwz0+1e-5; % 
Omega_delta_upz = upz0; %
Omega_v_low = sparse([v_low]);   % lower bound vector 

flag=0;      %  flag==1 means the problem is feasible
Omega_k=1;   %  Omega_k==0 terminates branch-bound procedure

iter_time=0;
 
initial=0;   %  initial==0 means the first time to solve sub-problem


while (Omega_k>0 && abs(v_opt-v_low)>eps0)
      branch_feasible=1;
    if initial==0  % the first time to solve sub-problem
        initial=1;
        i_minv=1;
         lb(n+1:n+num) = Omega_delta_lwz;
         ub(n+1:n+num) = Omega_delta_upz;
         lwz= Omega_delta_lwz ;
         [x,y,ef] = baron(fun,A,rl,ru,lb,ub,NonlinearConstraints,cl,cu,[],x0,options) ;    
         if ef~=1
            branch_feasible=0;
         end
    else           % update c and At for sub-problem for sequence sub-problem 
        min_v_low=1e+10; 
        for j=1:Omega_k
            if min_v_low>Omega_v_low(j)
                min_v_low=Omega_v_low(j);
                i_minv=j; 
            end
        end
        v_low=min_v_low;
        lwz=Omega_delta_lwz(i_minv,:)';
        upz=Omega_delta_upz(i_minv,:)';
        
        lb(n+1:n+num) = lwz;
        ub(n+1:n+num) = upz;
        [x,y,ef] = baron(fun,A,rl,ru,lb,ub,NonlinearConstraints,cl,cu,ctype,[],options) ;    
      
         if ef~=1
            branch_feasible=0;
         end
 
    end
     
   
    %x'
    %inf
    %Omega_delta_lwz
    %Omega_delta_upz
    %Omega_v_low
    
    if branch_feasible==0    % if the branch-problem is infeasibe
        if  Omega_k>i_minv
            Omega_delta_lwz(i_minv,:)=Omega_delta_lwz(Omega_k,:);
            Omega_delta_upz(i_minv,:)=Omega_delta_upz(Omega_k,:);
            Omega_v_low(i_minv)=Omega_v_low(Omega_k);
        end
        Omega_k=Omega_k-1; 
        if Omega_k>0
            temp_Omega_delta_lwz=Omega_delta_lwz(1:Omega_k,:);
            temp_Omega_delta_upz=Omega_delta_upz(1:Omega_k,:);
            temp_Omega_v_low=Omega_v_low(1:Omega_k);
            
            Omega_delta_lwz=temp_Omega_delta_lwz;
            Omega_delta_upz=temp_Omega_delta_upz;
            Omega_v_low=temp_Omega_v_low;
        end
    else                % (inf.dinf==0) if the sub-problem is feasibe
        v_opt_i=y;
      
        if v_opt_i>v_opt % cut this branch
            if  Omega_k>i_minv
                Omega_delta_lwz(i_minv,:)=Omega_delta_lwz(Omega_k,:);
                Omega_delta_upz(i_minv,:)=Omega_delta_upz(Omega_k,:);
                Omega_v_low(i_minv)=Omega_v_low(Omega_k);
            end
            Omega_k=Omega_k-1; 
            if Omega_k>0
                temp_Omega_delta_lwz=Omega_delta_lwz(1:Omega_k,:);
                temp_Omega_delta_upz=Omega_delta_upz(1:Omega_k,:);
                temp_Omega_v_low=Omega_v_low(1:Omega_k);
            
                Omega_delta_lwz=temp_Omega_delta_lwz;
                Omega_delta_upz=temp_Omega_delta_upz;
                Omega_v_low=temp_Omega_v_low;
            end
        else        % try to find an improved solution to the original problem
            xpi=x(1:n);
            sub_feasible=1;
            f=a+B'*xpi;
                  ddd_P=0;
                 for k=1:num
                       X=-(f0+f'*model.mu(k,:)')/sqrt(f'*model.Sigma(:,:,k)*f);
                       ddd_P=ddd_P+model.PComponents(k)*normcdf(X);                  
                 end  
                if 1-alpha-ddd_P>eps
                    sub_feasible=0; 
                end
   
          
            if sub_feasible==1 % feasbile solution to original problem
               flag=1;
                x_opt=x(1:n);
                v_opt=v_opt_i;
                if  Omega_k>i_minv
                    Omega_delta_lwz(i_minv,:)=Omega_delta_lwz(Omega_k,:);
                    Omega_delta_upz(i_minv,:)=Omega_delta_upz(Omega_k,:);
                    Omega_v_low(i_minv)=Omega_v_low(Omega_k);
                end
                Omega_k=Omega_k-1; 
                if Omega_k>0
                    temp_Omega_delta_lwz=Omega_delta_lwz(1:Omega_k,:);
                    temp_Omega_delta_upz=Omega_delta_upz(1:Omega_k,:);
                    temp_Omega_v_low=Omega_v_low(1:Omega_k);
            
                    Omega_delta_lwz=temp_Omega_delta_lwz;
                    Omega_delta_upz=temp_Omega_delta_upz;
                    Omega_v_low=temp_Omega_v_low;
                end

                
                h=0;
                glag=0;
                
                         
                for i=1:Omega_k   % cut branches
                    if Omega_v_low(i)<v_opt-eps 
                        if glag==0
                            temp_Omega_delta_lwz=Omega_delta_lwz(i,:);
                            temp_Omega_delta_upz=Omega_delta_upz(i,:);
                            temp_Omega_v_low=Omega_v_low(i);
                            glag=1;
                        else
                            temp_Omega_delta_lwz=sparse([temp_Omega_delta_lwz;...
                                                         Omega_delta_lwz(i,:)]);
                            temp_Omega_delta_upz=sparse([temp_Omega_delta_upz;...
                                                         Omega_delta_upz(i,:)]);
                            temp_Omega_v_low=sparse([temp_Omega_v_low;...
                                                     Omega_v_low(i)]);
                        end
                        h=h+1;
                    end
                end
                Omega_k=h;
                if Omega_k>0
                    Omega_delta_lwz=temp_Omega_delta_lwz;
                    Omega_delta_upz=temp_Omega_delta_upz;
                    Omega_v_low=temp_Omega_v_low;
                end
            else
          % subdivide rectangle \Delta^i
                def_z=Omega_delta_upz(i_minv,:)-Omega_delta_lwz(i_minv,:);
                max_d=0;
                max_i=0;
                for i=1:num
                    if def_z(i)>max_d
                        max_d=def_z(i);
                        max_i=i;
                    end
                end

                %def_z
               
                %Omega_delta_lwz(i_minv,:)
                %Omega_delta_upz(i_minv,:)
                
                Omega_delta_lwz=sparse([Omega_delta_lwz; Omega_delta_lwz(i_minv,:)]);
                Omega_delta_lwz=sparse([Omega_delta_lwz; Omega_delta_lwz(i_minv,:)]);
                                                         
                Omega_delta_upz=sparse([Omega_delta_upz; Omega_delta_upz(i_minv,:)]);
                Omega_delta_upz=sparse([Omega_delta_upz; Omega_delta_upz(i_minv,:)]);
                
                Omega_v_low=sparse([Omega_v_low; v_opt_i]);
                Omega_v_low=sparse([Omega_v_low; v_opt_i]);   
                
                Omega_k=Omega_k+2;
                
               
                %i_minv
                %max_i
                
                aa=(Omega_delta_upz(i_minv,max_i)+Omega_delta_lwz(i_minv,max_i))/2;                               
                                                     
                Omega_delta_lwz(Omega_k,max_i)=aa;              
                Omega_delta_upz(Omega_k-1,max_i)=aa;
                
                Omega_delta_lwz(i_minv,:)=Omega_delta_lwz(Omega_k,:);
                Omega_delta_upz(i_minv,:)=Omega_delta_upz(Omega_k,:);
                Omega_v_low(i_minv)=Omega_v_low(Omega_k);
         
                Omega_k=Omega_k-1; 

                temp_Omega_delta_lwz=Omega_delta_lwz(1:Omega_k,:);
                temp_Omega_delta_upz=Omega_delta_upz(1:Omega_k,:);
                temp_Omega_v_low=Omega_v_low(1:Omega_k);
            
                Omega_delta_lwz=temp_Omega_delta_lwz;
                Omega_delta_upz=temp_Omega_delta_upz;
                Omega_v_low=temp_Omega_v_low;
            end
        end
    end
    
    iter_time=iter_time+1  
  if iter_time>300
      iter_time
  end
    %for k=1:m 
    %    real_marginal_risks(k)=lam_gam(1,k)*xpi'*U(:,k)*U(:,k)'*xpi-lam_gam(2,k)*xpi'*V(:,k)*V(:,k)'*xpi;
    %end
    %marginal_risk_rho=rho
    %real_marginal_risks
    
    v_opt
    v_low
    solution_opt=x_opt(1:n)';
    current_objective_value=y
    Omega_v_low;
    flag
    Omega_k
end

ccp=0;
if flag==1
   for k=1:num 
     ccp=ccp+model.PComponents(k)*normcdf(-(f0+(a+B'*solution_opt')'*model.mu(k,:)')/sqrt((a+B'*solution_opt')'*model.Sigma(:,:,k)*(a+B'*solution_opt')));
   end    
    ccp  
    solution_opt=x_opt(1:n)'
    objective_value=sum(model.PComponents*(model.mu*solution_opt'))    
else
    flag=0;%infeasible
    objective_value=-Inf
end
endtime=cputime;  

spend_time=endtime-initialtime

sub_problems=iter_time
 