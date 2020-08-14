% ---------------------------------------------
% Simulation code for varifing the convergence of the iterative SCA
% algorithm for maximizing the energy efficience(EE) of multi-user MIMO
% interference networks under perfect and imperfect csi scenario. 
% Code written by Sixing Lu , email:562237037@qq.com
%-----------------------------------------------
clear
Nr =2;      %Number of receive antennas
Nt =2;      %Number of receive antennas
dk =2;
Pd =5;        %Maximum transmit power for each user 
Ps =0.1;        % Static power comsumption
K =2;          % number of users
var_noise =1;   %variance of additive gaussian noise

max_iter =20;   
nsim =10;       

%----------imperfect channel-------------------
rho1 = 0.4;%receive correlation
for i=1:Nr
    for j=1:Nr
        R_r(i,j) = rho1^abs(i-j);
    end
end

rho2 = 0.2;%transmit correlation
for i=1:Nt
    for j=1:Nt
        R_t(i,j) = rho2^abs(i-j);
    end
end

var_CHN_error=0.1;  %chanel estimation error rate

Rate_set1=zeros(1,max_iter);EE_set1=zeros(1,max_iter);
E_set1=zeros(1,max_iter);Rate_set2=zeros(1,max_iter);
EE_set2=zeros(1,max_iter);E_set2=zeros(1,max_iter);

for s=1:nsim
    Q_int=zeros(Nt,Nt,K);
    for i=1:K                                      %initialize covariance of precoding matrix with tr(Q_k)<=Pd
        tmp=sqrt(Pd/(Nt*dk))*randn(Nt,dk) ;        %
        Q_int(:,:,i)=tmp*tmp';                     %
    end                                              %
    Q1=Q_int;
    Q2=Q_int;
    h_real=zeros(Nr,Nt,K,K);
    h_esti = h_real;%estimtion
    for n=1:K
        for m=1:K
            h_real(:, :, n, m) =  R_r^0.5*(randn(Nr, Nt)+randn(Nr, Nt)*1i)*R_t^0.5;
            h_real(:, :, n, m) =h_real(:, :, n, m)/norm(h_real(:, :, n, m) ,'fro');   
            Ew = R_r^0.5*sqrt(var_CHN_error)*(randn(Nr, Nt)+randn(Nr, Nt)*1i)*R_t^0.5;
            h_esti(:, :, n, m) = h_real(:, :, n, m) - Ew;          %channel estimation
            h_esti(:, :, n, m) =h_esti(:, :, n, m)/norm(h_esti(:, :, n, m) ,'fro');   
        end
    end
    
    [R1,E1,EE1] =status(h_esti,Q_int,var_noise,Ps); %initialize
    lambda1 =min(EE1);                            %
    lambda2 =lambda1;                             %
    
    EE_set_tmp1 =[];EE_set_tmp2 =[];               
    E_set_tmp1 =[]; E_set_tmp2 =[];               
    Rate_set_tmp1 =[]; Rate_set_tmp2 =[];
    
    iter=1;
    
    while (iter<=max_iter)
        [Q1] = SCA_rbst(h_esti,lambda1,Q1,var_noise,Pd,Ps,0,R_r,R_t);  % update precoding covariance matrix using sca algorithm
        [~,~,EE1] =status(h_esti,Q1,var_noise,Ps);        %estimate the new EE 
        lambda1 =min(EE1);                              %
        [R1,E1,~] =status(h_real,Q1,var_noise,Ps);      %compute the real rate and EE
        EE_set_tmp1 =[EE_set_tmp1 sum(R1)/sum(E1)];   %store the results
        Rate_set_tmp1 =[Rate_set_tmp1 sum(R1)];       %
        E_set_tmp1 =[E_set_tmp1 sum(E1)];             %
        
        [Q2] =SCA_rbst(h_esti,lambda2,Q2,var_noise,Pd,Ps,var_CHN_error,R_r,R_t);   %take channel uncertinty into consideration
        [~,~,EE2] =status(h_esti,Q2,var_noise,Ps);
        lambda2 =min(EE2);
        [R2,E2,~] =status(h_real,Q2,var_noise,Ps);
        EE_set_tmp2 =[EE_set_tmp2 sum(R2)/sum(E2)];
        Rate_set_tmp2 =[Rate_set_tmp2 sum(R2)];
        E_set_tmp2 =[E_set_tmp2 sum(E2)];
        
        iter =iter+1;
    end
    Rate_set1=Rate_set1+Rate_set_tmp1;
    EE_set1=EE_set1+EE_set_tmp1;
    E_set1=E_set1+E_set_tmp1;
    Rate_set2=Rate_set2+Rate_set_tmp2;
    EE_set2=EE_set2+EE_set_tmp2;
    E_set2=E_set2+E_set_tmp2;
    
    s
end

Rate_set1=Rate_set1/(nsim);   
EE_set1=EE_set1/(nsim);
E_set1=E_set1/(nsim);
Rate_set2=Rate_set2/(nsim);
EE_set2=EE_set2/(nsim);
E_set2=E_set2/(nsim);
save

figure(1)
clf
set(gcf, 'color', 'white')
plot( EE_set1, '-ro')
hold on
plot( EE_set2, '-bs')
xlabel('Iteration')
ylabel('Energy Efficiency(bit/Hz/J)')
legend('SCA-non-Robust ', 'SCA-Robust')

figure(2)
clf
set(gcf, 'color', 'white')
plot( Rate_set1, '-ro')
hold on
plot( Rate_set2, '-bs')
xlabel('Iteration')
ylabel('Sumrate(bit/Hz)')
legend('SCA-non-Robust ', 'SCA-Robust')

figure(3)
clf
set(gcf, 'color', 'white')
plot( E_set1, '-ro')
hold on
plot( E_set2, '-bs')
xlabel('Iteration')
ylabel('Total Energy coumsuption (J)')
legend('SCA-non-Robust ', 'SCA-Robust')

grid on
