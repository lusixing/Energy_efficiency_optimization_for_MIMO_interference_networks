% ---------------------------------------------
% Simulation code for varifing the convergence of the iterative SDP-WMMSE 
% algorithm for maximizing the energy efficience(EE) of multi-user MIMO
% interference networks under perfect and imperfect csi scenario. 
% Code written by Sixing Lu , email:562237037@qq.com
%-----------------------------------------------
clear
Nr = 2;  %Number of receive antennas
Nt = 2;  %Number of receive antennas
dk = 2;
Pd = 5;     %Maximum transmit power for each user 
Ps = 0.1;   % Static power comsumption
K = 2;      % number of users
var_noise=1;  %variance of additive gaussian noise

maxiter =10; %maximum iterations
nsim =10;   

v1 = zeros(Nt,dk,K); 
u1 = zeros(dk,Nr,K); 
m1 = zeros(dk,dk,K); 
w1 = zeros(dk,dk,K); 
v2 = zeros(Nt,dk,K); 
u2 = zeros(dk,Nr,K);
m2 = zeros(dk,dk,K); 
w2 = zeros(dk,dk,K); 

%----------imperfect channel-------------------
rho1 = 0.4;%receive correlation factor
for i=1:Nr
    for j=1:Nr
        R_r(i,j) = rho1^abs(i-j);
    end
end

rho2 = 0.2;%transmit correlation factor
for i=1:Nt
    for j=1:Nt
        R_t(i,j) = rho2^abs(i-j);
    end
end

var_CHN_error =0.1; %chanel estimation error rate

EE_set1=zeros(1,maxiter); EE_set2=zeros(1,maxiter); % store the network energy efficiency
R_set1=zeros(1,maxiter);  R_set2=zeros(1,maxiter);  % store the network sum rate
E_set1= zeros(1,maxiter); E_set2=zeros(1,maxiter);  % store the total power comsumption

for n=1:nsim
    h_real=zeros(Nr,Nt,K,K);
    h = h_real;%estimtion
    for l=1:K
        for m=1:K
            h_real(:, :, l, m) =R_r^0.5*(abs(randn(Nr, Nt))+randn(Nr, Nt)*1i)*R_t^0.5 ; %real channel
            h_real(:, :, l, m) =h_real(:, :, l, m)/norm(h_real(:, :, l, m) ,'fro');    %channel normalization
            Ew = R_r^0.5*var_CHN_error^0.5*(randn(Nr, Nt)+randn(Nr, Nt)*1i)*R_t^0.5;  
            h(:, :, l, m) = h_real(:, :, l, m) - Ew;                                  %channel estimation
            h(:, :, l, m) =h(:, :, l, m)/norm(h(:, :, l, m) ,'fro');                 
        end
    end
    
    v1=sqrt(Pd/(Nt*dk))*randn(Nt,dk,K);  
    for k=1:K
        J_k=var_noise*eye(Nr);
        for j=1:K
            if j~=k
                J_k=J_k+h(:,:,k,j)*v1(:,:,j)*v1(:,:,j)'*h(:,:,k,j)';
            end
        end
        
        u1(:,:,k)=v1(:,:,k)'*h(:,:,k,k)'*(J_k+h(:,:,k,k)*v1(:,:,k)*v1(:,:,k)'*h(:,:,k,k)')^-1;
        m1(:,:,k)=eye(dk)-u1(:,:,k)*h(:,:,k,k)*v1(:,:,k);
        w1(:,:,k)=m1(:,:,k)^-1;
    end
    v2=v1;
    u2=u1;
    m2=m1;
    w2=w1;
    
    %-------compute the initial EE for all users------
    [~,~,EE]=status(h,v1,var_noise,Ps);
    %-------------------------------------------------
    [lambda1]=min(EE);
    lambda2=lambda1;
    
    iter=1;
    while(1)
        %--------Run the BCD-WMMSE Algotithm -----------------------------
       %%  non-robust situation
        [v1] =SDP_rbst(lambda1,h,u1,v1,w1,var_noise,Pd,Ps,0,0,0);        %update the covariance matrix of each user using SDP algorithm
        [~,~,EE1] =status(h,v1,var_noise,Ps);                            %compute the estimated EE
        [lambda1] =min(EE1) ;                                          
        [R1,E1,~] =status(h_real,v1,var_noise,Ps);                       %compute the actual EE and SE
        
        for k=1:K                                                                                  % WMMSE parameter update
            J_k1=var_noise*eye(Nr);                                                                 % 
            for j=1:K                                                                                
                if j~=k                                                                              
                    J_k1=J_k1+h(:,:,k,j)*v1(:,:,j)*v1(:,:,j)'*h(:,:,k,j)';                           
                end                                                                                  
            end                                                                                      
            u1(:,:,k)=v1(:,:,k)'*h(:,:,k,k)'*(J_k1+h(:,:,k,k)*v1(:,:,k)*v1(:,:,k)'*h(:,:,k,k)')^-1;
            m1(:,:,k)=eye(dk)-u1(:,:,k)*h(:,:,k,k)*v1(:,:,k);
            w1(:,:,k)=m1(:,:,k)^-1;
        end
        
       %%  robust situation
        [v2] = SDP_rbst(lambda2,h,u2,v2,w2,var_noise,Pd,Ps,var_CHN_error,R_t,R_r);    %same as above, but take channel estimation error
        [~,~,EE2] =status(h,v2,var_noise,Ps);                                         % into consideration
        [lambda2] =min(EE2) ;                                                          %
        [R2,E2,~] =status(h_real,v2,var_noise,Ps);                                       %
        
        for k=1:K                                                                                       %WMMSE parameter update
            J_k2=var_noise*eye(Nr);                                                                      %
            for j=1:K                                                                                      %
                if j~=k                                                                                    %
                    J_k2=J_k2+h(:,:,k,j)*v2(:,:,j)*v2(:,:,j)'*h(:,:,k,j)';                                 %
                end                                                                                        
            end                                                                                            
            N_k=0;
            for j=1:K
                N_k=N_k+var_CHN_error*trace(R_t*v2(:,:,j)*v2(:,:,j)')*R_r;
            end
            
            u2(:,:,k)=v2(:,:,k)'*h(:,:,k,k)'*(J_k2+h(:,:,k,k)*v2(:,:,k)*v2(:,:,k)'*h(:,:,k,k)'+N_k)^-1;
            m2(:,:,k)=eye(dk)-u2(:,:,k)*h(:,:,k,k)*v2(:,:,k);
            w2(:,:,k)=m2(:,:,k)^-1;
        end      
        %---------------------------------------------
        if  iter>=maxiter
            break;
        end
        iter=iter+1;
        EE_set1(iter)=EE_set1(iter)+ sum(R1)/sum(E1);
        EE_set2(iter)=EE_set2(iter)+ sum(R2)/sum(E2);
        
        R_set1(iter)=R_set1(iter)+sum(R1);
        R_set2(iter)=R_set2(iter)+sum(R2);
        
        E_set1(iter)=E_set1(iter)+sum(E1);
        E_set2(iter)=E_set2(iter)+sum(E2);
        
    end
    n
end

EE_set1=EE_set1/(nsim);    %averaved the simulation results
EE_set2=EE_set2/(nsim);    
R_set1=R_set1/(nsim);      
R_set2=R_set2/(nsim);      
E_set1=E_set1/(nsim);      
E_set2=E_set2/(nsim);      

figure(1)
clf
plot (EE_set1, '-ro')
hold on
plot (EE_set2, '-bs')
legend('non-Robust ', 'Robust')
xlabel('Iteration')
ylabel('Energy Efficiency(bit/Hz/J)')

figure(2)
clf
plot (R_set1, '-ro')
hold on
plot (R_set2, '-bs')
legend('non-Robust ', 'Robust')
xlabel('Iteration')
ylabel('Sum Rate(bit/Hz)')

figure(3)
clf
plot (E_set1, '-ro')
hold on
plot (E_set2, '-bs')
legend('non-Robust ', 'Robust')
xlabel('Iteration')
ylabel('Energy comsumption (J/Hz)')