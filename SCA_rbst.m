% ---------------------------------------------
% Successive convex approximation (SCA) algorithm updating the percoding 
% covariance matrix of each user in the mimo interfrence networks, 
% with presence of channel estimation uncertainty.
% Code written by Sixing Lu, email:562237037@qq.com
% ---------------------------------------------
function [Q_opti] = SCA_rbst(h,lambda,Q_n,p_noise,Pd,Ps,var_ce,R_r,R_t)
Nr =size(h,1);
Nt =size(h,2);
K =size(h,3);
alpha =1;

cvx_begin quiet sdp
    variable Q_new(Nt,Nt,K) complex 
    variable tau
    maximize ( tau )
    subject to
      for k =1:K
          Jk=var_ce/p_noise*R_r*trace(R_t*Q_new(:,:,k));
          Jk_n=var_ce/p_noise*R_r*trace(R_t*Q_n(:,:,k));
          for j=1:K
              if j~=k
               Jk =Jk+ 1/p_noise*(h(:,:,k,j)*Q_new(:,:,j)*h(:,:,k,j)'+var_ce*R_r*trace(R_t*Q_new(:,:,j))  );
               Jk_n =Jk_n+ 1/p_noise*(h(:,:,k,j)*Q_n(:,:,j)*h(:,:,k,j)'+var_ce*R_r*trace(R_t*Q_n(:,:,j))  );
              end
          end
                 
        -log_det(eye(Nr)+1/p_noise*h(:,:,k,k)*Q_new(:,:,k)*h(:,:,k,k)'+Jk)+real( real(log(det(eye(Nr)+Jk_n)))+trace( (eye(Nr)+Jk_n)^-1*(Jk-Jk_n) )  +lambda*alpha*(trace(Q_new(:,:,k))+Ps ))  +tau <=0
        real( trace(Q_new(:,:,k)))<=Pd
        sym_cvx( Q_new(:,:,k))>=0
        %Q_new(:,:,k)>=0
      end   
 cvx_end  

Q_opti=Q_new;
end
