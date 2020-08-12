% ---------------------------------------------
% Semi-definite programming(SDP) part for the iterative SDP-WMMSE algorithm
% update the the precoding matrix for each user
% Code written by Sixing Lu , email:562237037@qq.com
% ---------------------------------------------

function [v_n] = SDP_rbst(lambda,h,u,v,w,p_noise,Pd,Ps,chn_err,Rt,Rr)
Nr=size(h,1);
Nt=size(h,2);
K=size(h,3);
dk=size(v,2);

H=zeros(Nr*K,Nt*K);
for i=1:K
    for j=1:K
        H((i-1)*Nt+1:i*Nr,(j-1)*Nt+1:j*Nt)=h(:,:,i,j);   
    end
end

%-----------robust------------------
if chn_err~=0
   R=zeros(Nt*K,Nt*K);
   for i=1:K
     R((i-1)*Nt+1:i*Nt,(i-1)*Nt+1:i*Nt)=Rt; 
   end
   R_c=chol(R);
else
   R_c=eye(Nt*K); 
end
%------------------------------------
cvx_begin sdp quiet
    variable tau
    variable v_n(Nt,dk,K) complex
    expression V
    V=v_n(:,:,1);                    % block-diag matrix
    for i=2:K                        % 
           V=blkdiag(V,v_n(:,:,i));  %
    end
    
    maximize(tau )
    subject to
     for k=1:K
          w(:,:,k)=w(:,:,k)-1i*imag(diag(diag(w(:,:,k))));    
         [w_c]=chol(w(:,:,k),'lower');                        %cholesky decomposition
         
         Theta=zeros(Nr,Nr*K);                                
         Xi=zeros(dk,dk*K);                                   
         Theta(:,((k-1)*Nr+1):(k*Nr))=eye(Nr);                
         Xi(:,((k-1)*dk+1):(k*dk))=eye(dk);                    
         
         T1=sqrt(chn_err*trace(w(:,:,k)*u(:,:,k)*Rr*u(:,:,k)'))*vec(R_c*V) ; % rubust term
         T2= lambda^0.5*vec(v_n(:,:,k));                                     
         T3=kron(eye(dk*K),w_c*u(:,:,k)*Theta*H)*vec(V)-vec(w_c*Xi);         
         
        [eye(dk*Nt) vec(v_n(:,:,k));vec(v_n(:,:,k))' Pd]'>=0                                                           
        sym_cvx([         eye(length(T1)+length(T2)+length(T3))         [T1;T2;T3] ;                                   % ensure the symmetry of the £¬                   
        [T1;T2;T3]'   -tau+log(det(w(:,:,k)))+dk-lambda*Ps-p_noise*trace( (w_c*u(:,:,k))*(w_c*u(:,:,k))' )] )   >= 0   % constrain matrix
     end                                                                                                               
     
cvx_end  
%-----------------------
return

