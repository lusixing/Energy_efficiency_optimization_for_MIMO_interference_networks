%-------------------------------------------------------------
%compute the sum-rate energy efficience and power comsumption
%of the MIMO interference network
%-----------------------------------------------------
function [R,E,EE]=status(h,B,p_noise,Ps)
Nr=size(h,1);
Nt=size(h,2);
K=size(h,3);
Q=zeros(Nt,Nt,K);
for i=1:K
    Q(:,:,i)=B(:,:,i)*B(:,:,i)';
end
R=zeros(K,1);
E=zeros(K,1);
EE=zeros(K,1);

for k=1:K
  for j=1:K  
      Jk=eye(Nr)*p_noise;
      if j~=k
          Jk=Jk+ h(:,:,k,j)*Q(:,:,j)*h(:,:,k,j)';
      end    
  end
  R(k)=real(log(det(eye(Nr)+( h(:,:,k)*Q(:,:,k)*h(:,:,k)')*(Jk+p_noise*eye(Nr)  )^-1 )));
  E(k)=real(trace(Q(:,:,k)))+Ps;
  EE(k)=R(k)/E(k);
end