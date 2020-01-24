function [P,D] = polydiffA(z,H)
M=numel(z); n=size(H,2); 
P=ones(M,n+1);
D=zeros(M,n+1);
for k=1:n
    p=z(:).*P(:,k);
    p=p-P(:,1:k)*H(1:k,k);
    P(:,k+1)=p/H(k+1,k);
    
    d=z(:).*D(:,k)+P(:,k);
    d=d-D(:,1:k)*H(1:k,k);
    D(:,k+1)=d/H(k+1,k);
end
end