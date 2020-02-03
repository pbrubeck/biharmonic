function [P,D] = polydiffV(z,H)
n=size(H,2);
%z=z-H(1); c=H(2);
w0=(5+1i)/2; w0=0.5;
c=abs((5+1i)-w0); c=1;
z=z-w0;
P=(z(:)/c).^(0:n);
D=zeros(size(P));
D(:,2:end)=(1/c)*P(:,1:end-1)*diag(1:n);
end