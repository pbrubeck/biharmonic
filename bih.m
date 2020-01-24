function [u,res,dofs] = bih(n,z,un,u0,u1,pol,w)
% Solves the biharmonic equation with the lightning method
if(nargin<5), pol=[]; end
if(nargin<6), w=[]; end

if(isempty(w))
    num=ones(numel(pol),1);
else
    Dwp=w(:)-pol(:).';
    [~,ii]=min(abs(Dwp),[],1);
    id=sub2ind(size(Dwp), ii, 1:size(Dwp,2));
    num=-Dwp(id).^2/sqrt(numel(z));
end
%num=ones(numel(pol),1);
NUM=spdiag(num);
D=z(:)-pol(:).';
C0=(1./D)*NUM;
C1=-C0./D;
[~,H]=polyfitA(z,z,n);
%H(1)=mean(z); H(2)=max(abs(z(:)-H(1)));
[V0,V1]=polydiffA(z,H);

A0=[V0, C0];
A1=spdiag(conj(un))*[V1, C1];
B0=spdiag(conj(z))*A0;
B1=spdiag(conj(z))*A1+spdiag(un)*A0;
A=[real(A0),imag(A0),real(B0),imag(B0);
   real(A1),imag(A1),real(B1),imag(B1)];
b=[u0; real(conj(un).*u1(:))];

if(isempty(D))
    W=speye(size(A,1))/norm(b);
else
    Dzw=z(:)-w(:).';
    wt=min(abs(Dzw),[],2);
    W=spdiag([wt;wt]);
    W=W/norm(W*b);
end
x=(W*A)\(W*b);

res=norm(W*(A*x-b));
dofs=length(x);

m=numel(pol);
nn=n+m;
f=x(1:nn+1)-1i*x(nn+2:2*nn+2);
g=x(2*nn+3:3*nn+3)-1i*x(3*nn+4:end);
u=@goursat;

function [u,du]=goursat(s)
Dsp=s(:)-pol(:).';
T=zeros(numel(s),m+n+1,2);
[T(:,1:size(H,1),1),T(:,1:size(H,1),2)]=polydiffA(s,H);
T(:,1+size(H,1):end,1)=(1./Dsp)*NUM;
T(:,1+size(H,1):end,2)=-T(:,1+size(H,1):end,1)./Dsp;
u =reshape(real(T(:,:,1)*f+conj(s(:)).*(T(:,:,1)*g)),size(s));
du=reshape(T(:,:,2)*f+conj(s(:)).*(T(:,:,2)*g)+conj(T(:,:,1)*g),size(s));
end
end