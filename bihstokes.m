function [psi,dofs,r] = bihstokes(n,z,un,bctype,bcdata,w,pol,hol,mass)
% Solves the biharmonic equation with the lightning method
if(nargin<6), w=[]; end
if(nargin<7), pol=[]; end
if(nargin<8), hol=[]; end
if(nargin<9), mass=[]; end

% singular basis, normalized to unity in the nearest corner
if(isempty(w))
    num=ones(numel(pol),1);
else
    Dwp=w(:)-pol(:).';
    [~,ii]=min(abs(Dwp),[],1);
    id=sub2ind(size(Dwp), ii, 1:size(Dwp,2));
    num=1i*Dwp(id);
end
NUM=spdiag(num);
D=z(:)-pol(:).';
C0=(1./D)*NUM;
C1=-C0./D;

% polynomial basis, discretely orthogonal
[~,H]=polyfitA(z,z,n);
[V0,V1]=polydiffA(z,H);
R0=[V0, C0]; 
R1=[V1, C1];

% logarithmic basis
HZ=z(:)-hol(:).'; XZ=real(HZ); YZ=imag(HZ);
QZ=abs(HZ).^2;
LZ=log(QZ)/2;
RZ=1./HZ;
L0=[LZ, XZ.*LZ, YZ.*LZ, QZ.*LZ];
L1=1i*[RZ, XZ.*RZ+LZ, YZ.*RZ-1i*LZ, QZ.*RZ.*(2*LZ+1)];
ZH=zeros(size(L0));


G0=R0;
F0=spdiag(conj(z))*G0;
G1=spdiag(un)*R1;
i=bctype==-1;
G1(i,:)=0;
F1=spdiag(conj(z))*G1+spdiag(conj(un))*G0;
F1(i,:)=F1(i,:);

AR=zeros(numel(z),4*size(R0,2)+size(L0,2));
AI=zeros(numel(z),4*size(R0,2)+size(L0,2));
i=(bctype(:)<=0);

L1(i,:)=spdiag(un(i))*L1(i,:);

AR(i,:)=[imag(G1(i,:)), real(G1(i,:)),imag(F1(i,:)), real(F1(i,:)), imag(L1(i,:))];
AI(i,:)=[imag(G0(i,:)), real(G0(i,:)),imag(F0(i,:)), real(F0(i,:)), real(L0(i,:))];

G=R1;
F=zeros(size(R0,1),size(R0,2),2);
F(:,:,1)=spdiag(conj(z))*R1-R0;
F(:,:,2)=spdiag(conj(z))*R1+R0;

i=(bctype(:)==2);
F(i,:,1)=F(i,:,1)+4*R0(i,:); % Add presssure integral (4*f^*)
F(i,:,2)=F(i,:,2)-4*R0(i,:);

i=(bctype(:)==3);
F(i,:,1)=F(i,:,1)+2*R0(i,:); % Free boundary
F(i,:,2)=F(i,:,2)-2*R0(i,:);

i=(bctype(:)>0);
AR(i,:)=[real(G(i,:)),-imag(G(i,:)),real(F(i,:,1)),-imag(F(i,:,1)), real(L1(i,:))];
AI(i,:)=[imag(G(i,:)), real(G(i,:)),imag(F(i,:,2)), real(F(i,:,2)), imag(L1(i,:))];


% Obstacle in differential form
i=(bctype(:)==4);
G=zeros(size(R0,1),size(R0,2),2);
F=zeros(size(R0,1),size(R0,2),2);
G(:,:,1)=spdiag(un)*R1;
F(:,:,1)=spdiag(un.*conj(z))*R1;
F(:,:,2)=spdiag(1i*conj(un))*R0;
AR(i,:)=[real(G(i,:,1)),-imag(G(i,:,1)),real(F(i,:,1)),-imag(F(i,:,1)), ZH(i,:)];
AI(i,:)=[imag(G(i,:,2)), real(G(i,:,2)),imag(F(i,:,2)), real(F(i,:,2)), ZH(i,:)];

A=[AR; AI];
b=[real(bcdata(:)); -imag(bcdata(:))];

%D=[];
%mass=[];
if(isempty(mass))
if(isempty(D))
    W=speye(size(A,1))/norm(b);
else
    Dzw=z(:)-w(:).';
    [wt,id]=min(abs(Dzw),[],2);
    for k=1:size(Dzw,2)
        dd=wt(id==k);
        wt(id==k)=dd./max(dd);
    end
    wt=sqrt(wt);
    W=spdiag([wt;wt]);
end
else
    mass=sqrt(mass);
    W=spdiag([mass;mass]);
end

% Left preconditioning
A=W*A;
r=W*b; 
x=zeros(size(A,2),1);

% Right preconditioning
a=max(abs(A),[],1);
ja=find(a>0);

%a(:)=1;
P=spdiag(1./a(ja));
A=A(:,ja)*P;
z=A\r;
x(ja)=P*z;

r=(r-A*z)/norm(r);
dofs=length(z);
nb=size(R0,2);
h=x(4*nb+1:end);
x=reshape(x(1:4*nb),[],4);
g=x(:,1)+1i*x(:,2);
f=x(:,3)+1i*x(:,4);


% [f,g]
% h
% f(:)=0; g(:)=0; h(:)=0;
% g(1)=7.50936;
% f(2)=-6.50936;
% h(1)=8.44817;
% h(4)=4.57055;
% g(1)=1i*g(1);
% f(2)=4i/2.53*f(2);
% 
% [f,g]
% h

psi=@goursat;
figure(32); semilogy(1:nb,abs(f),1:nb,abs(g));

function [psi,u,p]=goursat(s)
nh=size(H,1);
nr=numel(pol);
T=zeros(numel(s),length(f),2);

j1=1:nh;
j2=(j1(end)+1):(j1(end)+nr);
[T(:,j1,1),T(:,j1,2)]=polydiffA(s,H);
T(:,j2,1)=(1./(s(:)-pol(:).'))*NUM;
T(:,j2,2)=-T(:,j2,1)./(s(:)-pol(:).');

L=zeros(numel(s),length(h),2);
HS=s(:)-hol(:).'; XS=real(HS); YS=imag(HS);
QS=abs(HS).^2;
LS=log(QS)/2;
RS=1./HS;
L(:,:,1)=[LS, XS.*LS, YS.*LS, QS.*LS];
L(:,:,2)=1i*[RS, XS.*RS+LS, YS.*RS-1i*LS, QS.*RS.*(2*LS+1)];

f0=T(:,:,1)*f; f1=T(:,:,2)*f;
g0=T(:,:,1)*g; g1=T(:,:,2)*g;
h0=L(:,:,1)*h; h1=L(:,:,2)*h;

psi=reshape(-1i*(g0+conj(s(:)).*f0+1i*h0),size(s));
u=reshape(conj(g1+conj(s(:)).*f1-conj(f0)+h1),size(s));
p=reshape(conj(4*f1),size(s));
end
end