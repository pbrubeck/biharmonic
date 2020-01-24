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

% logarithmic basis
L0=log(z(:)-hol(:).');
L1=1./(z(:)-hol(:).');
ZL0=(z(:)-hol(:).').*log(z(:)-hol(:).');
ZL1=log(z(:)-hol(:).')+1;

R0=[V0, C0, L0, ZL0]; 
R1=[V1, C1, L1, ZL1];


G0=R0;
F0=spdiag(conj(z))*G0;
G1=spdiag(un)*R1;
i=bctype==-1;
G1(i,:)=0;
F1=spdiag(conj(z))*G1+spdiag(conj(un))*G0;
F1(i,:)=1i*F1(i,:);

AR=zeros(numel(z),4*size(R0,2));
AI=zeros(numel(z),4*size(R0,2));
i=(bctype(:)<=0);
AR(i,:)=[imag(G1(i,:)), real(G1(i,:)),imag(F1(i,:)), real(F1(i,:))];
AI(i,:)=[imag(G0(i,:)), real(G0(i,:)),imag(F0(i,:)), real(F0(i,:))];

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
AR(i,:)=[real(G(i,:)),-imag(G(i,:)),real(F(i,:,1)),-imag(F(i,:,1))];
AI(i,:)=[imag(G(i,:)), real(G(i,:)),imag(F(i,:,2)), real(F(i,:,2))];


% Obstacle in differential form
i=(bctype(:)==4);
G=zeros(size(R0,1),size(R0,2),2);
F=zeros(size(R0,1),size(R0,2),2);
G(:,:,1)=spdiag(un)*R1;
F(:,:,1)=spdiag(un.*conj(z))*R1;
F(:,:,2)=spdiag(1i*conj(un))*R0;
AR(i,:)=[real(G(i,:,1)),-imag(G(i,:,1)),real(F(i,:,1)),-imag(F(i,:,1))];
AI(i,:)=[imag(G(i,:,2)), real(G(i,:,2)),imag(F(i,:,2)), real(F(i,:,2))];

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
ja=find(a>eps);

%a(:)=1;
P=spdiag(1./a(ja));
A=A(:,ja)*P;

z=A\r;
%z=(A'*A)\(A'*r);
%[Q,R]=qr(A,0);z=R\(Q'*r);


x(ja)=P*z;

r=(r-A*z)/norm(r);
dofs=length(z);
nb=size(R0,2);

x=reshape(x,[],4);
g=x(:,1)+1i*x(:,2);
f=x(:,3)+1i*x(:,4);
psi=@goursat;
figure(32); semilogy(1:nb,abs(f),1:nb,abs(g));

function [psi,u,p]=goursat(s)
nh=size(H,1);
nr=numel(pol);
nl=numel(hol);
T=zeros(numel(s),length(f),2);

j1=1:nh;
j2=(j1(end)+1):(j1(end)+nr);
if(nl>0)
    j3=(j2(end)+1):(j2(end)+nl);
    j4=(j3(end)+1):(j3(end)+nl);
else
    j3=[]; j4=[];
end

[T(:,j1,1),T(:,j1,2)]=polydiffA(s,H);
T(:,j2,1)=(1./(s(:)-pol(:).'))*NUM;
T(:,j2,2)=-T(:,j2,1)./(s(:)-pol(:).');
T(:,j3,1)=log(s(:)-hol(:).');
T(:,j3,2)=1./(s(:)-hol(:).');
T(:,j4,1)=(s(:)-hol(:).').*log(s(:)-hol(:).');
T(:,j4,2)=log(s(:)-hol(:).')-1;

f0=T(:,:,1)*f; f1=T(:,:,2)*f;
g0=T(:,:,1)*g; g1=T(:,:,2)*g;

psi=reshape(-1i*(g0+conj(s(:)).*f0),size(s));
u=reshape(conj(g1+conj(s(:)).*f1-conj(f0)),size(s));
p=reshape(conj(4*f1),size(s));
% f0=R0*f; f1=R1*f; g0=R0*g; g1=R1*g;
% [aaapsi,aaau,aaap]=aaa_goursat(z(:),f0,f1,g0,g1);
% psi=reshape(aaapsi(s),size(s));
% u=reshape(aaau(s),size(s));
% p=reshape(aaap(s),size(s));
end

end



function [psi,u,p] = aaa_goursat(Z,f0,f1,g0,g1,tol)
if(nargin<6), tol=1E-6; end
N = numel(Z);
[aaaf0,polf0] = aaa(f0,Z,'mmax',N/2,'tol',0.1*tol,'lawson',0,'cleanup',0);
[aaaf1,polf1] = aaa(f1,Z,'mmax',N/2,'tol',0.1*tol,'lawson',0,'cleanup',0);
[aaag0,polg0] = aaa(g0,Z,'mmax',N/2,'tol',0.1*tol,'lawson',0,'cleanup',0);
[aaag1,polg1] = aaa(g1,Z,'mmax',N/2,'tol',0.1*tol,'lawson',0,'cleanup',0);

psi=@(z) imag(aaag0(z)+conj(z).*aaaf0(z));
u=@(z) conj(aaag1(z)+conj(z).*aaaf1(z)-conj(aaaf0(z)));
p=@(z) conj(4*aaaf1(z));
return
inpolygonc = @(z,w) inpolygon(real(z), ...
            imag(z),real(w),imag(w));  
polaaa = [polf0(:);polf1(:);polg0(:);polg1(:)];
        
if isempty(find(inpolygonc(polaaa,w),1)) % AAA successful
else                                     % AAA unsuccess.: pole in region
   %badpol = polaaa(find(inpolygonc(polaaa,ww)));% poles in polygon
    warning('GOURSAT AAA compression unsuccessful; returning uncompressed solution.')
end
end