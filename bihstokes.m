function [psi,dofs,r] = bihstokes(n,z,un,bctype,bcdata,w,pol,hol,mass,ifp0)
% Solves the biharmonic equation with the lightning method
if(nargin<6), w=[]; end
if(nargin<7), pol=[]; end
if(nargin<8), hol=[]; end
if(nargin<9), mass=[]; end
if(nargin<10), ifp0=false; end
iflog=true;


% simple poles
[~,jj]=min(abs(w(:)-pol(:).'),[],1);
alpha=w(jj);
beta=1./(pol-alpha);
alpha=reshape(alpha,1,[]);
beta=reshape(beta,1,[]);

num1=min(1,beta);
num0=1./max(1,beta);
D=(z(:)-alpha).*beta-1;
C0=num0./D;
C1=-num1./(D.^2);

% polynomial basis, discretely orthogonal
if(n>=0)
    [~,H]=polyfitA(z,z,n);
    [V0,V1]=polydiffA(z,H);
else
    H=zeros(0,1);
    V0=zeros(numel(z),0); V1=V0;
end

% polynomials in holes, discretely orthogonal
RH=1./(z(:)./hol(:).'-1); 
nn=ceil(sqrt(numel(pol)/numel(w)));

HH=zeros(nn+1,nn,numel(hol));
HV0=zeros(numel(z),size(HH,1),numel(hol));
HV1=zeros(numel(z),size(HH,1),numel(hol));
for e=1:numel(hol)
    [~,HH(:,:,e)]=polyfitA(RH(:,e),RH(:,e),nn);
    [HV0(:,:,e),HV1(:,:,e)]=polydiffA(RH(:,e),HH(:,:,e));
    HV1(:,:,e)=spdiag(-RH(:,e).^2/hol(e))*HV1(:,:,e);
end
R0=[V0, reshape(HV0(:,2:end,:),numel(z),[]), C0]; 
R1=[V1, reshape(HV1(:,2:end,:),numel(z),[]), C1];

% logarithmic basis
HZ=z(:)-hol(:).';
if(~iflog), HZ=zeros(numel(z),0); end
XZ=real(HZ); YZ=imag(HZ);
QZ=abs(HZ).^2;
LZ=log(QZ)/2;
RZ=1./HZ;
L0=[LZ, XZ.*LZ, YZ.*LZ, QZ.*LZ];
L1=1i*[RZ, XZ.*RZ+LZ, YZ.*RZ-1i*LZ, QZ.*RZ.*(2*LZ+1)];
ZH=zeros(size(L0));

G0=R0;
F0=spdiag(conj(z))*R0;
G1=spdiag(un)*R1;
i=bctype==-1;
G1(i,:)=0;
F1=spdiag(conj(z).*un)*R1+spdiag(conj(un))*R0;
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

if(isempty(mass))
if(isempty(w))
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
    wt=sqrt(mass);
    W=spdiag([wt;wt]);
end

% Left preconditioning
A=W*A;
r=W*b;

%ifp0=false; %ifp0=true;
nsp=0; GB=0;
if(ifp0)
    w0=w(end); % Some arbitrary point
    [~,i]=min(abs(z(:)-w0));
    nsp=5-any(bctype==0);
    C=zeros(nsp,size(A,2));
    
    G1=R0(i,:);  F1=z(i)'*R0(i,:); % set psi(w0)=0
    C(end,:)=[imag(G1),real(G1),imag(F1),real(F1), imag(L0(i,:))];
    
    F1=mass'*R1; G1=zeros(size(F1)); % set mean(p)=0
    C(end-1,:)=[real(G1),-imag(G1),real(F1),-imag(F1), imag(L0(i,:))];
    
    G1=R1(i,:);  F1=zeros(size(G1)); % set g'(w0)=0
    C(1,:)=[real(G1),-imag(G1),real(F1),-imag(F1), imag(L0(i,:))];
    C(2,:)=[imag(G1), real(G1),imag(F1), real(F1), imag(L0(i,:))];
    
    G1=R0(i,:);  F1=zeros(size(G1)); % set real(g(w0))=0
    C(3,:)=[real(G1),-imag(G1),real(F1),-imag(F1), imag(L0(i,:))];

    [GB,rd]=rref(C); % Choose pivots
    kd=setdiff(1:size(C,2),rd);
    GB = -GB(:,kd);
    A = A(:,kd)+A(:,rd)*GB;
end

% Right preconditioning
a=sqrt(sum(A.^2,1));
%a=max(abs(A),[],1);
%a(:)=1;

ja=find(a>0); 
P=spdiag(1./a(ja));
A=A(:,ja)*P;

y=A\r;
r=(r-A*y)/norm(r);
dofs=length(y);

x=zeros(size(A,2),1);
x(ja)=P*y;
if(nsp>0)
    xx=zeros(numel(x)+size(GB,1),1);
    xx(kd)=x; xx(rd)=GB*x; x=xx;
    C*x
end

nb=size(R0,2);
h=x(4*nb+1:end);
x=reshape(x(1:4*nb),[],4);
g=x(:,1)+1i*x(:,2);
f=x(:,3)+1i*x(:,4);

psi=@goursat;
return

if(ifaaa)
    atol=1E-5;
    [af0,af1,ag0,ag1] = aaa_goursat(z,R0*f,R1*f,R0*g,R1*g,atol);
    psi=@(z) eval_goursat(af0,af1,ag0,ag1,z);
end

%figure(32); semilogy(1:nb,abs(f),1:nb,abs(g));


function [psi,u,p]=goursat(s)
nh=size(H,1);
nt=numel(pol);
ns=numel(hol);
RS=1./(s(:)./hol(:).'-1);
T=zeros(numel(s),nb,2);
if(nh>0)
    [T(:,1:nh,1),T(:,1:nh,2)]=polydiffA(s,H);
end
for j=1:size(HH,3)
    j1=1+nh+(j-1)*(size(HH,1)-1); j2=nh+j*(size(HH,1)-1);
    [SH0,SH1]=polydiffA(RS(:,j),HH(:,:,j));
    T(:,j1:j2,1)=SH0(:,2:end);
    T(:,j1:j2,2)=spdiag(-RS(:,j).^2/hol(j))*SH1(:,2:end);
end
nh=size(H,1)+(size(HH,1)-1)*size(HH,3);
T(:,1+nh:nt+nh,1)=num0./((s(:)-alpha).*beta-1);
T(:,1+nh:nt+nh,2)=-num1.*((s(:)-alpha).*beta-1).^-2;

L=zeros(numel(s),length(h),2);
HS=s(:)-hol(:).'; 
if(~iflog), HS=zeros(numel(s),0); end
XS=real(HS); YS=imag(HS);
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



function [af0,af1,ag0,ag1] = aaa_goursat(Z,f0,f1,g0,g1,tol)
if(nargin<6), tol=1E-6; end
N = numel(Z);
[af0,polf0] = aaa(f0,Z,'mmax',N/2,'tol',tol,'lawson',0,'cleanup',0);
[af1,polf1] = aaa(f1,Z,'mmax',N/2,'tol',tol,'lawson',0,'cleanup',0);
[ag0,polg0] = aaa(g0,Z,'mmax',N/2,'tol',tol,'lawson',0,'cleanup',0);
[ag1,polg1] = aaa(g1,Z,'mmax',N/2,'tol',tol,'lawson',0,'cleanup',0);

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



function [psi,u,p]=eval_goursat(af0,af1,ag0,ag1,s)
L=zeros(numel(s),length(h),2);
HS=s(:)-hol(:).'; 
if(~iflog), HS=zeros(numel(s),0); end
XS=real(HS); YS=imag(HS);
QS=abs(HS).^2;
LS=log(QS)/2;
RS=1./HS;
L(:,:,1)=[LS, XS.*LS, YS.*LS, QS.*LS];
L(:,:,2)=1i*[RS, XS.*RS+LS, YS.*RS-1i*LS, QS.*RS.*(2*LS+1)];

f0=af0(s(:));  f1=af1(s(:));
g0=ag0(s(:));  g1=ag1(s(:)); 
h0=L(:,:,1)*h; h1=L(:,:,2)*h;

psi=reshape(-1i*(g0+conj(s(:)).*f0+1i*h0),size(s));
u=reshape(conj(g1+conj(s(:)).*f1-conj(f0)+h1),size(s));
p=reshape(conj(4*f1),size(s));
end




end