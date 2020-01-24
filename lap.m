function [u,dofs,r] = lap(n,z,un,bctype,bcdata,w,pol,hol,mass)
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
    num=Dwp(id);
    %num=ones(numel(pol),1);
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

F=zeros(size(R0));
L=zeros(size(R0,1),numel(hol));

i=(bctype(:)==0);
F(i,:)=R0(i,:);
L(i,:)=log(z(i)-hol(:).');

i=(bctype(:)==1);
F(i,:)=spdiags(conj(un(i)))*R1(i,:);
L(i,:)=spdiags(conj(un(i)))*(1./(z(i)-hol(:).'));

A=[real(F),-imag(F),real(L)];
b=bcdata(:);

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
    W=spdiag(wt);
end
else
    W=spdiag(sqrt(mass));
end
% Left preconditioning
A=W*A;
r=W*b;
x=zeros(size(A,2),1);

% Right preconditioning
a=max(abs(A),[],1);
ja=a>eps;
P=spdiag(1./a(ja));
A=A(:,ja)*P;

z=A\r;
x(ja)=P*z;

r=(r-A*z)/norm(r);

dofs=length(x);
nb=size(R0,2);
h=x(2*nb+1:end);
x=reshape(x(1:2*nb),[],2);
f=x(:,1)+1i*x(:,2);
u=@eval_rat;
figure(32); semilogy(1:nb,abs(f));

function [f0,f1]=eval_rat(s)
nh=size(H,1);
nt=numel(pol);
ns=numel(hol);

T=zeros(numel(s),nb,2);
[T(:,1:nh,1),T(:,1:nh,2)]=polydiffA(s,H);
T(:,1+nh:nt+nh,1)=(1./(s(:)-pol(:).'))*NUM;
T(:,1+nh:nt+nh,2)=-T(:,1+nh:nt+nh,1)./(s(:)-pol(:).');

S=zeros(numel(s),ns,2);
S(:,:,1)=log(s(:)-hol(:).');
S(:,:,2)=1./(s(:)-hol(:).');

f0=reshape(T(:,:,1)*f+S(:,:,1)*h,size(s));
f1=reshape(T(:,:,2)*f+S(:,:,2)*h,size(s));
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