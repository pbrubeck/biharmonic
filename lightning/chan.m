%close all;

L = 3;
w = [1i+L,1i,1i-L,inf,-L,0,-1i,-1i+L,inf];


s1 = 1/5;  
% i = side (- or +)
% j = fun (f or g)
% k = diff (f or f')

a=[0,0,2,-4/3; 1/3,1/2,0,-1/6]; % cubics in y

deg=size(a,2)-1;
coef=zeros(deg+1,2,2,2);
for i=1:2
    coef(:,i,1,1)=[0;0.5i*a(i,3);0.75*a(i,4);0]; % f
    coef(:,i,2,1)=[1i*a(i,1);a(i,2);-0.5i*a(i,3);-0.25*a(i,4)]; % g
end
coef=reshape(coef,size(coef,1),4,2);
coef(1:end-1,:,2)=spdiags((1:deg)',0,deg,deg)*coef(2:end,:,1);    
coef=reshape(coef,[],2,2,2);
f=@(i,j,k,z) polyval(flipud(coef(:,i,j,k)),z); 
pfun=@(i,z) (1+(2*i-3)*tanh(s1*z))/2;
qfun=@(i,z) s1*((2*i-3)*(sech(s1*z)).^2)/2;
ff = @(i,z) pfun(i,z).*(conj(z).*f(i,1,1,z) + f(i,2,1,z));
fu = @(i,z) qfun(i,z).*(conj(z).*f(i,1,1,z) + f(i,2,1,z))+...
            pfun(i,z).*(conj(z).*f(i,1,2,z) + f(i,2,2,z))+...
            -conj(pfun(i,z).*f(i,1,1,z));
fp = @(i,z) pfun(i,z).*f(i,1,2,z)+qfun(i,z).*f(i,1,1,z);

f0 = @(z) -1i*(ff(1,z)+ff(2,z));        
u0 = @(z) conj(fu(1,z)+fu(2,z));
p0 = @(z) 4*conj(fp(1,z)+fp(2,z));

ubc=@(z) -u0(z);
%ubc=@(z) nan -1i*real(f0(z)) + 1i*(2/3)*(imag(z)>0.5);

g=cell(size(w));
[g{:}] = deal(ubc);
figure(1);
[u, maxerr, f, ZLS] = stokes(w,g,'tol',1E-10,'rel');

%return
nplt = 31;
H = 9;
x = linspace(-H,H,H*nplt);
y = linspace(-1,1,nplt);
[xx,yy]=ndgrid(x,y);
Z = xx+1i*yy;
Z(xx<0 & yy<0) = nan;

uz = u(Z)   + u0(Z);
fz = f(0,Z) + f0(Z);

figure(2);
%quiver(real(Z),imag(Z),real(uz),imag(uz)); 
contour(real(Z),imag(Z),real(fz),linspace(0,2/3,21)); 
%contour(real(Z),imag(Z),real(fz),21); 
%surf(real(Z),imag(Z),real(fz));
%surf(real(Z),imag(Z),abs(uz));
hold on;
plot(ZLS,'.k');
hold off; xlim([-H,H]); ylim([-5,5]);
dst = max(real(ZLS))-min(real(ZLS));
title(sprintf("Channel length = %E",dst));