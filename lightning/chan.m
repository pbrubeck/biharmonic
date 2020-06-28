close all;

L = 3;
w = [1i+L,1i,1i-L,inf,-L,0,-1i,-1i+L,inf];

sfun = @tanh; % sigmoid function
s1 = 5;  s2 = 5;

f0 = @(z) -1i*(((z+1i)/2).*(1+sfun(z/s1))/2  + z.*(1-sfun(z/s2))/2 );
gfun=@(z) (imag(z)>0.5) - real(f0(z));
g=cell(size(w));
[g{:}] = deal(gfun);
figure(1);
[u, maxerr, f, ZLS] = lap(w,g,'tol',1E-10,'rel');

%return
nplt = 63;
H = 9;
x = linspace(-H,H,H*nplt);
y = linspace(-1,1,nplt);
[xx,yy]=ndgrid(x,y);
Z = xx+1i*yy;
Z(xx<0&yy<0) = nan;

fz = f(Z) + f0(Z);
figure(2);
contour(real(Z),imag(Z),real(fz),linspace(0,1,11),'b'); hold on;
contour(real(Z),imag(Z),imag(fz),7*H,'r'); 
plot(ZLS,'.k');
hold off; xlim([-H,H]); ylim([-5,5]);
dst = max(real(ZLS))-min(real(ZLS));
title(sprintf("Channel length = %E",dst));