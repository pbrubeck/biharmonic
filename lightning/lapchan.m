%close all;


w = [1i+3,1i,1i-3,inf,-3,0,-1i,-1i+3,inf];
g = [1,1,nan,nan,0,0,0,nan,nan];
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

fz = f(Z);
figure(2);
contour(real(Z),imag(Z),real(fz),linspace(0,1,11),'b'); hold on;
contour(real(Z),imag(Z),imag(fz),7*H,'r'); 
plot(ZLS,'.k');
hold off; xlim([-H,H]); ylim([-5,5]);
dst = max(real(ZLS))-min(real(ZLS));
title(sprintf("Channel length = %E",dst));