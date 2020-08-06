%close all;

L = 3;
w = [1i+L,1i,1i-L,inf,-L,0,-1i,-1i+L,inf];
g = [nan+1i/1.5;nan+1i/1.5;nan;nan;nan;nan;nan;nan;nan];

nw = length(w);
winf = find(isinf(w));

fbkg = @(z) 0*z;
ubkg = @(z) 0*z;
H = [1 -1 1 -1; 1 1 1 1; 0 1 -2 3; 0 1 2 3];
for k=winf   
    % Set up  a cubic
    kpp = mod(k-3,nw)+1;
    kp = mod(k-2,nw)+1;
    kn = mod(k,nw)+1;
    w1 = w(kp);
    w2 = w(kn);
    w0 = (w2+w1)/2;
    sk = w1-w(kpp);
    w2 = w2 - real((w2-w1)*conj(sk))/conj(sk); 
    a0 = (w2-w1)/2i;

    s = 1;
    coef = H\[imag(g(kpp)); imag(g(kn)); 0; 0];
    fbkg = @(z) fbkg(z) + imag(goursat_cubic(coef,s,(z-w0)/a0));
    ubkg = @(z) ubkg(z) + (1/a0)*velocity_cubic(coef,s,(z-w0)/a0);
end

npts = 65;
H = 2;
x=linspace(-H*L,H*L,H*npts);
y=linspace(-1,1,npts);
[xx,yy]=ndgrid(x,y);
zz = xx+1i*yy;
zz(yy<0 & xx<0)=nan;
surf(xx,yy,real(fbkg(zz)));
axis equal;


function [f]=goursat_cubic(c, s, z)
    z2 = imag(z).*z;
    z3 = 0.25*(3*conj(z)-z).*z.^2;
    f = 1i*c(1)+z.*c(2)+z2.*c(3)+z3.*c(4);
    f = f.*(1+tanh(s*z))/2;
end

function [u]=velocity_cubic(c, s, z)
    z2 = imag(z).*z;
    z3 = 0.25*(3*conj(z)-z).*z.^2;
    f = 1i*c(1)+z.*c(2)+z2.*c(3)+z3.*c(4);
	th = (1+tanh(s*z))/2;
    u1 = c(2) + c(3)*(2*z-conj(z))/2i + c(4)*(0.75*z.*(2*conj(z)-z));
    u2 =        c(3)*conj(z)/2i + c(4)*(0.75*conj(z).^2);
    u = (0.5*s*sech(s*z).^2).*f + th.*u1 - conj(th).*u2;
end