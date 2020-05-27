function [res] = moffatt(alpha)
if(nargin<1)
    alpha=pi/18;
end
ifprint=false;
prefix='moff';
lw='linewidth';
ms='markersize';
kmax=12;
ne=10; % number of eddies
nc=5;  % number of contours per eddy
np=128;

sa=sin(2*alpha)/(2*alpha);
g=@(x,y) [sin(x).*cosh(y)+sa*x; cos(x).*sinh(y)+sa*y];
x=[4.21; 2.26]; 
x=NewtonRaphson(g,x);
lambda=1+([1,1i]*x)/(2*alpha);
disp(lambda-1);

% Antisymmetric
A= cos((lambda-2)*alpha);
B=-cos(lambda*alpha);
es=angle(-lambda*(A+B));
r0=exp(es/imag(lambda));
rho=exp(pi/imag(lambda));
U=real((r0.^lambda)/(lambda*(A+B)));
rmax=r0*exp(-(es+pi*(0:nc*ne)/nc)/imag(lambda));

f=@(th) A*cos(lambda*th)+B*cos((lambda-2)*th);
ff=@(th) -lambda*A*sin(lambda*th)-(lambda-2)*B*sin((lambda-2)*th);
fex=@(r,th) U*real(((r/r0).^lambda).*f(th));
uex=@(r,th) conj((1i*U/r0)*exp(-1i*th).*(real((r/r0).^(lambda-1).*f(th)*lambda) + ...
                              -1i*real((r/r0).^(lambda-1).*ff(th))));

R=1;
w=[R*exp(1i*alpha*[-1;1]);0]; 
rad=inf(size(w));
rad(1)=R;
if(all(isinf(rad)))
    R=R/cos(alpha);
end
w=[R*exp(1i*alpha*[-1;1]);0]; 

res=zeros(kmax,1);
dofs=zeros(kmax,1);
for k=1:kmax
    n=4*k;
    N=2*(n+1);
    x=chebpts(N,1);
    %scl=atanh(1-exp(-4*(sqrt(n/8)-1)));
    %x=tanh(scl*linspace(-1,1,N));
    [zs,un]=polypts(w,rad,x);

    % Boundary data
    tol=2*eps; %tol=inf;
    rs=abs(zs);
    ts=angle(zs);
    b1=abs(rs/abs(w(1))-1)<=tol;
    
    mass=ones(size(zs));
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    
    bctype(b1)=1;
    bcdata(b1)=uex(rs(b1),ts(b1));
    
    [ufun,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,[],[],mass);
    res(k)=norm(r);
    dofs(k)=n;
end


r=linspace(0,R,np);
th=linspace(-alpha,alpha,np);
[rr,tt]=ndgrid(r,th);
xx=rr.*cos(tt);
yy=rr.*sin(tt);
zz=rr.*exp(1i*tt);
pout=(xx>1 & all(isinf(rad)));
[psih,uh]=ufun(zz);

bnd=R*exp(1i*th);
if(all(isinf(rad)))
    bnd=1+1i*sin([-1,1]*alpha);
end

% Plot stream function
cs=fex(rmax,0);
figure(1); clf;

psi=real(fex(rr,tt)); psi(pout)=nan;
subplot(1,2,1); contour(xx,yy,psi,cs,'k'); axis equal; xlim([0,R]);
hold on; plot(bnd,'--k'); hold off; grid off; axis off; axis tight;
title('Exact');

cs=real(ufun(rmax));
psi=real(psih); psi(1,:)=0; psi(end,[1,end])=0; psi(pout)=nan;
subplot(1,2,2); contour(xx,yy,psi,cs,'k'); axis equal; xlim([0,R]);
hold on; plot(bnd,'--k'); hold off; grid off; axis off; axis tight;
title('Numerical');

if(ifprint), print('-depsc',sprintf('%s_%s',prefix,'psi')); end

figure(2); clf
semilogy(dofs,res,'.-k',lw,2,ms,30);
ylim([1E-15,1]);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('Polynomial order, $N$'); title('Residual');
ylim([1E-20,1E0]);
if(ifprint)
    print('-depsc',sprintf('%s_%s',prefix,'res'));
end

figure(3); clf;
[~,duh]=ufun(rmax);
loglog(rmax,abs(duh),'.k',rmax,abs(uex(rmax,0)),'-k',lw,2,ms,20);
title('Velocity');
xlabel('Distance from corner, $x$'); legend('Numerical','Exact','Location','Best');
grid on; set(gca,'xminorgrid','off','yminorgrid','off'); 
if(ifprint)
    print('-depsc',sprintf('%s_%s',prefix,'vel'));
end


figure(4);
eddy_hunter(ufun,w([2,3,1]),6);
end