function [] = moffatt(alpha)
if(nargin<1)
    alpha=pi/18;
end
ifprint=false;
prefix='va';
ms='markersize';
kmax=12;
ne=10; % number of eddies
nc=5;  % number of contours per eddy

sa=sin(2*alpha)/(2*alpha);
g=@(x,y) [sin(x).*cosh(y)+sa*x; cos(x).*sinh(y)+sa*y];
x=[4.21; 2.26]; 
x=NewtonRaphson(g,x);
lambda=1+([1,1i]*x)/(2*alpha);

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
uex=@(r,th) U*real(((r/r0).^lambda).*f(th));
grad=@(r,th) (U/r0)*exp(-1i*th).*(real((r/r0).^(lambda-1).*f(th)*lambda) + ...
                              -1i*real((r/r0).^(lambda-1).*ff(th)));

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
    %[pol,zs,un]=newpoles(N,w,rad);
    pol=[];[zs,un]=polypts(w,rad,chebpts(N,1));

    % Boundary data
    tol=1E-4; %tol=inf;
    rs=abs(zs);
    ts=angle(zs);
    b1=abs(rs-abs(w(1)))<tol; 
    u0=zeros(size(zs));
    u1=zeros(size(zs));
    u0(b1)=uex(rs(b1),ts(b1));
    u1(b1)=grad(rs(b1),ts(b1));
    
    [u,res(k),dofs(k)]=bih(n,zs,un,u0,u1,pol);
    dofs(k)=n;
end

np=256;
r=linspace(0,R,np);
th=linspace(-alpha,alpha,np);
[rr,tt]=ndgrid(r,th);
xx=rr.*cos(tt);
yy=rr.*sin(tt);
zz=rr.*exp(1i*tt);
pout=(xx>1 & all(isinf(rad)));
[uh,duh]=u(zz);

bnd=R*exp(1i*th);
if(all(isinf(rad)))
    bnd=1+1i*sin([-1,1]*alpha);
end

% Plot stream function
cs=uex(rmax,0);
figure(1); clf;

uu=real(uex(rr,tt)); uu(pout)=nan;
subplot(2,1,1); contour(xx,yy,uu,cs,'k'); axis equal; xlim([0,R]);
hold on; plot(bnd,'--k'); hold off; 

uu=real(uh); uu(1,:)=0; uu(end,[1,end])=0; uu(pout)=nan;
subplot(2,1,2); contour(xx,yy,uu,cs,'k'); axis equal; xlim([0,R]);
hold on; plot(bnd,'--k'); hold off;
if(ifprint)
    print('-depsc',sprintf('%s_%s',prefix,'psi'));
end

figure(2); clf
semilogy(dofs,res,'.-k',ms,20);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('Polynomial order $n$');
ylim([1E-15,1E0]);
if(ifprint)
    print('-depsc',sprintf('%s_%s',prefix,'res'));
end

figure(3); clf;
[~,duh]=u(rmax);
loglog(rmax,abs(duh),'.k',rmax,abs(grad(rmax,0)),'-k',ms,20);
xlabel('$x$'); legend('Numerical','Exact','Location','Best');
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
if(ifprint)
    print('-depsc',sprintf('%s_%s',prefix,'vel'));
end
end