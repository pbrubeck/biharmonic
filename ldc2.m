function ldc2(kmax)
ifprint=false;
npoly=4; 
w=exp((2i*pi/npoly)*(1:npoly)')/(exp(2i*pi/npoly)-1); 
w0=0;
w=w-w0; 
w=w./max(max(imag(w)),max(real(w)));

c0=1E-10;
%A=5; v1=2; v2= 1; cs=[c0,0.36,0.066,-0.001,-0.0003,-0.18,-0.038,0.0005,0.00016];
%A=5; v1=1; v2= 0; cs=[c0,0.18,0.028,-0.0005,-0.0001,1.4E-6 , 2E-7 ,-4E-9 ,-1E-9];
%A=5; v1=1; v2=-1; cs=[c0,0.18,0.038,-0.00052,-0.0002,-3.73E-5];

%A=1; v1=2; v2= 1; cs=[c0,0.36,0.194,-0.088,-0.15];
A=1; v1=1; v2= 0; cs=[c0,0.19,0.1,0.016];
%A=1; v1=1; v2=-1; cs=[c0,0.23579,0.162,0.0124];

w=real(w)+1i*A*imag(w);
res=zeros(kmax,1);
dofs=zeros(kmax,1);
for k=1:kmax
    N=k^2;
    n=ceil(N/2);
    [pol,zs,un]=newpoles(N,w);
    tol=1E-4;
    u0=zeros(size(zs));
    u1=zeros(size(zs));
    b1=abs(un-1i)<=tol;
    b2=abs(un+1i)<=tol;
    u1(b1)=v1*un(b1);
    u1(b2)=-v2*un(b2);
    [u,res(k),dofs(k)]=bih(n,zs,un,u0,u1,pol,w);
end

% Plotting
nplt=128;
xx=linspace(min(real(w)),max(real(w)),nplt);
yy=linspace(min(imag(w)),max(imag(w)),nplt);
[xx,yy]=ndgrid(xx,yy); zz=xx+1i*yy;
[uu,du]=u(zz);
du([1,end],[1,end])=0;

nc=3;
zc=(-1-A*1i)+0.05*(1+A*1i)*(1:nc)/nc;
cs=[cs,u(zc)];

lw='Linewidth';
ms='markersize';

figure(1); clf;
set(gcf,'Renderer', 'Painters');
pcolor(real(zz),imag(zz),abs(du)); hold on;  
colormap(jet(256));  shading interp; alpha(0.8);
contour(real(zz),imag(zz),real(uu),cs,'k',lw,1); hold off;
axis equal; grid off;
cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint)
    print('-depsc','ldc_soln');
    %export_fig('ldc_soln','-eps','-transparent');
end

figure(2); clf;
subplot(1,2,1); semilogy(sqrt(dofs),res,'.-k',lw,1,ms,20);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-15,1E0]);

subplot(1,2,2); plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10); hold off; axis equal;
if(ifprint)
    print('-depsc','ldc_conv');
end
end