function [res]=bstep(kmax)
ifprint=false;
ifslow=false;

nplt=64;
L1=1; L2=5; 
H1=1; H2=1;
%w=[L2+1i*H2; -L1+1i*H2; -L1-1i*H1; L2-1i*H1];
%w=[L2+1i*H2; -L1+1i*H2; -L1; 0; -1i*H1; L2-1i*H1];
w=[L2+1i*H2; 1i*H2; -L1+1i*H2; -L1; 0; -1i*H1; L2-1i*H1];
%w=[L2+1i*H2; -L1+1i*H2; -L1; 0; -1i*H1; L2-1i*H1];

a0=exp(1i*0.5); a0=1;
w=a0*w;

kmin=1; 
res=zeros(kmax,1);
dofs=zeros(kmax,1);
npol=repmat(kmin,size(w));

tic;
%kmin=kmax; npol(:)=kmin; npol(w==0)=10+(kmin-4);
for k=kmin:kmax
    if(ifslow), npol(:)=k; end
    n=numel(w)*k;
    
    %disp(npol');
    [pol,zs,un,id,mass]=adapt_poles(npol.^2,w);
    

    xs=real(zs/a0);
    ys=imag(zs/a0);
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    
    tol=1E-15;
    bot=abs(ys+H1)<=tol;
    top=abs(ys-H2)<=tol;
    left=abs(xs+L1)<=tol;
    right=abs(xs-L2)<=tol;
    
    bctype(right)=2;
    bctype(left)=1;
    bctype(top)=1;
    bctype(bot)=0;

    bcdata(left) = (1-(2*ys(left)/H2-1).^2)*a0;
    %bcdata(right)= (1-((2*ys(right)-(H2-H1))/(H2+H1)).^2)*a0/2;
    if(length(w)==4)
        bcdata(left) = (1-ys(left).^2)*a0;
        bcdata(right)= 2*bcdata(right)*a0;
    end

    [goursat,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,pol,[],mass);
    res(k)=norm(r);
    rtot=res(k)^2;
    
    
    ri=full(sparse([id;id],ones(size(r)),r.^2,numel(w),1));
    kk=(ri>rtot*0.5);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
    
    %full(sparse(id,ones(size(mass)),mass,max(id),1))
end
toc

% Plotting
x1=min(real(w)); x2=max(real(w)); dx=x2-x1;
y1=min(imag(w)); y2=max(imag(w)); dy=y2-y1;
h=min(dx,dy)/nplt;
nx=ceil(1+dx/h); 
ny=ceil(1+dy/h);
x=linspace(x1,x2,nx);
y=linspace(y1,y2,ny);
[xx,yy]=ndgrid(x,y); 
zz=xx+1i*yy;
[inp,onp]=inpolygon(real(zz),imag(zz),real(w),imag(w)); ib=(inp|onp);
psi=nan(size(zz)); uu=nan(size(zz)); pp=nan(size(zz));
[psi(ib),uu(ib),pp(ib)]=goursat(zz(ib));


lw='Linewidth'; ms='markersize'; ctol=0;
cs=[linspace(min(real(psi(:))),ctol,4),linspace(ctol,max(real(psi(:))),20)];


figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end

pcolor(real(zz),imag(zz),abs(uu)); hold on;
plot(w([1:end,1]),'-k',lw,1);
contour(real(zz),imag(zz),real(psi),cs,'k',lw,1);
%contour(real(zz),imag(zz),imag(psi),numel(cs),'k',lw,1);
% au=abs(uu(:,:,j));
% quiver(real(zz(:,:,j)),imag(zz(:,:,j)),real(uu(:,:,j))./au,imag(uu(:,:,j))./au,'k');

hold off; grid off;
colormap(jet(256)); shading interp; alpha(0.8); caxis([0,1]);
xlim([-L1,L2]); ylim([-H1,H2]); axis equal; 

cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint)
    print('-depsc','ldc_soln');
    %export_fig('ldc_soln','-eps','-transparent');
end


figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
surf(real(zz),imag(zz),real(pp)); hold on;
hold off; grid off;
xlim([-L1,L2]); ylim([-H1,H2]); 

zref=[-L1+1i*H1/2, L2];
[psiref,uref,pref]=goursat(zref);
daspect([1,1,abs(diff(real(pref)))]); caxis(sort(real(pref))); zlim([-1/3,4/3]*real(pref(1)));
colormap(jet(256)); %shading interp;
cb=colorbar(); cb.TickLabelInterpreter='latex';
title('Pressure');

if(ifprint)
    print('-depsc','ldc_pres');
    %export_fig('ldc_soln','-eps','-transparent');
end

figure(3); clf;
subplot(1,2,1); semilogy(sqrt(dofs),res,'.-k',lw,1,ms,20);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-15,1E0]);

subplot(1,2,2); plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10); hold off; axis equal;
if(ifprint)
    print('-depsc','ldc_conv');
end

figure(4);
for k=1:numel(w)
    plot(r(id==k)); hold on;
end
hold off; legend(num2str((1:numel(w))'));

return
[yq,wq]=gauleg(-H1,H2,512);
zx=L2+1i*yq(:);
[sx,ux,px,fx]=goursat(a0*zx);

figure(5); plot(imag(zx),abs(ux+4*fx),'r');
disp(wq*px/sum(wq))
end