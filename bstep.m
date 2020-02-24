function [res]=bstep(kmax)
ifprint=false;
ifslow=false;

nc=10;
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

    mass(:)=1;
    [ufun,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,pol,[],mass);
    res(k)=norm(r);
    rtot=res(k)^2;

    ri=full(sparse([id;id],ones(size(r)),r.^2,numel(w),1));
    kk=adapt_hist(ri);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
tsol=toc;


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
tic;
[psi(ib),uu(ib),pp(ib)]=ufun(zz(ib));
tval=toc*1E3/nnz(ib);
pp(zz==0)=-5;


lw='Linewidth'; ms='markersize'; fc='facecolor'; fs='fontsize';
cs=[linspace(0.9*min(real(psi(:))),0,4),linspace(0,max(real(psi(:))),nc)];
figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
pcolor(real(zz),imag(zz),abs(uu)); hold on;
contour(real(zz),imag(zz),real(psi),cs(cs>=0),'k',lw,2); hold on;
contour(real(zz),imag(zz),real(psi),cs(cs<=0),'y',lw,2); hold on;
plot(w([1:end,1]),'-k',lw,2);
colormap(jet(256)); shading interp; axis off; caxis([0,1]); 

hold off; grid off; axis equal; 
cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint), print('-depsc','step_soln'); end


figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
surf(real(zz),imag(zz),real(pp)); hold on;
hold off; grid off; axis off;
xlim([-L1,L2]); ylim([-H1,H2]); 

zref=[-L1+1i*H1/2, L2];
[psiref,uref,pref]=ufun(zref);
daspect([1,1,abs(diff(real(pref)))]); caxis(sort(real(pref))); zlim([-5,16])
colormap(jet(256)); %shading interp;
cb=colorbar(); cb.TickLabelInterpreter='latex';
title('Pressure');
%if(ifprint), print('-depsc','step_pres'); end

figure(3); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
subplot(1,2,2); plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10);
hold off; axis equal; axis square;

subplot(1,2,1); 
semilogy(sqrt(dofs),res,'.-k',lw,2,ms,30);
axis tight; grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-15,1E0]);
xlim([0,10*ceil(0.1*sqrt(dofs(end)))]); ylim([1E-15,1E0]);
text(1,1E-11,sprintf('Solve time %.2f sec',tsol),fs,20);
text(1,1E-13,sprintf('Eval time %.2f ms/gridpoint',tval),fs,20);
if(ifprint), print('-depsc','step_conv'); end


L=0.5; wedge=[1i*L; 0; L]-1i;
figure(4); clf; eddy_hunter(ufun,wedge,2,256);

return
figure(4);
for k=1:numel(w)
    plot(r(id==k)); hold on;
end
hold off; legend(num2str((1:numel(w))'));
end