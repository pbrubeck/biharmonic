function [res]=ldc(kmax)
ifprint=true;
ifslow=false;

w=[1+1i;-1+1i;-1-1i;1-1i];
%w=[1;1+1i;1i;-1+1i;-1;-1-1i;-1i;1-1i];

c0=0;
%A=5; v1=2; v2= 1; cs=[c0,-0.36,-0.066,0.001,0.0003,0.18,0.038,-0.0005,-0.00016];
%A=5; v1=1; v2= 0; cs=[c0,-0.18,-0.028,0.0005,0.0001,-1.4E-6 ,-2E-7 ,4E-9 ,1E-9];
%A=5; v1=1; v2=-1; cs=[c0,-0.18,-0.038,0.00052,0.0002,3.73E-5];

%A=1; v1=2; v2= 1; cs=[c0,-0.36,-0.194,0.088,0.15];
A=1; v1=1; v2= 0; cs=[c0,-0.19,-0.1,-0.016];
%A=1; v1=1; v2=-1; cs=[c0,-0.23579,-0.162,-0.0124];

w=real(w)+1i*A*imag(w);

res=zeros(kmax,1);
dofs=zeros(kmax,1);
kmin=1;
npol=ones(size(w));

%kmin=kmax; npol(:)=kmin; npol(imag(w)>0)=10;
tic;
for k=kmin:kmax
    n=k*numel(w);
    if(k>1), n=n+2*k; end

    if(ifslow), npol(:)=k; end
    %disp(npol.')
    [pol,zs,un,id,mass]=adapt_poles(npol.^2,w);
    
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    tol=1E-15;
    top=abs(un-1i)<=tol;
    bot=abs(un+1i)<=tol;
    bctype(top|bot)=1;
    bcdata(top)=v1;
    bcdata(bot)=v2;
    
    %mass(:)=1;
    [ufun,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,pol,[],mass);
    res(k)=norm(r);
    
    ri=full(sparse([id;id],ones(size(r)),r.^2,numel(w),1));
    kk=adapt_hist(ri);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
tsol=toc;

% Plotting
nplt=128;
xx=linspace(min(real(w)),max(real(w)),nplt);
yy=linspace(min(imag(w)),max(imag(w)),nplt);
[xx,yy]=ndgrid(xx,yy); zz=xx+1i*yy;

tic;
[psi,uu]=ufun(zz);
tval=toc*1E3/numel(zz);
uu([1,end],[1,end])=0;

nc=3;
zc=(-1-A*1i)+0.05*(1+A*1i)*(1:nc)/nc;
cs=[cs,real(ufun(zc))];
%cs=[-cs,linspace(min(max(0,psi(:))),max(max(0,psi(:))),5)];


lw='Linewidth'; ms='markersize'; fs='fontsize';

figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
pcolor(real(zz),imag(zz),abs(uu)); hold on;  
colormap(jet(256)); shading interp; caxis([0,1]); %alpha(0.8); 
contour(real(zz),imag(zz),real(psi),cs(cs>=0),'y',lw,1.5);
contour(real(zz),imag(zz),real(psi),cs(cs<=0),'k',lw,1.5); 
plot(w([1:end,1]),'k',lw,1.6);
hold off; grid off; axis equal; axis off; axis tight;
%cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint), print('-depsc','ldc_soln'); end

figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
subplot(1,2,2); plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10); 
hold off; axis equal; axis square;

subplot(1,2,1);
subplot(1,1,1); 
semilogy(sqrt(dofs),res,'.-k',lw,2,ms,30);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); 
xlim([0,10*ceil(0.1*sqrt(dofs(end)))]); ylim([1E-15,1E0]);
text(2,1E-11,sprintf('Solve time %.2f sec',tsol),fs,20);
text(2,1E-13,sprintf('Eval time %.2f ms/gridpoint',tval),fs,20);
if(ifprint), print('-depsc','ldc_conv'); end
end