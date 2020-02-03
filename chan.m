function [res]=chan(nstep,iprob)
ifprint=true;
ifslow=false;

nc=10;
nplt=64;
if(nargin<2),iprob=1;end

L1=5; H1=1;

% Corners of finite domain
w=[L1+1i*H1; 1i*H1; -L1+1i*H1; -L1-1i*H1; -1i*H1; L1-1i*H1]; 
if(iprob==1)
    w=[L1+1i*H1; 1i*H1; -L1+1i*H1; -L1; 0; -1i*H1; L1-1i*H1];
end


% Corners of unbounded sides
uw=w(abs(real(w))==L1);
uw=uw-real(uw)/2;
% Tangents
ut=sign(real(uw));


ftop=@(z) 0*z;utop=@(z) 0*z;
if(iprob==1)
    ftop=@(z) -1i*(((z+1i)/2).*(1+tanh(z))/2  + z.*(1-tanh(z))/2 );
    utop=@(z)   1*(((z+1i)/2).*(sech(z).^2)/2 + z.*(-sech(z).^2)/2 + (1+tanh(z))/4+ (1-tanh(z))/2);
end
ftop=@(z) 0*z;utop=@(z) 0*z;

kmin=1; 
res=zeros(nstep,1);
dofs=zeros(nstep,1);
npol=repmat(kmin,size(w));

tic;
%kmin=kmax; npol(:)=kmin; npol(w==0)=10;
for k=kmin:nstep
    if(ifslow), npol(:)=k; end 
    n=4*k;

    disp(npol');
    [pol1,z1,un1,id1,wq1]=adapt_poles(npol.^2,w);
    sigma=log(4);
    beta=H1/1.5;
    
    
    nub=k*k;
    h=1/3;
    kk=h:h:nub;
    tt=beta*exp(sigma*(sqrt(nub)-sqrt(kk)));
    z2=repmat(uw,1,length(tt))+ut*(tt-beta);
    un2=repmat(1i*(2*(imag(uw)>0)-1).*ut,1,length(tt));
    id2=repmat(numel(w)-numel(uw)+(1:numel(uw))',1,length(tt));
    wq2=ones(size(z2));
    
    pol2=1i*(H1+beta*exp(sigma*(k-sqrt(1:k^2))));
    pol2=[pol2;-pol2];

    bin=abs(real(z1))<=L1/2;
    zs=[z1(bin);z2(:)];
    un=[un1(bin);un2(:)];
    id=[id1(bin);id2(:)];
    mass=[wq1(bin);wq2(:)];
    mass(:)=1;
    
    if(iprob==1)
        pol=[pol1(abs(real(pol1))<L1&imag(pol1)<0); pol2(:)];
    else
        pol=pol2;
    end
    
    xs=real(zs);
    ys=imag(zs);
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    
    tol=1E-15;
    top=(abs(ys-H1)<=tol);
    bot=~top;

    if(iprob==1)
        bcdata(top)=-2i/3;
    else
        bcdata(top)=-4i/3;
    end
    bcdata(bot)=0;
    
    [ufun,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,pol,[],mass);
    res(k)=norm(r);
    ri=full(sparse([id;id],ones(size(r)),r.^2,max(id),1));
    kk=adapt_hist(ri);    
    
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
tsol=toc;


w=w+real(w);
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

psi=(1+1i)*nan(size(zz)); uu=psi; pp=psi;
tic;
[psi(ib),uu(ib),pp(ib)]=ufun(zz(ib));
tval=toc*1E3/nnz(ib);
psi=psi+ftop(zz);
pp(zz==0)=nan;
uu=uu+conj(utop(zz));


lw='Linewidth'; ms='markersize'; fc='facecolor'; fs='fontsize';
cs=[linspace(min(real(psi(:))),0,4),linspace(0,max(real(psi(:))),nc)];
figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
pcolor(real(zz),imag(zz),abs(uu)); hold on;
contour(real(zz),imag(zz),real(psi),cs(cs>=0),'k',lw,2);
contour(real(zz),imag(zz),real(psi),cs(cs<=0),'y',lw,2);
plot(w([1:end,1]),'-k',lw,2);
%contour(real(zz),imag(zz),imag(psi),numel(cs),'k',lw,1);
% au=abs(uu);
% quiver(real(zz),imag(zz),real(uu)./au,imag(uu)./au,'k');

hold off; grid off; axis off;
colormap(jet(256)); shading interp; caxis([0,1]);
xlim([-L1,L1]); ylim([-H1,H1]); axis equal; 

%cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint), print('-depsc','chan_soln'); end


figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
surf(real(zz),imag(zz),real(pp)); hold on;
hold off; grid off; axis off;
xlim([-L1,L1]); ylim([-H1,H1]); 

colormap(jet(256)); %shading interp;
cb=colorbar(); cb.TickLabelInterpreter='latex';
title('Pressure');
%if(ifprint), print('-depsc','chan_pres'); end

figure(3); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
subplot(1,2,2); plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10);
hold off; axis equal; xlim([x1,x2]); 

subplot(1,2,1); 
%if(ifprint), subplot(1,1,1); end
semilogy(sqrt(dofs),res,'.-k',lw,2,ms,30);
axis tight; grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-15,1E0]);
xlim([0,10*ceil(0.1*sqrt(dofs(end)))]); ylim([1E-15,1E0]);
text(1,1E-11,sprintf('Solve time %.2f sec',tsol),fs,20);
text(1,1E-13,sprintf('Eval time %.2f ms/gridpoint',tval),fs,20);
if(ifprint), print('-depsc','chan_conv'); end

return
figure(4);
for k=1:numel(w)
    plot(r(id==k)); hold on;
end
hold off; legend(num2str((1:numel(w))'));
end