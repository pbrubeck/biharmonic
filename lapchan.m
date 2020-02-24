function [res]=lapchan(kmax,iprob)
ifprint=false;
ifslow=false;

nc=10;
nplt=16;
if(nargin<2),iprob=1;end

L1=2; H1=1;
r=14; A=2/(r-1); z0=1i*(H1+A);

% Corners of finite domain
w=[0; -1i*H1; L1-1i*H1; L1+1i*H1; -L1+1i*H1; -L1];
% Corners of unbounded sides
uw=w(abs(real(w))==L1);
ut=sign(real(uw));
uw=uw-real(uw)./(1+(imag(uw)<=0)); 

s1=1/L1;
ftop=@(z) 0*z;utop=@(z) 0*z;
if(iprob==1)
    ftop=@(z) -1i*(((z+1i)/2).*(1+tanh(s1*z))/2  + z.*(1-tanh(s1*z))/2 );
    utop=@(z)   2*(((z+1i)/2).*(s1*sech(s1*z).^2)/2 + z.*(-s1*sech(s1*z).^2)/2 + (1+tanh(s1*z))/4+ (1-tanh(s1*z))/2);
end
%ftop=@(z) 0*z;utop=@(z) 0*z;

kmin=1; 
res=zeros(kmax,1);
dofs=zeros(kmax,1);
npol=repmat(kmin,4,1);

tic;
%kmin=kmax; npol(:)=kmin; npol=[11,5,5,5];
for k=kmin:kmax
    if(ifslow), npol(:)=k; end
    n=0;

    disp(npol');
    nub=max(npol(3:4))^2;
    sigma=4; beta=(1-exp(-sigma*sqrt(nub)))/(A+2); h=1/3;
    
    kk=h/2:h:nub;
    t=beta*exp(-sigma*(sqrt(nub)-sqrt(1:nub)));
    tt=(1/beta)*(exp(sigma*(sqrt(nub)-sqrt(kk)))-1);
%     tt=tt(tt<1E6);
    pol2=1i*t;
    pol2=[pol2;-pol2];
    pol2=z0+1./pol2;
    z2=repmat(uw,1,length(tt))+ut*tt;
    hol=-2i/(A+2);
    hol=z0+1./hol;
    
    un2=repmat(1i*(2*(imag(uw)>0)-1).*ut,1,length(tt));
    id2=repmat(3,length(uw),length(tt));
    wq2=ones(size(z2));
   

    [pol1,z1,un1,id1,wq1,nps] = adapt_poles(npol.^2,w([1:3,end]));
    bin=(id1<=2);
    zs=[z1(bin);z2(:)];
    un=[un1(bin);un2(:)];
    id=[id1(bin);id2(:)];
    mass=[wq1(bin);wq2(:)];
    if(iprob==1)
        pol=[pol1(1:sum(nps(1:2))); pol2(:)];
    else
        pol=pol2;
    end

    %plot(zs,'.k'); hold on; plot(pol,'.r'); hold off; axis equal; xlim([-5,5]); gg
    
    xs=real(zs);
    ys=imag(zs);
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    tol=1E-15;
    top=(abs(ys-H1)<=tol);
    bot=~top;

    bcdata(top)=exp(-xs(top).^2);
    if(iprob==1)
    bcdata(top)=1-real(ftop(zs(top)));
    bcdata(bot)=0-real(ftop(zs(bot)));
    end

    mass(:)=1;
    [ufun,dofs(k),r]=lap(n,zs,un,bctype,bcdata,w,pol,hol,mass);
    res(k)=norm(r);
    ri=full(sparse(id,ones(size(r)),r.^2,max(id),1));
    kk=adapt_hist(ri);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
tsol=toc;
wx=real(w); w=w+(10/max(wx)-1)*wx;
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

psi=(1+1i)*nan(size(zz)); uu=psi;
tic;
[psi(ib),uu(ib)]=ufun(zz(ib));
tval=toc*1E3/nnz(ib);
psi(ib)=psi(ib)+ftop(zz(ib));
uu(ib)=conj(2i*uu(ib)+utop(zz(ib)));

lw='Linewidth'; ms='markersize'; fc='facecolor'; fs='fontsize';
ctol=0;
cs=[linspace(min(real(psi(:))),ctol,nc),0,linspace(ctol,max(real(psi(:))),nc)];


figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
pcolor(real(zz),imag(zz),real(psi)); hold on;
%contour(real(zz),imag(zz),real(psi),nc,'k',lw,1); hold on;
%contour(real(zz),imag(zz),imag(psi),nc,'k',lw,1); hold on;
%au=abs(uu); quiver(real(zz),imag(zz),real(uu)./au,imag(uu)./au,'k');

hold off; grid off; axis off;
colormap(jet(256)); shading interp; caxis([0,1]); %alpha(0.8); 
xlim([x1,x2]); ylim([y1,y2]); axis equal; 
%cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint),print('-depsc','lapchan_soln'); end

figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
surf(real(zz),imag(zz),real(psi)); hold on;
hold off; grid off; shading interp;
xlim([x1,x2]); ylim([y1,y2]);
%daspect([1,1,1]); %zlim([-2/3,2/3]);

colormap(jet(256)); %shading interp;
%cb=colorbar(); cb.TickLabelInterpreter='latex';

%if(ifprint), print('-depsc','lapchan_pres'); end

figure(3); clf;
subplot(1,2,2); 
plot(w([1:end,1]),'-k'); hold on;
plot(real(hol),imag(hol),'ob'); 
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10); 
hold off; axis equal; xlim([x1,x2]); 

subplot(1,2,1); semilogy(sqrt(dofs),res,'.-k',lw,2,ms,30);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-15,1E0]);
text(1,1E-11,sprintf('Solve time %.2f sec',tsol),fs,20);
text(1,1E-13,sprintf('Eval time %.2f ms/gridpoint',tval),fs,20);
if(ifprint), print('-depsc','lapchan_conv'); end

figure(4);
for k=1:max(id)
    plot(r(id==k)); hold on;
end
hold off; legend(num2str((1:max(id))'));
end