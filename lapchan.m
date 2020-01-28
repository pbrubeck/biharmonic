function [res]=lapchan(kmax,iprob)
ifprint=false;
ifslow=true;

nc=10;
nplt=16;
if(nargin<2),iprob=1;end

L1=2; H1=1;

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


s1=1/4;
ftop=@(z) 0*z;utop=@(z) 0*z;
if(iprob==1)
    ftop=@(z) -1i*(((z+1i)/2).*(1+tanh(s1*z))/2  + z.*(1-tanh(s1*z))/2 );
    utop=@(z)   2*(((z+1i)/2).*(s1*sech(s1*z).^2)/2 + z.*(-s1*sech(s1*z).^2)/2 + (1+tanh(s1*z))/4+ (1-tanh(s1*z))/2);
end
ftop=@(z) 0*z;utop=@(z) 0*z;

kmin=1; 
res=zeros(kmax,1);
dofs=zeros(kmax,1);
npol=repmat(kmin,size(w));

tic;
%kmin=kmax; npol(:)=kmin; npol(abs(w)<1)=10;
for k=kmin:kmax
    if(ifslow), npol(:)=k; end
    n=2*k; %n=0;

    %disp(npol');
    [pol1,z1,un1,id1,wq1]=adapt_poles(npol.^2,w);
    sigma=log(4);
    beta=L1;
    
    nub=k*k;
    h=1/3;
    kk=h:h:nub;
    tt=beta*exp(sigma*(sqrt(nub)-sqrt(kk)));
    z2=repmat(uw,1,length(tt))+ut*(tt-tt(end));
    un2=1i*(2*(imag(uw)>0)-1).*ut;
    id2=repmat(numel(w)-numel(uw)+(1:numel(uw))',1,length(tt));
    wq2=ones(size(z2));
    
    beta=beta;
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

    bcdata(top)=exp(-xs(top).^2);
    if(iprob==1)
    bcdata(top)=1-real(ftop(zs(top)));
    bcdata(bot)=0-real(ftop(zs(bot)));
    end
    
    [ufun,dofs(k),r]=lap(n,zs,un,bctype,bcdata,w,pol,[],mass);
    res(k)=norm(r);
    rtot=res(k)^2;

    ri=full(sparse(id,ones(size(r)),r.^2,max(id),1));
    kk=(ri>rtot*0.5);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
toc


w=w+99*real(w);
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
[psi(ib),uu(ib)]=ufun(zz(ib));
psi(ib)=psi(ib)+ftop(zz(ib));
uu(ib)=conj(2i*uu(ib)+utop(zz(ib)));

lw='Linewidth'; ms='markersize'; ctol=0;
cs=[linspace(min(real(psi(:))),ctol,nc),0,linspace(ctol,max(real(psi(:))),nc)];


figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
pcolor(real(zz),imag(zz),abs(uu)); hold on;
contour(real(zz),imag(zz),real(psi),nc,'k',lw,1); hold on;
contour(real(zz),imag(zz),imag(psi),nc,'k',lw,1); hold on;
%au=abs(uu); quiver(real(zz),imag(zz),real(uu)./au,imag(uu)./au,'k');

hold off; grid off;
colormap(jet(256)); shading interp; alpha(0.8); caxis([0,1]);
xlim([x1,x2]); ylim([y1,y2]); axis equal; 

cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint),print('-depsc','cyl_soln'); end

figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
surf(real(zz),imag(zz),real(psi)); hold on;
hold off; grid off;
xlim([x1,x2]); ylim([y1,y2]);
%daspect([1,1,1]); %zlim([-2/3,2/3]);

colormap(jet(256)); %shading interp;
cb=colorbar(); cb.TickLabelInterpreter='latex';
title('Stream function');

if(ifprint), print('-depsc','cyl_pres'); end

figure(3); clf;
subplot(1,2,1); semilogy(sqrt(dofs),res,'.-k',lw,1,ms,20);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-15,1E0]);

subplot(1,2,2); 
plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10); hold off; axis equal;
xlim([x1,x2]); ylim(3*[y1,y2]);
if(ifprint)
    print('-depsc','ldc_conv');
end

figure(4);
for k=1:max(id)
    plot(r(id==k)); hold on;
end
hold off; legend(num2str((1:max(id))'));
end