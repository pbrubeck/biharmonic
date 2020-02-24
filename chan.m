function [res]=chan(kmax,iprob)
ifprint=false;
ifslow=false;

nc=10;
nplt=32; Lplt=5;
if(nargin<2), iprob=1; end

L1=5; H1=1;
r=14; A=2/(r-1); z0=1i*(H1+A);

% Corners of finite domain
w=[0; -1i*H1; L1-1i*H1; L1+1i*H1; 1i*H1; -L1+1i*H1; -L1];
% Corners of unbounded sides
uw=w(abs(real(w))==L1);
ut=sign(real(uw));
uw=uw-real(uw)/2;
s2=2/3; s1=1;
ftop=@(z) 0*z;utop=@(z) 0*z;
if(iprob==1)
    % i = side (- or +)
    % j = fun (f or g)
    % k = diff (f or f')

    a=[0,0,2,-4/3; 1/3,1/2,0,-1/6]; % cubics in y
    deg=size(a,2)-1;
    coef=zeros(deg+1,2,2,2);
    for i=1:2
        coef(:,i,1,1)=[0;0.5i*a(i,3);0.75*a(i,4);0]; % f
        coef(:,i,2,1)=[1i*a(i,1);a(i,2);-0.5i*a(i,3);-0.25*a(i,4)]; % g
    end
    coef=reshape(coef,size(coef,1),4,2);
    coef(1:end-1,:,2)=spdiag(1:deg)*coef(2:end,:,1);    
    coef=reshape(coef,[],2,2,2);
    f=@(i,j,k,z) polyval(flipud(coef(:,i,j,k)),z); 
    pfun=@(i,z) (1+(2*i-3)*tanh(s1*z))/2;
    qfun=@(i,z) s1*((2*i-3)*(sech(s1*z)).^2)/2;
    ff = @(i,z) pfun(i,z).*(conj(z).*f(i,1,1,z) + f(i,2,1,z));
    ftop = @(z) ff(1,z)+ff(2,z);
    fu = @(i,z) qfun(i,z).*(conj(z).*f(i,1,1,z)+f(i,2,1,z))+...
                pfun(i,z).*(conj(z).*f(i,1,2,z)+f(i,2,2,z))+...
                -conj(pfun(i,z).*f(i,1,1,z));
    utop = @(z) fu(1,z)+fu(2,z);
    fp = @(i,z) pfun(i,z).*f(i,1,2,z)+qfun(i,z).*f(i,1,1,z);
    ptop = @(z) 4*conj(fp(1,z)+fp(2,z));
end
%ftop=@(z) 0*z;utop=@(z) 0*z;
function [psi,uu,pp]=vfun(z)
    [psi,uu,pp]=ufun(z);
    psi=psi-1i*ftop(z);
    uu=uu+conj(utop(z));
    pp=pp+ptop(z);
end

kmin=1; 
res=zeros(kmax,1);
dofs=zeros(kmax,1);
npol=repmat(kmin,4,1);

tic;
%kmin=kmax; npol(:)=kmin; npol=[10,kmax,kmax,kmax];
for k=kmin:kmax
if(ifslow), npol(:)=k; end
    n=0;
    disp(npol');
    nub=ceil((npol(4))^2);
    sigma=4; beta=(1-exp(-sigma*sqrt(nub)))/(A+2.001); h=1/3;

    kk=h/2:h:nub;
    t=beta*exp(-sigma*(sqrt(nub)-sqrt(1:nub)));
    tt=(L1/2)*(exp(sigma*(sqrt(nub)-sqrt(kk)))-1);
    
    cutoff=1/sqrt(eps);
    %cutoff=1/eps;
    tt=tt(tt<cutoff);
    %t=t(1./t<cutoff);

    pol2=[1i*t; -1i*t];
    pol2=z0+1./pol2;
    %pol2=pol2(abs(pol2)<10*cutoff);
    
    z2=repmat(uw,1,length(tt))+ut*tt;
    hol=-2i/(A+2);
    hol=z0+1./hol;
    %hol=[hol;-hol];
    %hol=[];
    
    un2=repmat(1i*(2*(imag(uw)>0)-1).*ut,1,length(tt));
    id2=repmat(4,length(uw),length(tt));
    wq2=ones(size(z2));
    %wq2=min(abs(z2).^2)./abs(z2).^2;

    mpol=ones(size(w));
    mpol([1,2,5])=npol(1:3).^2;
    [pol1,z1,un1,id1,wq1,nps] = adapt_poles(mpol,w);
    bin=(id1<=2|id1==5);
    id1(id1==5)=3;

    zs=[z1(bin);z2(:)];
    un=[un1(bin);un2(:)];
    id=[id1(bin);id2(:)]; 
    mass=[wq1(bin);wq2(:)];
    
    nps=cumsum(nps);
    pol=[pol1([1:nps(2),1+nps(4):nps(5)]); pol2(:)];
    
    %plot(zs,'.k'); hold on; plot(pol,'.r'); hold off; axis equal; xlim([-5,5]); gg
    
    xs=real(zs);
    ys=imag(zs);
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    top=(un==1i);
    bot=~top;

    if(iprob==1)
    bcdata(top)=-1i*(s2-imag(ftop(zs(top))))-imag(un(top).*utop(zs(top)));
    bcdata(bot)=-1i*(0-imag(ftop(zs(bot))))-imag(un(bot).*utop(zs(bot)));
    
    bctype(:)=1; bcdata=-conj(utop(zs));
    end

    %mass(id<=3)=1;
    mass(:)=1; 
    ifp0=true;
    %ifp0=false;
    [ufun,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,pol,hol,mass,ifp0);
    res(k)=norm(r);
    ri=full(sparse([id;id],ones(size(r)),r.^2,max(id),1));
    kk=adapt_hist(ri);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
tsol=toc;
wx=real(w); w=w+(Lplt/max(wx)-1)*wx;

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
[psi(ib),uu(ib),pp(ib)]=vfun(zz(ib));
tval=toc*1E3/nnz(ib);
psi=psi-real(vfun(0));
pp(zz==0)=nan;


%[psi0,u0,p0]=ufun(zp0);p0

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
xlim([x1,x2]); ylim([y1,y2]); axis equal; 

cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint), print('-depsc','chan_soln'); end


figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
surf(real(zz),imag(zz),real(psi)); hold on;
hold off; grid off; %axis equal; %axis off;
xlim([x1,x2]); ylim([y1,y2]); %caxis([0,s2]);

colormap(jet(256)); %shading interp;
cb=colorbar(); cb.TickLabelInterpreter='latex';
title('Stream function');
%if(ifprint), print('-depsc','chan_pres'); end

figure(3); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
subplot(1,2,2); 
plot(w([1:end,1]),'-k'); hold on;
plot(real(hol),imag(hol),'ob'); 
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


figure(5);
plot(real(zz(imag(zz)==1)),real(psi(imag(zz)==1))); hold on;
plot(real(zz(imag(zz)==0)),real(psi(imag(zz)==0))); hold on;
plot(real(zz(imag(zz)==-1)),real(psi(imag(zz)==-1))); hold off;


return;
L=0.5; wedge=[1i*L; 0; L]-1i;
figure(4); clf; eddy_hunter(@vfun,wedge,1,32);
end