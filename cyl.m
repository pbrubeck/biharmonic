function []=cyl(kmax)
ifprint=false;
ifslow=true;

nc=10;
nplt=3*32;

L1=2;
L2=1/3; 
%L2=exp(0.25i*pi)*sqrt(2)/6;
%L2=(1/3)/imag(exp(1i*pi/3));


%w0=[L1+1i; 1i; -L1+1i]; w0=[L1; w0; -L1; conj(w0(end:-1:1))];
w0=[L1+1i; -L1+1i]; w0=[w0; conj(w0(end:-1:1))];
%w0=[L1+1i; 1i; -L1+1i]; w0=[w0; conj(w0(end:-1:1))];


w1=L2*[1+1i;-1+1i;-1-1i;1-1i];
%w1=L2*exp(1i*pi*(1:6)/3);
L2=abs(L2);


a0=exp(1i*0.5); a0=1;
x1=min(real(w0)); x2=max(real(w0));
y1=min(imag(w0)); y2=max(imag(w0));
w0=a0*w0; w1=a0*w1;


i0=1:length(w0);
i1=i0(end)+(1:length(w1));
w=[w0(:);w1(:)];

kmin=1; 
res=zeros(kmax,1);
dofs=zeros(kmax,1);
npol=repmat(kmin,size(w));

tic;
%kmin=kmax; npol(:)=kmin; npol(abs(w)<1)=10;

uex=@(z) conj(1-L2.^2./z.^2);
psiex=@(z) imag(z+L2.^2./z);

hol=mean(w(i1));
%w(i1)=[]; i1=[];
for k=kmin:kmax
    if(ifslow), npol(:)=k; end
    n=numel(w)*k;
    
    %disp(npol');
    [pol1,z1,un1,id1] = adapt_poles(npol(i0).^2,w(i0));
    if(isempty(i1))
        nh=ceil(n);
        pol2=hol; z2=exp(2i*pi*(1:nh)'/nh)*L2; un2=z2./abs(z2); id2=ones(size(z2));  
    else
        [pol2,z2,un2,id2] = adapt_poles(npol(i1).^2,w(i1(end:-1:1)));
    end
    pol=[pol1;pol2]; zs=[z1;z2]; un=[un1;un2]; id=[id1;id2+numel(i0)];

    xs=real(zs/a0);
    ys=imag(zs/a0);
    bctype=ones(size(zs));
    bcdata=zeros(size(zs));
    
    tol=1E-15;
    bot=abs(ys-y1)<=tol;
    top=abs(ys-y2)<=tol;
    left=abs(xs-x1)<=tol;
    right=abs(xs-x2)<=tol;
    wall=top|bot|left|right;
    hole=~wall;
    
    bctype(hole)=0;
    bctype(wall)=1;

    bcdata(wall)=uex(zs(wall));  

    
    bctype(wall)=1; bcdata(wall)=1; bcdata(wall)=1-abs(ys(wall)).^2;
    %bctype(right)=2; bcdata(right)=0;
    %bctype(top|bot)=4; bcdata(top|bot)=0;
    
    [goursat,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,pol,hol);
    res(k)=norm(r);
    rtot=res(k)^2;

    ri=full(sparse([id;id],ones(size(r)),r.^2,max(id),1));
    kk=(ri>rtot*0.5);    
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
    ri
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

[inp,onp]=inpolygon(real(zz),imag(zz),real(w(i0)),imag(w(i0))); ib0=(inp|onp);
if(isempty(i1))
    ib1=abs(zz-hol)<L2;
else
    [inp,onp]=inpolygon(real(zz),imag(zz),real(w(i1)),imag(w(i1))); ib1=(inp&~onp);
end
ib=ib0&~ib1;
%ib=logical(conv2(ib,ones(3),'same')-ib);

psi=(1+1i)*nan(size(zz)); uu=psi; pp=psi;
[psi(ib),uu(ib),pp(ib)]=goursat(zz(ib));

uu=uu.*min(3./abs(uu),1);
pp=pp.*min(10./abs(pp),1);

%zz(~ib)=nan; psi=imag(zz+L2.^2./zz); uu=(1-L2.^2./zz.^2);
%uu=uu-uex(zz); psi=psi-psiex(zz);

lw='Linewidth'; ms='markersize'; ctol=1E-8;
cs=[linspace(min(real(psi(:))),ctol,nc),0,linspace(ctol,max(real(psi(:))),nc)];
cs=linspace(-2/3,2/3,2*nc+1);

figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
pcolor(real(zz),imag(zz),abs(uu)); hold on;
contour(real(zz),imag(zz),real(psi),numel(cs),'k',lw,1); hold on;
%contour(real(zz),imag(zz),imag(psi),numel(cs),'k',lw,1); hold on;
if(isempty(i1))
    plot(L2*exp(2i*pi*(0:64)/64),'-k',lw,1); 
else
    plot(w(i1([1:end,1])),'-k',lw,1); 
end
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
%daspect([1,1,4/3]); %zlim([-2/3,2/3]);

colormap(jet(256)); %shading interp;
cb=colorbar(); cb.TickLabelInterpreter='latex';
title('Stream function');

if(ifprint), print('-depsc','cyl_pres'); end

figure(3); clf;
subplot(1,2,1); semilogy(sqrt(dofs),res,'.-k',lw,1,ms,20);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-15,1E0]);

subplot(1,2,2); 
plot(w(i0([1:end,1])),'-k'); hold on;
if(~isempty(i1)), plot(w(i1([1:end,1])),'-k'); end
plot(real(hol),imag(hol),'ob');
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10); hold off; axis equal;
if(ifprint)
    print('-depsc','ldc_conv');
end

figure(4);
for k=1:max(id)
    plot(r(id==k)); hold on;
end
hold off; legend(num2str((1:max(id))'));
end