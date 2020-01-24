function [res]=laphole(kmax)
ifprint=false;
ifslow=true;

nc=10;
nplt=15*16;

L1=1; H1=1;
L2=1/15; H2=L2;
hol=-1/3;

w1=[L1+1i*H1; 1i*H1; -L1+1i*H1]; w1=[L1; w1; -L1; conj(w1(end:-1:1))];
w1=[L1+1i*H1; -L1+1i*H1]; w1=[w1; conj(w1(end:-1:1))];

w2=[L2; L2+1i*H2; -L2+1i*H2; -L2; -L2-1i*H2; L2-1i*H2];
w2=[L2+1i*H2; -L2+1i*H2; -L2-1i*H2; L2-1i*H2];


w3=w2-conj(hol);
w2=w2+hol;

a0=exp(1i*0.5); a0=1;
x1=min(real(w1)); x2=max(real(w1));
y1=min(imag(w1)); y2=max(imag(w1));
w1=a0*w1; w2=a0*w2; w3=a0*w3;

i1=1:length(w1);
i2=i1(end)+(1:length(w2));
i3=i2(end)+(1:length(w2));
w=[w1(:); w2(:); w3(:)];

kmin=1; 
res=zeros(kmax,1);
dofs=zeros(kmax,1);
npol=repmat(kmin,size(w));

tic;
%kmin=kmax; npol(:)=kmin; npol(abs(w)<1)=10;
hol=[mean(w(i2)); mean(w(i3))];
for k=kmin:kmax
    if(ifslow), npol(:)=k; end
    n=numel(w)*k;

    %disp(npol');
    [pol1,z1,un1,id1,wq1] = adapt_poles(npol(i1).^2,w(i1));
    [pol2,z2,un2,id2,wq2] = adapt_poles(npol(i2).^2,w(i2(end:-1:1)));
    [pol3,z3,un3,id3,wq3] = adapt_poles(npol(i3).^2,w(i3(end:-1:1)));
    
    pol=[pol1;pol2;pol3]; zs=[z1;z2;z3]; un=[un1;un2;un3]; 
    id=[id1;id2+numel(i1);id3+numel(i1)+numel(i2)]; mass=[wq1;wq2;wq3];

    xs=real(zs/a0);
    ys=imag(zs/a0);
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    
    tol=1E-15;
    bot=abs(ys-y1)<=tol;
    top=abs(ys-y2)<=tol;
    left=abs(xs-x1)<=tol;
    right=abs(xs-x2)<=tol;
    wall=top|bot|left|right;
    hole=~wall;
    
    bctype(hole)=0;
    bctype(wall)=0;
    bcdata(hole)=sign(xs(hole));
    %bcdata(hole)=1;

    [ufun,dofs(k),r]=lap(n,zs,un,bctype,bcdata,w,pol,hol,mass);
    res(k)=norm(r);
    rtot=res(k)^2;

    ri=full(sparse(id,ones(size(r)),r.^2,max(id),1));
    kk=(ri>rtot*0.5);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
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

[inp,onp]=inpolygon(real(zz),imag(zz),real(w(i1)),imag(w(i1))); ib1=(inp|onp);
[inp,onp]=inpolygon(real(zz),imag(zz),real(w(i2)),imag(w(i2))); ib2=(inp&~onp);
[inp,onp]=inpolygon(real(zz),imag(zz),real(w(i3)),imag(w(i2))); ib3=(inp&~onp);

ib=ib1&~ib2&~ib3;
psi=(1+1i)*nan(size(zz)); uu=psi;
[psi(ib),uu(ib)]=ufun(zz(ib));
uu=uu.*min(10./abs(uu),1);

lw='Linewidth'; ms='markersize'; ctol=0;
cs=[linspace(min(real(psi(:))),ctol,nc),0,linspace(ctol,max(real(psi(:))),nc)];


figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
%pcolor(real(zz),imag(zz),abs(uu)); hold on;
contour(real(zz),imag(zz),real(psi),nc,'k',lw,1); hold on;
contour(real(zz),imag(zz),imag(psi),nc,'k',lw,1); hold on;

plot(w(i2([1:end,1])),'-k',lw,1); hold on;
plot(w(i3([1:end,1])),'-k',lw,1); hold on;


%au=abs(uu); quiver(real(zz),imag(zz),real(uu)./au,imag(uu)./au,'k');

hold off; grid off;
colormap(jet(256)); shading interp; alpha(0.8); caxis([0,max(abs(uu(ib)))]);
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
plot(w(i1([1:end,1])),'-k'); hold on;
if(~isempty(i2)), plot(w(i2([1:end,1])),'-k'); end
if(~isempty(i3)), plot(w(i3([1:end,1])),'-k'); end
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