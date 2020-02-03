function [res]=annulus(kmax)
ifprint=false;
ifslow=true;

nc=10;
nplt=16*8;

R1=1; R2=2;

kmin=1; 
res=zeros(kmax,1);
dofs=zeros(kmax,1);
npol=0;

tic;
%kmin=kmax; npol(:)=kmin; npol(abs(w)<1)=10;

hol=0;
w=[];
for k=kmin:kmax
    if(ifslow), npol(:)=k; end
    n=2*k;
    n=2;

    nn=3*k^2; 
    pol1=[]; pol2=[];
    z1=R1*exp(2i*pi*(1:nn)'/nn); z2=R2*exp(2i*pi*(1:nn)'/nn);
    un1=z1./abs(z1); un2=z2./abs(z2); id1=repmat(1,size(z1)); id2=repmat(2,size(z2)); 
    wq1=ones(size(z1)); wq2=ones(size(z2));
    
    pol=[pol1;pol2]; zs=[z1;z2]; un=[un1;un2]; 
    id=[id1;id2]; mass=[wq1;wq2];

    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    wall=abs(zs)>(R1+R2)/2;
    hole=~wall;
    
    bctype(hole)=0;
    bctype(wall)=0;
    
    bcdata(hole)=0-1i;
    bcdata(wall)=0;
    
    %bcdata=imag(un.*uex(zs))-1i*fex(zs);
    %bcdata(hole)=1i*mean(imag(bcdata(hole)));
    
    mass(:)=1;
    [ufun,dofs(k),r]=bihstokes(n,zs,un,bctype,bcdata,w,pol,hol,mass);
    %[ufun,dofs(k),r]=bih2(n,zs,un,bctype,bcdata,w,pol,hol,mass);
    res(k)=norm(r);
    rtot=res(k)^2;
        
    ri=sparse([id;id],ones(size(r)),r.^2,max(id),1);
    kk=(ri>rtot*0.5);
    %npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
toc

% Plotting
rr=linspace(R1,R2,nplt);
tt=linspace(0,2*pi,nplt);
[rr,tt]=ndgrid(rr,tt);
zz=rr.*exp(1i*tt);
[psi,uu,pp]=ufun(zz);

%psi=psi-fex(zz); uu=uu-conj(uex(zz));
%uu=uu.*min(1./abs(uu),1);

%zz(~ib)=nan; psi=imag(zz+L2.^2./zz); uu=(1-L2.^2./zz.^2);
%uu=uu-uex(zz); psi=psi-psiex(zz);

lw='Linewidth'; ms='markersize'; ctol=0E-8;
cs=[linspace(min(real(psi(:))),ctol,nc),0,linspace(ctol,max(real(psi(:))),nc)];
%cs=linspace(-2/3,2/3,2*nc+1);


figure(1); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
pcolor(real(zz),imag(zz),abs(uu)); hold on;
contour(real(zz),imag(zz),real(psi),numel(cs),'k',lw,1);
%contour(real(zz),imag(zz),imag(psi),numel(cs),'k',lw,1);
%au=abs(uu); quiver(real(zz),imag(zz),real(uu)./au,imag(uu)./au,'k');

hold off; grid off;
colormap(jet(256)); shading interp; alpha(0.8); caxis([0,max(abs(uu(:)))]);
axis equal; 

cb=colorbar(); cb.TickLabelInterpreter='latex';
if(ifprint),print('-depsc','cyl_soln'); end

figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
surf(real(zz),imag(zz),real(psi)); hold on;
hold off; grid off;
%daspect([1,1,4/3]); %zlim([-2/3,2/3]);

colormap(jet(256)); %shading interp;
cb=colorbar(); cb.TickLabelInterpreter='latex';
title('Stream function');

if(ifprint), print('-depsc','cyl_pres'); end

figure(3); clf;
semilogy(sqrt(dofs),res,'.-k',lw,1,ms,20);
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual'); ylim([1E-16,1E0]);
if(ifprint), print('-depsc','ldc_conv'); end
end