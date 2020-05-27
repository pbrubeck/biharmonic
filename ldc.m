function [res]=ldc(kmax,ifprint)
if(nargin<2), ifprint=false; end
ifslow=false;
myratio = 1+(1+sqrt(5))/2;
if(ifprint),close all; end
function []=myprint(filename,rat)
    if(nargin>1)
        pos = get(gcf,'Position');
        if(rat>1)
            pos(3) = rat*pos(4);
        else
            pos(4) = (1/rat)*pos(3);
        end
        set(gcf,'Position',pos);
    end
    if(ifprint)
        drawnow;
        set(gcf,'Renderer','painters'); 
        print('-depsc',filename); 
        set(gcf,'Renderer','opengl'); 
    end
end


w=[1+1i;-1+1i;-1-1i;1-1i];
%w=[1;1+1i;1i;-1+1i;-1;-1-1i;-1i;1-1i];

c0=-5E-11;
%A=5; v1=2; v2= 1; cs=[-0.36,-0.066,0.001,0.0003,0.18,0.038,-0.0005,-0.00016];
%A=5; v1=1; v2= 0; cs=[-0.18,-0.028,0.0005,0.0001,-1.4E-6 ,-2E-7 ,4E-9 ,1E-9];
%A=5; v1=1; v2=-1; cs=[-0.18,-0.038,0.00052,0.0002,3.73E-5];

%A=1; v1=2; v2= 1; cs=[-0.36,-0.194,0.088,0.15];
A=1; v1=1; v2= 0; cs=[c0,-0.19,-0.1,-0.016];
%A=1; v1=1; v2=-1; cs=[-0.23579,-0.162,-0.0124];

w=real(w)+1i*A*imag(w);

res=zeros(kmax,1);
dofs=zeros(kmax,1);
kmin=1;
npol=ones(size(w));

%kmin=kmax; npol(:)=kmin; npol(imag(w)>0)=10;
tic;
for k=kmin:kmax
    n=k*numel(w);

    %if(k>1), n=n+2*k; end
        
    if(ifslow), npol(:)=k; end
    %disp(npol.')
    
    x=chebpts(k+1,1);
    rad = inf(size(w));
    [zs1,un1,id1,mass1] = polypts(w,rad,x);
    mass1(:)=1;
    
    mpol = npol.^2;
    [pol,zs,un,id,mass]=adapt_poles(mpol,w);
    
    %pol = pol(imag(pol)<0);
    kd = id>0;
    zs=[zs(kd);zs1];
    un=[un(kd);un1];
    id=[id(kd);id1];
    mass=[mass(kd);mass1];
    
    
    bctype=zeros(size(zs));
    bcdata=zeros(size(zs));
    tol=1E-15;
    top=abs(un-1i)<=tol;
    bot=abs(un+1i)<=tol;
    
    bctype(top|bot)=1;
    bcdata(top)=v1-0*real(zs(top)).^4;
    bcdata(bot)=v2;
    
    
    mass(id>2)=1;
    [ufun,dofs(k),r,pol]=bihstokes(n,zs,un,bctype,bcdata,w,pol,[],mass);
    res(k)=norm(r);

    ri=full(sparse([id;id],ones(size(r)),r.^2,numel(w),1));
    kk=adapt_hist(ri);
    npol(kk)=npol(kk)+ceil(sqrt(npol(kk)));
    npol=npol+1;
end
tsol=toc;
lw='Linewidth'; ms='markersize'; fc='facecolor'; fs='fontsize';

figure(3);
eddy_asymptotic(w([2,3,4,1]),ufun);
myprint('ldc_loglog', myratio);
%return

figure(2); clf; if(ifprint), set(gcf,'Renderer', 'Painters'); end
L=0.2*A; wedge=[1i*L; 0; L]+w(3); ne=0; nplt=128; nc=4; stol=res(end)/10;
subplot(1,2,1); cs1=eddy_hunter(ufun,wedge,ne,nplt,nc,stol); title('Second Eddy');
%subplot(1,2,2); eddy_hunter(ufun,-conj(wedge),ne,nplt,nc,stol);
L=0.0125*A; wedge=[1i*L; 0; L]+w(3);
subplot(1,2,2); eddy_hunter(ufun,wedge,ne,nplt,nc,stol); title('Third Eddy');
myprint('ldc_eddy', myratio);


% Plotting
nplt=256;
xx=linspace(min(real(w)),max(real(w)),nplt);
yy=linspace(min(imag(w)),max(imag(w)),nplt);
[xx,yy]=ndgrid(xx,yy); zz=xx+1i*yy;

tic;
[psi,uu,pp]=ufun(zz);
tval=toc*1E3/numel(zz);
uu([1,end],[1,end])=0;

%cs1=[];

cs=[cs(:);reshape(cs1(abs(cs1)>stol),[],1)];
tc=1E-3;


figure(1); clf;
subplot(1,2,1);
pcolor(real(zz),imag(zz),abs(uu)); hold on;
contour(real(zz),imag(zz),real(psi),cs(abs(cs)>=tc),'k',lw,1.5); hold on;
contour(real(zz),imag(zz),real(psi),cs(abs(cs)<=tc),'y',lw,1.5); hold on;
%plot(w([1:end,1]),'-k',lw,2);
colormap(jet(256)); shading interp; caxis([0,1]); 
%colorbar;

plot(real(pol),imag(pol),'.r',ms,10); 
hold off; grid off; axis equal; axis off;
LR = 1.5;
xlim([-LR,LR]); ylim([-LR,LR]);
%cf=colorbar(); cf.TickLabelInterpreter='latex';
set(gca,'Position',get(gca,'OuterPosition'));


%return
subplot(1,2,2); 
semilogy(sqrt(dofs),res,'.-k',lw,2,ms,30);
axis square; grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('$\sqrt{4N}$'); title('Weighted residual'); ylim([1E-15,1E0]);
xlim([0,10*ceil(0.1*sqrt(dofs(end)))]); ylim([1E-15,1E0]);
text(1,1E-09,sprintf('dim($A$) = %d$\\times$%d',numel(r),dofs(end)),fs,14);
text(1,1E-11,sprintf('Solve time %.2f sec',tsol),fs,14);
text(1,1E-13,sprintf('Eval time %.2f ms/gridpoint',tval),fs,14);
myprint('ldc', myratio);

return
subplot(1,2,1); plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k',lw,1,ms,10); plot(pol,'.r',lw,1,ms,10);
hold off; axis equal; axis square;
end