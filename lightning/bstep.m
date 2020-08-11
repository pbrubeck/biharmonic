function [res]=bstep(ifprint)
if(nargin<1), ifprint=false; end
ifslow=false;
ifstats=true;
myratio = 1+(1+sqrt(5))/2;
%if(ifprint),close all; end
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

stol=1E-10;
nplt=128;
L1=1; L2=5; 
H1=1; H2=1;

top = @(z) nan(size(z))+2i/3;
bot = @(z) nan(size(z));
top = @(z) zeros(size(z));
bot = @(z) zeros(size(z));

inlet = @(z) 1-(2*imag(z)-1).^2;
outlet = @(z) (1-imag(z).^2)/2;

w=[L2+1i*H2; -L1+1i*H2; -L1; 0; -1i*H1; L2-1i*H1];
ubc = {top; inlet; bot; bot; bot; outlet};

%w=[L2+1i*H2; 1i*H2; -L1+1i*H2; -L1; 0; -1i*H1; L2-1i*H1];
%ubc = {top; top; inlet; bot; bot; bot; outlet};

figure(1);
[ufun, maxerr, f, Z, Zplot, A, pol] = stokes(w, ubc,'rel', 'tol', stol);
f0 = min(real(f(0,Z)));
ffun = @(z) real(f(0,z))-f0;

cs1=[0,0];
figure(2); clf;
[~,ic]=min(abs(w+1i*H2));
L=0.5; wedge=[1i*L; 0; L]+w(ic); ne=0; nce=4; 
subplot(1,2,1); cs1=eddy_hunter(ffun,wedge,ne,nplt,nce,stol); title('First Eddy');
L=1/36; wedge=[1i*L; 0; L]+w(ic);
subplot(1,2,2); eddy_hunter(ffun,wedge,ne,nplt,nce,stol); title('Second Eddy');
myprint('step_eddy',myratio); 

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
psi=nan(size(zz)); 
uu=nan(size(zz)); 
psi(ib)=ffun(zz(ib));
uu(ib)=ufun(zz(ib));

nc=10;
cs=[linspace(stol,max(real(psi(:))),nc)';cs1(abs(cs1)>stol)]; tc=1E-3;
lw='Linewidth'; ms='markersize'; fc='facecolor'; fs='fontsize';

figure(1); clf; 
pcolor(real(zz),imag(zz),abs(uu));  hold on;
contour(real(zz),imag(zz),real(psi),cs(abs(cs)>=tc),'k',lw,1.5); hold on;
contour(real(zz),imag(zz),real(psi),cs(abs(cs)<=tc),'y',lw,1.5); hold on;
plot([0,5],[-1,-1],'-y','linewidth',1.5)
%plot(w([1:end,1]),'-k',lw,2);
colormap(parula(256)); shading interp; axis off; caxis([0,1]); 

plot(real(pol),imag(pol),'.r',ms,10);
hold off; grid off; axis equal; 
%cf=colorbar(); cf.TickLabelInterpreter='latex';

xlim([-1.5, 5.5]);
ylim([-1.5, 1.5]);
myprint('step',myratio); 
end