function [cs] = eddy_hunter(ufun,w,ne,nplt,nc,stol)
% EDDY_HUNTER plot Moffatt eddies
if(nargin<3), ne=0; end
if(nargin<4), nplt=64; end
if(nargin<5), nc=4; end
if(nargin<6), stol=1E-11; end

va='VerticalAlignment';
ha='HorizontalAlignment';
lw='LineWidth';

dw=w([1,3])-w(2);
un=mean(dw./abs(dw));
%un=un/abs(un);
tha=angle(un);
R=mean(abs(dw)); 
a=R*un;

b = log(1/200);
map=@(x) w(2)+a*exp(b*(1-x));
map=@(x) w(2)+a*x;
f0=real(ufun(map(0))); 
f0=0;

cf=chebfun(@(x) real(ufun(map(x))-f0), [0,1],32); % low frequency
[fe,re]=minandmax(cf,'local');
[ct,kk]=sort(abs(fe));
kk=kk(ct>stol);
if(isempty(kk))
    disp('Could not find an eddy!');
    return
end
rmin=re(max(kk(1)-1,1));
%fe=fe(re>=rmin); re=re(re>=rmin);
rz=roots(cf);
rz=rz(rz>=rmin);
if(ne>0)
    rz=rz(end-ne+1:end);
end
rz=[0;rz];
% ensure one extrema between each root
rk=zeros(length(rz)-1,1);
for j=1:length(rz)-1
    rj=re(re>=rz(j)&re<=rz(j+1));
    if(~isempty(rj))
        [~,k]=max(abs(cf(rj)));
        rk(j)=rj(k);
    else
        rk(j)=mean(rz(j:j+1));
    end
end
re=rk; 
% ensure plot for the whole wedge
if(rz(end)<0.8)
    rz=[rz;1]; 
else
    rz(end)=1;
end

xt=[0.05,0.95];
fe=real(ufun(map(re(:))));
cs=fe(:)*linspace(xt(1),xt(2),nc);
%semilogy(abs(cf),'b'); hold on; semilogy(rz,rz*0,'.r'); semilogy(re,abs(cf(re)),'.b'); hold off; return


r=R*(diff(rz(:))*linspace(0,1,nplt)+rz(1:end-1))';
da=angle(dw(1)./dw(2))/2;
%th=linspace(tha-da,tha+da,2*nplt+1);
th=tha+da*chebpts(2*nplt+1)';
%th=tha+da*tanh(linspace(-5,5,2*nplt+1));


% Plot the eddy one by one
zz=w(2)+exp(1i*th(:))*r(:)';
zz=reshape(zz,length(th),[],numel(rz)-1);
psi=real(ufun(zz));
for e=1:size(cs,1)
    cse=cs(e,:); 
    cse=cse(abs(cse)>stol); 
    cse=cse(:); 
    ee = 1:e;
    ze=reshape(zz(:,:,ee),size(zz,1),[]);
    psie=reshape(psi(:,:,ee),size(psi,1),[]);
    %[min(abs(ze-w(2)),[],'all'),max(abs(ze-w(2)),[],'all'),max(abs(cse))]
    if(length(cse)>1)
        contour(real(ze),imag(ze),psie,cse,lw,2); hold on;
    end
end

% Plot the zero contour
cse=[f0,f0];
ee = 2:size(cs,1);
kk=floor(size(zz,2)/2):size(zz,2);
ze=reshape(zz(:,kk,ee),size(zz,1),[]);
psie=reshape(psi(:,kk,ee),size(psi,1),[]);
contour(real(ze),imag(ze),psie,cse,'k',lw,2); hold on;

if numel(re)==0
    return
end

% Add a label
if(ne==0), ne=nnz(abs(fe)>stol); end
for e=1:0
    ze=map(re(end-e+1));
    text(real(ze),imag(ze),sprintf('%1.1E',fe(end-e+1)),ha,'center',va,'middle','fontsize',16);
end
plot(w,'k',lw,2); 
hold off; axis equal;

cmap=zeros(nc+1,3);
cmap(2:nc+1,:)=magma(nc);

[~,id]=max(abs(fe));
if(fe(id)<0)
    cmap = flipud(cmap);
    ydir = 'reverse';
else
    ydir = 'normal';
end

colormap(gca,cmap);
dc = abs(cs(end,1)-cs(end,2));
zt = [min(abs(cs(end,:))), max(abs(cs(end,:)))];
zt = sign(fe(id))* (zt([1,1,2]) + [-3,-2,1]*dc/2);
caxis(sort(zt([1,end])));
xtk = sort([zt(2), cs(end,:)]);
xlb = sort([0, cs(end,:)]);
valexp = ceil(log(min(abs(cs(end,:))))/log(10));
lab = num2str(xlb(:)/10^valexp,"%.1f");
cbar=colorbar('TickLabelInterpreter', 'latex', 'Fontsize', 20, ...
                'YDir', ydir,'XTick', xtk,'XTickLabel',lab);
title(cbar,sprintf('$\\times 10^{%d}$',valexp),'Interpreter','latex');


xt=xticks();
xl=xlim();
yl=ylim();
z = (xt-xl(1))/(xl(2)-xl(1));
yt= (yl(2)-yl(1))*z + yl(1);
yticks(yt);
end