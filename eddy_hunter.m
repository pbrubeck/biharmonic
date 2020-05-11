function [cs] = eddy_hunter(ufun,w,ne,nplt,nc,stol)
% EDDY_HUNTER plot Moffatt eddies
if(nargin<3), ne=0; end
if(nargin<4), nplt=64; end
if(nargin<5), nc=4; end
if(nargin<6), stol=1E-11; end

va='VerticalAlignment';
ha='HorizontalAlignment';

dw=w([1,3])-w(2);
un=mean(dw./abs(dw));
%un=un/abs(un);
tha=angle(un);
R=mean(abs(dw)); 
a=R*un;


f0=real(ufun(w(2))); 
f0
%f0=0;
cf=chebfun(@(x) real(ufun(w(2)+a*x)-f0), [0,1], 64); % low frequency
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
re=rk; fe=cf(re);  
% ensure plot for the whole wedge
if(rz(end)<0.8)
    rz=[rz;1]; 
else
    rz(end)=1;
end



fe=fe(:)+f0;
fe=real(ufun(w(2)+a*re(:)));
cs=fe(:)*linspace(0.05,0.95,nc);
%plot(cf,'b'); hold on; plot(rz,rz*0,'.r'); plot(re,cf(re),'.b'); hold off; return


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
        contour(real(ze),imag(ze),psie,cse,'k'); hold on;
    end
end

% Plot the zero contour
cse=[f0,f0];
ee = 2:size(cs,1);
ze=reshape(zz(:,:,ee),size(zz,1),[]);
psie=reshape(psi(:,:,ee),size(psi,1),[]);
contour(real(ze),imag(ze),psie,cse,'k'); hold on;



% Add a label
if(ne==0), ne=nnz(abs(fe)>stol); end
for e=1:1
    ze=w(2)+a*re(end-e+1);
    text(real(ze),imag(ze),sprintf('%1.1E',fe(end-e+1)),ha,'center',va,'middle');
end
plot(w,'k'); 
hold off; axis equal;

xt=xticks();
xl=xlim();
yl=ylim();
z = (xt-xl(1))/(xl(2)-xl(1));
yt=(yl(2)-yl(1))*z + yl(1);
yticks(yt);
end