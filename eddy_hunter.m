function [] = eddy_hunter(ufun,w,ne,nplt,nc)
% EDDY_HUNTER plot Moffatt eddies
if(nargin<3), ne=0; end
if(nargin<4), nplt=64; end
if(nargin<5), nc=4; end

va='VerticalAlignment';
ha='HorizontalAlignment';

dw=w([1,3])-w(2);
un=mean(dw./abs(dw));
tha=angle(un);
R=mean(abs(dw)); 
a=R*un;


f0=real(ufun(w(2)));
f=@(t) real(ufun(w(2)+a*t)-f0);
cf=chebfun(f,[0,1],'splitting','on');


x0=[roots(cf);1];
xe=roots(diff(cf));
ne=min(ne,length(x0)-1);
if(ne>0)
    x0=x0(max(1,end-ne-1):end); 
end

tol=0;
%x0=x0(x0>tol);
%xe=xe(xe>tol);

t=linspace(0,1,nplt);
r=R*(diff(x0)*t+x0(1:end-1)*(1+0*t));
r=reshape(r',1,[]);


da=angle(dw(1)./dw(2))/2;
th=linspace(tha-da,tha+da,nplt);
zz=w(2)+r(:)*exp(1i*th);

tol=1E-15;
q=linspace(0.05,0.9,nc);
fe=f(xe); 
%be=abs(fe)>tol; xe=xe(be); fe=fe(be);
cs=f0+fe*q;
cs=[cs(:);f0];


psi=ufun(zz); 
contour(real(zz),imag(zz),real(psi),cs,'k'); hold on;
if(ne==0), ne=length(xe); end
for e=1:ne
    ze=w(2)+a*xe(end-e+1);
    text(real(ze),imag(ze),sprintf('%1.1E',fe(end-e+1)),ha,'center',va,'middle');
end
plot(w,'k'); hold off;
axis equal;


end