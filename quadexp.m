function []=quadexp(kmax)

p=3;
sigma=3;


err=zeros(kmax,1);
for k=1:kmax
    n=k^2;


    
    
    h=1/p;
    sj=sqrt(h:h:n)';
    sj=sj(end:-1:1);
    beta=exp(-sigma*(sqrt(n)-sj));
    
    zq=1-beta;
    wq=0.5*h*sigma*beta./sj;
    zq=[-zq(end:-1:2);zq];
    wq=wq([end:-1:2,1:end]);
    
    
    nn=p*n;
    j=(-nn-1:nn)+1/2;
    dt=sqrt(1/nn)/2;
    t=dt*j;
    s=(pi/2)*sinh(t);
    zq=tanh(s);
    wq=(pi*dt/2)*cosh(t)./cosh(s).^2;
    
    f=ones(size(zq));
    err(k)=abs((wq(:)'*f(:))-2);
end

nn=2*p*(1:kmax).^2;
hh=1./nn;

C=err(end)/hh(end);

err=max(err,eps);

ms='MarkerSize';
lw='LineWidth';
figure(1);
subplot(1,2,1);loglog(nn,err,'.-k',nn,C*hh,'k',lw,1,ms,16);title('Error in integral of 1');
subplot(1,2,2);plot(zq,wq,'.-k',lw,1,ms,16);title('Quadrature');

end