function [z,un,id] = polypts(w,r,x)

w=reshape(w,[],1);
x=reshape(x,1,[]);
t=linspace(-1,1,numel(x));

m=numel(w);
dw=w([2:m,1])-w;
un=-1i*dw./abs(dw);
omega=asin(abs(dw)./(2*r));
zr=r.*un;
zc=w-zr.*exp(-1i*omega);

b1=find(r==inf);
b2=find(r<inf);
zc(b1)=w(b1);
un=repmat(un,1,numel(x));
z=repmat(zc,1,numel(x));
z(b1,:)=z(b1,:)+dw(b1)*(1+x)/2;
z(b2,:)=z(b2,:)+repmat(zr(b2),1,numel(t)).*exp(1i*omega(b2)*t);
un(b2,:)=repmat(zr(b2)./r(b2),1,numel(t)).*exp(1i*omega(b2)*t);
id=repmat((1:m)',1,numel(x));
z=z(:);
un=un(:);
id=id(:);
end