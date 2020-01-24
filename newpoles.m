function [pol,zs,un,id] = newpoles(n,w,rad)
% returns n exponentially clustered poles for each vertex w
w=w(:);
if(nargin<3)
    rad=inf(size(w));
end
p=3;
m=length(w);
left=[m,1:m-1];
right=[2:m,1];

dw=w-w(right);
omega=asin(abs(dw)./(2*rad));
zr=1i*rad.*(dw./abs(dw));
zc=w-zr.*exp(-1i*omega);

dir=(w(right)-w)./abs(w(right)-w)+(w(left)-w)./abs(w(left)-w);
dir=dir.*(2*(angle((w(right)-w)./(w-w(left)))>0)-1);
b=(dir==0);
dir(b)=-1i*dw(b);
dir=dir./abs(dir);

dir=dir*min(max(real(w))-min(real(w)),max(imag(w))-min(imag(w)));

sigma=4; % FIXME
beta=-exp(-sigma*(sqrt(n)-sqrt(1:n)));
pol=repmat(w,1,n)+dir*beta; pol=pol(:);

% Sample points in the boundary
h=1/p;
betah=exp(-sigma*(sqrt(n)-sqrt(h:h:n)));
delta=abs(dir*betah);

v=zeros(m,n*p,2);
v(:,:,1)=repmat(w(right)-w,1,n*p);
v(:,:,2)=repmat(w(left)-w,1,n*p);
v=v./min(abs(v(:,:,1)),abs(v(:,:,2)));
v=v/2;

zs=repmat(w,1,n*p,2);
id=repmat((1:m)',1,n*p,2);
bnd=zeros(m,n,p,2);
for k=1:size(zs,3)
    zs(:,:,k)=zs(:,:,k)+delta.*v(:,:,k);
    bnd(:,:,:,k)=repmat(mod((1-k)+(0:m-1),m)'+1,1,n,size(bnd,3));
end

% Normals
un=-1i*(w([2:end,1])-w);
un=un./abs(un);
un=un(bnd);

% Circular boundaries
for k=1:size(zs,4)
    for j=1:size(zs,3)
        for i=1:size(zs,2)
            for l=1:size(zs,1)
                s=mod(l-k,m)+1;
                if(rad(s)<inf)
                    t=(j/(size(zs,3)))*beta(i)*(3-2*k);
                    t=sign(t)-t;
                    dz=exp(1i*omega(s)*t)*zr(s);
                    zs(l,i,j,k)=zc(s)+dz;
                    un(l,i,j,k)=dz/rad(s);
                end
            end
        end
    end
end

zs=zs(:);
un=un(:);
id=id(:);
end