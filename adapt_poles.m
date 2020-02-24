function [pol,z,un,id,mass,npol] = adapt_poles(npol,w,rad)
% returns n exponentially clustered poles for each vertex w
w=w(:);
if(nargin<3)
    rad=inf(size(w));
end

p=3; 
h=1/p;
sigma=4;

m=numel(w);
npol(1:m)=npol;
left=[m,1:m-1];
right=[2:m,1];

dw=zeros(m,2);
dw(:,1)=w(left)-w;
dw(:,2)=w(right)-w;

omega=asin(abs(dw(:,2))./(2*rad));
zr=-1i*rad.*(dw(:,2)./abs(dw(:,2)));
zc=w-zr.*exp(-1i*omega);

dir=sum(dw./abs(dw),2);
dir=(2*(angle(dw(:,2)./dw(:,1))>0)-1).*dir;
b=(dir==0);
dir(b)=-1i*dw(b,2);
dir=-dir./abs(dir);

if(sum(angle(dw(:,1)./dw(:,2)))>=0)
    [dmax,kmax]=max(abs(dw(:,1)));
    ww=(dw(kmax,1)/dmax)*w;
    wr=sort(real(ww)); wr=wr([1 end]);
    wi=sort(imag(ww)); wi=wi([1 end]);
    scl=max([diff(wr),diff(wi)])*sqrt(2);
else
    dst=abs((repmat(w(:),1,m)-repmat(w(:).',m,1)));
    scl=min(dst+max(dst(:))*eye(m),[],2);
    scl=scl/sqrt(2);
end
dir=dir.*scl;

tol=eps;
pol=[];z=[];un=[];id=[];mass=[];
for i=1:m
    n=npol(i);
    beta=-exp(-sigma*(sqrt(n)-sqrt(1:n)'));
    beta=beta(abs(beta)>tol);
    n=length(beta);
    npol(i)=n;
    poli=w(i)+dir(i)*beta;

    % Sample points at the boundary
    sj=sqrt(h/2:h:n)';
    zq=0.5*exp(-sigma*(sqrt(n)-sj));
    wq=0.5*h*sigma*zq./sj; wq(end)=wq(end)/2;
    wq=zq;

    kd=abs(zq)>tol;
    zq=zq(kd);
    wq=wq(kd);
    np=length(zq);
    
    zi=repmat(dw(i,:),np,1);
    zi(:,1)=w(i)+zq.*zi(:,1);
    zi(:,2)=w(i)+zq.*zi(:,2);

    wi=repmat(abs(dw(i,:)),np,1);
    wi(:,1)=wq.*wi(:,1);
    wi(:,2)=wq.*wi(:,2);

    % Normals
    uni=-1i*dw(i,:);
    uni=uni./abs(uni);
    uni=[repmat(-uni(1),np,1),repmat(uni(2),np,1)];
    
    % Circular boundaries
    [zi(:,1),uni(:,1)]=circ(rad(left(i)),w(left(i)),w(i),zi(:,1),uni(:,1));
    [zi(:,2),uni(:,2)]=circ(rad(i),w(i),w(right(i)),zi(:,2),uni(:,2));
    
    pol=[pol;poli(:)];
    z=[z;zi(:)];
    un=[un;uni(:)];
    id=[id;repmat(i,numel(zi),1)];
    mass=[mass;wi(:)];
end

function [z,un]=circ(r,a,b,z,un)
    if(~isinf(r))
        da=asin(abs(b-a)/(2*r));
        R=(b-a)/(2i*sin(da));
        c=b-R*exp(1i*da);
        t=2*(z-a)/(b-a)-1;
        un=(R/r)*exp(1i*da*t);
        z=c+r*un;
    end
end
end