function [pol,z,un,id,mass] = adapt_poles(npol,w,rad)
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
dir=dir*min(max(real(w))-min(real(w)),max(imag(w))-min(imag(w)));

pol=[];z=[];un=[];id=[];mass=[];
for i=1:m
    n=npol(i);
    beta=-exp(-sigma*(sqrt(n)-sqrt(1:n)'));
    beta=beta(abs(beta)>1E-15);
    n=length(beta);
    poli=w(i)+dir(i)*beta;

    % Sample points at the boundary
    zq=0.5*exp(-sigma*(sqrt(n)-sqrt(h:h:n)'));
    wq=0.5*h*sigma*zq./sqrt(h:h:n)';
    wq(end)=wq(end)/2;

    
    nn=p*n;
    j=(0:nn)'+1/2;
    dt=sqrt(1/nn)/2;
    t=dt*j;
    s=(pi/2)*sinh(t);
    %zq=(1-tanh(s))/2;
    %wq=(dt*pi/4)*cosh(t)./cosh(s).^2;


    kd=abs(zq)>1E-15;
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
    
    pol=[pol;poli(:)];
    z=[z;zi(:)];
    un=[un;uni(:)];
    id=[id;repmat(i,numel(zi),1)];
    mass=[mass;wi(:)];
end

% Circular boundaries % TODO
% for k=1:size(zs,4)
%     for j=1:size(zs,3)
%         for i=1:size(zs,2)
%             for l=1:size(zs,1)
%                 s=mod(l-k,m)+1;
%                 if(rad(s)<inf)
%                     t=(j/(size(zs,3)))*beta(i)*(3-2*k);
%                     t=sign(t)-t;
%                     dz=exp(1i*omega(s)*t)*zr(s);
%                     zs(l,i,j,k)=zc(s)+dz;
%                     un(l,i,j,k)=dz/rad(s);
%                 end
%             end
%         end
%     end
% end
end