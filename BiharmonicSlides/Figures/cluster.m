function []=cluster(n)
ms='MarkerSize';
mc='MarkerFaceColor';
lw='LineWidth';
fs='FontSize';
ha='HorizontalAlignment';
va='VerticalAlignment';


h=1/3;
sigma=4;
t=exp(-sigma*(sqrt(n)-sqrt(1:n)'));
th=exp(-sigma*(sqrt(n)-sqrt(h:h:n)'));

w=[-1;0;-1i];
dw=w(2)-w([1,3]);

un=sum(dw./abs(dw));
un=un./abs(un);

s=zeros(length(th),2);
s(:,1)=w(2)+(w(1)-w(2))*th;
s(:,2)=w(2)+(w(3)-w(2))*th;
pol=w(2)+un(1)*t;

figure(1); set(gcf,'Renderer','Painters');
plot(w,'-k',lw,1.5); hold on;
plot(s,'.r',ms,30);
plot(pol,'.b',ms,30);

za=(s(end-1,1)+s(end-2,1))/2+0.3i*(w(2)-w(1));
zb=s(end-2,2)+0.3i*(w(3)-w(2));
zc=pol(end-1)-0.4i*un;
zd=(w(1)+w(3))*0.3/2;
text(-0.5,-0.5,'$\Omega$',fs,30,'color','k',ha,'center',va,'middle');
text(real(za),imag(za),'$\partial\Omega$',fs,30,'color','k',ha,'center',va,'top');
text(real(zb),imag(zb),'$s_k$',fs,30,'color','r',ha,'right',va,'middle');
text(real(zc),imag(zc),'$z_j$',fs,30,'color','b',ha,'right',va,'bottom');
text(real(zd),imag(zd),'$w$',fs,30,'color','k',ha,'center',va,'middle');
xlim([-1.5,1.5]); ylim([-1.5,1.5]);
hold off; axis equal; axis off; axis tight;
print('-depsc','cluster1');


h=1/3;
sigma=log(8);
t=exp(sigma*(sqrt(n)-sqrt(1:n)'));
th=exp(sigma*(sqrt(n)-sqrt(h:h:n)'));
w=[1+1i;-1+1i;-1-1i;1-1i];
un=[1;-1;-1;1];
pol=zeros(length(t),2);
s=zeros(length(th),length(w));
for k=1:length(w)
    s(:,k)=w(k)+un(k)*(th-1);
end
pol(:,1)=mean(w(1:2))-1i*diff(w(1:2)).*t;
pol(:,2)=mean(w(3:4))-1i*diff(w(3:4)).*t;
L=1.25*max(abs(real(s(:))));
bnd=reshape(w+(L-1)*real(w),[],2);
figure(2); set(gcf,'Renderer','Painters');
plot(bnd,'-k',lw,1.5); hold on;
plot(bnd([1,4]),'>k',lw,1.5,mc,'k');
plot(bnd([2,3]),'<k',lw,1.5,mc,'k');
plot(s,'.r',ms,30);
plot(pol,'.b',ms,30);

za=(s(1,2)+s(2,2))/2-0.5i*(w(2)-w(1));
zb=s(2,1)-0.5i*(w(2)-w(1));
zc=pol(2,1)+1.5;
text(real(za),imag(za),'$\partial\Omega$',fs,30,'color','k',ha,'center',va,'bottom');
text(real(zb),imag(zb),'$s_k$',fs,30,'color','r',ha,'center',va,'bottom');
text(real(zc),imag(zc),'$z_j$',fs,30,'color','b',ha,'left',va,'middle');

hold off; axis equal; axis off; axis tight;
print('-depsc','cluster2');
end