function [lambda] = eddy_asymptotic(w,ufun,ne)
if(nargin<3), ne=4; end

ifxtick = 1;

alpha = 0.5*angle((w(1)-w(2))/(w(3)-w(2)));
L=abs(w(4)-w(2));


sa=sin(2*alpha)/(2*alpha);
g=@(x,y) [sin(x).*cosh(y)+sa*x; cos(x).*sinh(y)+sa*y];
x=[4.21; 2.26]; 
x=NewtonRaphson(g,x);
lambda=1+([1,1i]*x)/(2*alpha);
de = cos(lambda*alpha)+cos((lambda-2)*alpha);
display(lambda);

beta = imag(lambda);
t=exp(-beta/pi);
%t = 10^floor(log(t)/log(10));


ns=32*ne;
re=exp(linspace(log(L*t^(ne)),log(L*t),ns)');
ze=w(2)+exp(1i*angle(w(4)-w(2)))*re;
moff = re(:).^(lambda);
bs = real(ufun(ze(:)));
As = [real(moff), imag(moff)];
xs = As\bs(:);
ce = xs(1)-1i*xs(2);




phi = atan(real(ce)/imag(ce));

np=128;
xp = (1+chebpts(np))/2;
yp = -ne:0;
[xp,yp] = ndgrid(xp,yp);
rho = phi + pi*(xp+yp);
rp = exp(rho/beta);
zp = w(2)+exp(1i*angle(w(4)-w(2)))*rp;

up = abs(real(ce.*rp.^lambda));
up([1,end],:)=1E-20;


lw = 'LineWidth'; fs = 'FontSize'; ha= 'HorizontalAlignment';
rlim=[exp((phi+(0.5-ne)*pi)/beta) , 1.01*L];

base=10;
%base=exp(pi/beta); phi=0;

xt = rp(end,:)/L;
if(ifxtick)
rp = base.^((beta*log(rp)-phi)/pi);
rlim=base.^((beta*log(rlim)-phi)/pi);
end


loglog(rp(:),up(:),'--r', rp(:),abs(real(ufun(zp(:)))),'k', lw,1);
xlabel('Relative distance from the corner, $r/L$');  title('Stream function');
legend("Asymptotic $|\textrm{Re}\{ C r^{\lambda}\}|$ ","Numerical $|\psi|$","Location","northwest");

xlim(rlim); 
ylim([1E-20,1E0]); 
grid on;
set(gca,'xminorgrid','off','yminorgrid','off', ...
'ytick',[1E-20,1E-15,1E-10,1E-5,1E0]);

if(ifxtick)
xt2 = floor(log(xt)/log(10));
xt1 = xt.*(10.^(-xt2));
set(gca, 'XTickLabel', num2str([xt1(:),xt2(:)],'$%.2f \\times 10^{%d}$'));
end
end

