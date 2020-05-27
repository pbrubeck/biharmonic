function [lambda] = eddy_asymptotic(w,ufun,ne)
if(nargin<3), ne=4; end
alpha = 0.5*angle((w(1)-w(2))/(w(3)-w(2)));

sa=sin(2*alpha)/(2*alpha);
g=@(x,y) [sin(x).*cosh(y)+sa*x; cos(x).*sinh(y)+sa*y];
x=[4.21; 2.26]; 
x=NewtonRaphson(g,x);
lambda=1+([1,1i]*x)/(2*alpha);
de = cos(lambda*alpha)+cos((lambda-2)*alpha);

beta = imag(lambda);
L=abs(w(4)-w(2));

ns=128;
re=exp(linspace(log(1E-4*L),log(1E-1*L),ns)');
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
rlim=[exp((phi+(0.5-ne)*pi)/beta) ,L];

base=10;
%base=exp(pi/beta); phi=0;
if(false)
rp = base.^((beta*log(rp)-phi)/pi);
rlim=base.^((beta*log(rlim)-phi)/pi);
end

loglog(rp(:),up(:),'--r', rp(:),abs(real(ufun(zp(:)))),'k', lw,1);
xlabel('Distance from the corner, $r$');  title('Stream function');
legend("Asymptotic $|\textrm{Re}\{ C r^{\lambda}\}|$ ","Numerical $|\psi|$","Location","northwest");

xlim(rlim); 
ylim([1E-20,1E0]); 
grid on;
set(gca,'xminorgrid','off','yminorgrid','off', ...
        'ytick',[1E-20,1E-15,1E-10,1E-5,1E0]);

return
xt = 1-ne:0;
set(gca,'xminorgrid','off','yminorgrid','off','XTickLabel', num2str(xt(:),'$b^{%d}$'));
text(1/sqrt(base),1E-6,sprintf('$b=%.2f$',exp(pi/beta)),fs,20,ha,'center');
end

