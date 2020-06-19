%% wedges
% Pablo Brubeck, Jowett Walk, 19 June 2020


%%
% Here I have attempted to reproduce some of Taneda's experiments 
% in a more realistic way.

%%
% Here's a flow over a v-shaped cavity:
a = (pi/180)*(28.5/2); om = exp(1i*a);
wedge = [1i*om 0 1i/om];
wedge = 2*(wedge -1i*imag(wedge(1)));
step = [3+1i -3+1i -3 wedge 3];
zbc = @(z) nan(size(z));
ubc = {@(z) 0*z; @(z) 1-(2*imag(z)-1).^2; zbc; zbc; zbc; zbc; @(z) (1-(2*imag(z)-1).^2)};
stokes(step, ubc);


%%
% Here's a v-shaped cavity driven by a rotating cylinder:
a = (pi/180)*(28.5/2); om = exp(1i*a);
rad = -1.5/2;
w = {[1i/om, rad]; 1i*om; 0};
b = asin(sin(a)/abs(rad)); 
z0 = 1i*(abs(rad)*cos(b)+cos(a));
zbc = @(z) nan(size(z));
ubc = {@(z) 1i*(z-z0)./abs(z-z0); zbc; zbc};
stokes(w, ubc);