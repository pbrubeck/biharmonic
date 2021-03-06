%% Stokes flow on an infinite strip
% Pablo Brubeck, 7 Sep 2020

%%
L=3;
w=[1i+L, 1i-L, inf, -1i-L, -1i+L, inf]; % simplest unbounded domain
g=[nan+1i,nan+1i,nan,nan,nan,nan+1i];   % no-slip bcs, psi=1 on the top wall
stol = 1E-8;

% For unbounded domains, the code automatically produces a background
% solution psi_bg satisfying the bcs asymptotically as z->inf.

% psi_bg is a cubic scaled by a sigmoid function centered on the boundary 
% of the finite central portion of the domain (defined by ignoring the inf's)


%%

% solve with Mobius transform
% psi = psi_bg + (abs(z-z0)^2)*phi_hat
% In this case the error is weighted by the distance to the corners

stokes(w,g,'tol',stol);


%%

% solve without Mobius transform
% psi = psi_bg + phi_hat
% In this case the error is weighted by (abs(z-z0)^-2)*(distance to corners)

stokes(w,g,'tol',stol,'nomobius');


%%

% Solve with Mobius transform but without background solution
% psi = (abs(z-z0)^2)*phi_hat

stokes(w,g,'tol',stol,'nobg');