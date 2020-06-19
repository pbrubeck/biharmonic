%% Stokesplay1
% Nick Trefethen, Iffley, 19 June 2020

%%
% Pablo has sent me stokes.m and lap.m codes, and I am
% very happily playing with them. I won't list them here in this
% draft version as this will no doubt evolve fast.

%% 
% Here's a driven cavity.   (Maybe there's a sign error in
% the stream function?)
stokes([1+1i -1+1i -1-1i 1-1i],[1 nan 0 0]);

%%
% If we put in a fake vertex where none is needed,
% some poles accumulate there much as at the other
% lower corners, where the singularities are very weak.
% I'm not sure this is a problem in practice, though
% it's a shame cosmetically since one might imagine
% that stokes.m would reveal the presence or absence of singularities.
% One sees similar effects with laplace.m, but I think
% not as pronounced.
stokes([1+1i -1+1i -1-1i -1i 1-1i],[1 nan 0 0 0]);

%%
% Here we have a cavity on the cavity:
stokes([1+1i -1+1i -1-1i -.5-1i -.5-2i .5-2i .5-1i 1-1i],[1 nan nan nan nan nan 0 0]);

%% 
% Here's the beautiful triangular wedge.
a = (pi/180)*(28.5/2); om = exp(1i*a);
stokes([1i/om 1i*om 0],[1 nan 0]);

%%
% And here's the flow over a step.
step = [3+1i -3+1i -3 0 -1i 3-1i];
zbc = @(z) nan(size(z));
ubc = {@(z) 0*z; @(z) 1-(2*imag(z)-1).^2; zbc; zbc; zbc; @(z) (1-imag(z).^2)/2};
stokes(step,ubc);

%%
% Adding in a fake vertex makes a modest difference.
step = [3+1i; 1i; -3+1i; -3; 0; -1i; 3-1i];
ubc = {@(z) 0*z; @(z) 0*z; @(z) 1-(2*imag(z)-1).^2; zbc; zbc; zbc; @(z) (1-imag(z).^2)/2};
stokes(step,ubc);

%%
% Here's a gentler step.
step = [3+1i -3+1i -3 -1 -1i 3-1i];
ubc = {@(z) 0*z; @(z) 1-(2*imag(z)-1).^2; zbc; zbc; zbc; @(z) (1-imag(z).^2)/2};
stokes(step,ubc);

%%
% Here's a flow over a cavity:
step = [3+1i -3+1i -3 -1 -1-1i 1-1i 1 3];
zbc = @(z) nan(size(z));
ubc = {@(z) 0*z; @(z) 1-(2*imag(z)-1).^2; zbc; zbc; zbc; zbc; zbc; @(z) (1-(2*imag(z)-1).^2)};
stokes(step,ubc);

%%
% Some possible adjustments to the code:
%
% 1.  As mentioned, maybe there is a sign error in the stream function.
%
% 2.  We could get rid of the black dot marking wc.  That made sense
% in the pre-Arnoldi days, but now wc is much less significant.
%
% 3.  The help text still talks about laplace.
%
% 4.  If a user just specifies velocities, the picture should
% come out with reasonable stream function values.  I think there
% should be a default procedure for such cases: e.g., set
% the minimum stream function value to zero.
%
% 5.  Perhaps we should make poles accumulate more slowly at
% corners that are nonsingular or nearly nonsingular.  (Stokes
% problems tend to have a lot of very weak singularities,
% I suspect.)
% Or indeed, perhaps we should add something to this code (and
% likewise to laplace) to estimate a posteriori how strong the
% singularities are?
%
% 6. This last idea is not far from an a posteriori application of aaa,
% something I've not pursued much but should perhaps return to thinking
% about.

