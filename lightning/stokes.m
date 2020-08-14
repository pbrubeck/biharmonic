function [u, maxerr, f, Z, Zplot, A, pol] = stokes(P, varargin) 
%STOKES  Lightning Stokes solver.
%         U = STOKES(P,G) solves the Stokes problem via the stream function
%         formulation with slip or no-slip boundary data on the simply-connected
%         region Omega bounded by P, which may be a polygon or circular polygon,
%         with computation also of the harmonic conjugate.
%         v1, (c) Pablo D. Brubeck, U. of Oxford, June 2020. 
%         For info see https://people.maths.ox.ac.uk/trefethen/lightning.html.
%         Please send corrections or suggestions to brubeckmarti@maths.ox.ac.uk.
%
%  Inputs:
%      P = vector of corners as complex numbers z = x+iy in counterclockwise
%              order to specify a polygon
%          or cell array of corners v and pairs [v r] to specify a circular
%              polygon: r = radius of curvature of arc from this v to the next
%          or 'pent'[agon], 'snow'[flake], 'iso'[spectral], 'L', or 'circleL'
%          or integer >= 3, the number of corners of a random polygon
%          or integer <= -3, -1 x no. of corners of a random circular polygon
%
%      g = function handle for complex-valued velocity boundary data u+1i*v
%          or cell array of function handles for sides P1-P2, P2-P3,...
%          or vector of constant values for these sides
%          (default @(z) real(z).^2)).  If g = nan+1i*C_0 on any side, the 
%          BC imposed there is no-slip, i.e., u_tangential = 0, psi = C_0.
%
%  Further optional inputs:
%    'tol', tol = tolerance for maximal absolute error (default 1e-6)
%    'noplots' to suppress plotting
%    'steps' for step-by-step plots of errors on boundary and poles
%    'rel' to weight error by scaled dist. to corner (automatic if g discont.)
%    'slow' to turn off adaptive mode for cleaner root-exp convergence curves
%    'noarnoldi' to turn off Arnoldi stabilization of polynomial term
%    'aaa' for AAA compression of result  (to be supported).  Requires Chebfun in path.
%
%  Outputs for [U,MAXERR,F,Z,ZPLOT,A] = stokes(P,G)
%       u = function handle for the velocity u(z)+i v(z)=2i*d/dz psi 
%       associated with the solution of lap^2(psi) = 0, u+iv = g, on boundary
%  maxerr = upper bound on maximal error, even near singularities at corners
%             (or error weighted by distances to corners, if 'rel' is specified)
%       f = function handle for biharmonic function f = psi - iA
%       Z = sample points on boundary
%       Zplot = set of points for plotting boundary
%       A = final rectangular matrix used for least-squares problem
%
% Reference: P. Brubeck and L. N. Trefethen, Lightning Stokes solver.
%
% Examples:
%
%   stokes('ldc');                         % lid driven cavity
%   stokes('step');                        % flow over a step
%   stokes([0 1 1+1i 1i],[0 1i 0 0]);      % rotated lid driven cavity
%   stokes('infchan');                     % infinite channel
%   stokes([3+1i -3+1i inf -3 0 -1i 3-1i inf],[nan+2i/3 nan+2i/3 0 0 0 0 0 nan+2i/3]);  % same as above

%% Set up the problem
[g, w, ww, pt, dw, tol, steps, plots, ...        % parse inputs
    slow, rel, arnoldi, aaaflag, mobflag, ubkg, fbkg] = ...
    parseinputs(P,varargin{:});
Zplot = ww;                                        
nw = length(w);                                  % number of corners
wr = sort(real(ww)); wr = wr([1 end]);
wi = sort(imag(ww)); wi = wi([1 end]);
wc = mean(wr+1i*wi);                             % for scale- and transl-invariance
scl = max([diff(wr),diff(wi)]);
q = .5; if slow == 1, q = 0; end                 % sets which corners get more poles
inpolygonc = @(z,w) inpolygon(real(z), ...
            imag(z),real(w),imag(w));            % complex variant of "inpolygon"
        
sig = zeros(nw,1); sig(:)=4;
outward = zeros(nw,1);
corners = find(~isinf(w))';

for k = corners
   kn = mod(k,nw)+1;
   kp = mod(k-2,nw)+1;
   if(isinf(w(kp)) || isinf(w(kn)))
       j = k; L = dw(j); 
       forward = pt{j}(.51*L) - pt{j}(.5*L);            % small step toward next corner                            
       backward = pt{j}(.49*L) - pt{j}(.5*L);           % small step toward last corner
       tmp = 1i*backward*sqrt(-forward/backward);
       outward(k) = (scl)*tmp/abs(tmp);
       sig(k) = -sig(k);
   else
       forward = pt{k}(.01*dw(k)) - w(k);            % small step toward next corner
       j = mod(k-2,nw)+1;                             
       backward = pt{j}(.99*dw(j)) - w(k);           % small step toward last corner
       tmp = 1i*backward*sqrt(-forward/backward);
       outward(k) = tmp/abs(tmp);                    % outward direction from corner
   end
end
warn = warning('off','MATLAB:rankDeficientMatrix');  % matrices are ill-conditioned

if(mobflag) % move wc outside omega
    j = find(sig<0 & dw>0,1); L = dw(j);
    tmp = pt{j}(.51*L) - pt{j}(.49*L); tmp=1i*tmp/abs(tmp);
    wc = wc - max(2*real(conj(tmp)*(wc-pt{j}(.5*L))),scl/2)*tmp ;
end

%% Set up for plots
if plots
   LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
   fs = 9; PO = 'position'; FW = 'fontweight'; NO = 'normal';
   sx = linspace(wr(1),wr(2),100); sy = linspace(wi(1),wi(2),100);
   [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
   ax = [wr(1:2); wi(1:2)] + .2*scl*[-1 1 -1 1]';
   axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';
end
if steps
   clf, subplot(1,2,1), shg, plot(w([1:end 1]),'k',LW,1)
   grid on, axis equal, axis(axwide)
end

%% Main loop: increase number of poles until convergence ==============================
Nvec = []; errvec = []; tic
errk = ones(nw,1);                               % max error near each corner
nkv = zeros(nw,1);                               % no. of poles at each corner
maxstepno = 30; err0 = Inf;
for stepno = 1:maxstepno
   % Fix poles and sample pts on bndry.  Side k means side from corner k to k+1.
   Z = [];           % col vector of sample points on boundary
   G = [];           % col vector of boundary values at these points
   T = [];           % col vector of unit tangent vectors at these points
   pol = [];         % row vector of poles of the rational approximation
   d = [];           % row vector of distances from poles to their corners
   tt = cell(nw,1);  % cell array of distances of sample points along each side

   for k = 1:nw
      nk = nkv(k); % no. of poles at this corner
      dk = sqrt(1:nk) - sqrt(nk);
      dk = scl*exp(abs(sig(k))*dk);             % stronger clustering near corner
      dk = dk(dk>1e-15*scl);                    % remove poles too close to corner
      if(~isinf(w(k)))
          if(sig(k)>0)
            polk = w(k) + outward(k)*dk;    % poles near this corner
          else
            polk = w(k) + outward(k)./dk;      
          end
          ii = find(inpolygonc(polk(dk>1e-12*scl),ww),1); % work around inaccuracy
          if ~isempty(ii)                               % don't allow poles in Omega
              dk = dk(1:ii-2); polk = polk(1:ii-2);
          end
          pol = [pol polk]; d = [d dk];
      end
      dvec = [(1/3)*dk (2/3)*dk dk];                % finer pts for bndry sampling
      j = mod(k-2,nw)+1;
      if(~isinf(w(j)) && ~isinf(w(mod(k,nw)+1)))
      tt{k} = [tt{k} dvec(dvec<dw(k))];             % add clustered pts near corner
      tt{j} = [tt{j} dw(j)-dvec(dvec<dw(j))];       % likewise in other direction
      end
   end

   for k = 1:nw
      ntt = max(30, max(nkv));

      tt{k} = [tt{k} (dw(k)/(ntt))*((1:ntt)-1/2)];      % additional pts along side
      tt{k} = sort(tt{k}(:));
      tk = tt{k}; pk = pt{k};                       % abbrevations 
      Z = [Z; pk(tk)];                              % sample pts on side k
      G = [G; g{k}(pk(tk))];                        % boundary data at these pts
      h = 1e-4;                                     % 4-pt trapezoidal rule
      T = [T; (pk(tk+h)-1i*pk(tk+1i*h) ...
             - pk(tk-h)+1i*pk(tk-1i*h))/(4*h);];    % unnormalized tangent vector
   end
   T = T./abs(T);                                   % normalize tangent vectors
   II = isnan(G);                                   % Neumann indices

   % Solve the Stokes problem
   n = 4*stepno;                                    % degree of polynomial term
   
   Np = length(pol);
   [wt,Kj] = build_wt(Z,w,wc,scl,rel,mobflag);
   [A,H] = build_ls(n,Z,wc,pol,d,scl,T,II,arnoldi,mobflag,wt);
   [M,N] = size(A);                                  % no. of cols = 2n+1+2Np
   
   Gn = G; Gn(II) = -1i*imag(Gn(II));  % set Neumann vals to 0
   Gn = [real(Gn); -imag(Gn); zeros(M-2*length(Gn),1)];
   if(~any(II))
      Gn(end-1)=-(wt'*fbkg(Z))*(1/sum(wt));
   end

   wtt = ones(M,1);
   wtt(1:2*length(wt)) = [wt;wt];
   W = spdiags(sqrt(wtt),0,M,M);      % weighting for case 'rel'
   WA = W*A;
   aa = sqrt(sum(WA.^2,1));  % column scaling
   %aa = max(abs(WA),[],1);
   ja = find(aa>0);
   PA = spdiags(reshape(1./aa(ja),[],1),0,numel(ja),numel(ja));
   c = zeros(size(A,2),1);
   c(ja) = PA*((WA(:,ja)*PA)\(W*Gn));
   cc = c(1:end/2) + 1i*c(end/2+1:end);
   f = @(fld,z) reshape(fzeval(fld,z(:),wc,...       % vector and matrix inputs
              cc,H,pol,d,arnoldi,mobflag,scl,n),size(z));    % to u and f both allowed
   psi= @(z) real(f(0,z));
   u  = @(z) f(1,z);
   for k = 1:nw
      K1 = find(Kj==k);
      K2 = K1+length(Z);
      errk(k) = norm(wtt(K1).*((A(K1,:)*c-Gn(K1))+1i*(A(K2,:)*c-Gn(K2))),inf); % error near corner k
   end
   err = norm(errk,inf);                     % global error

   ii = find(isinf(w));
   ii = unique([mod(ii-2,nw)+1, mod(ii,nw)+1]);
   %errk(ii) = max(errk(ii));
   
   polmax = 100;
   for k = 1:nw
      if (errk(k) > q*err) && (nkv(k) < polmax)
          nkv(k) = nkv(k)+ceil(1+sqrt(nkv(k)));      % increase no. poles
      else
          nkv(k) = max(nkv(k),ceil(stepno/2));
         %nkv(k) = min(polmax,nkv(k)+1);
      end  
   end
   %nkv(ii) = min(nkv(ii));
   nkv(isinf(w)) = max(nkv(ii));
   nkv
   
   if steps                                           % plot error on bndry
      subplot(1,2,1), plot(ww,'k',LW,1), grid on
      axis equal, axis(axwide), hold on
      plot(pol,'.r',MS,7), hold off
      subplot(1,2,2), semilogy([-pi pi],err*[1 1],'--b',LW,1)
      hold on, axis([-pi pi 1e-16 100]), grid on
      semilogy(angle(Z-wc),wt.*abs(u(Z)-G),'.k',MS,4)
      semilogy(angle(pol-wc),d,'.r',MS,7), hold off
      set(gca,'ytick',10.^(-16:4:0))
      set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
      set(gca,FS,fs-1), xlabel('angle on boundary wrt wc',FS,fs)
      title('bndry err (black) & poles (red)',FS,fs)
      disp('Press <enter> for next plot'); pause
   end
   errvec = [errvec err]; Nvec = [Nvec; N];
   if err < .5*tol, break, end                        % convergence success
   if err < err0                                      % save the best so far
      u0 = u; f0 = f; Z0 = Z; G0 = G; A0 = A; M0 = M;
      N0 = N; err0 = err; pol0 = pol; wt0 = wt;
   end
   if (N > 1200) || (stepno == maxstepno) || (Np >= polmax*nw)  % failure
      u = u0; f = f0; Z = Z0; G = G0; A = A0; M = M0;
      N = N0; err = err0; pol = pol0; wt = wt0;
      warning('STOKES failure.  Loosen tolerance or add corners?')
      break
   end
end
warning(warn.state,'MATLAB:rankDeficientMatrix')       % back to original state
tsolve = toc;  % =========== end of main loop =========================================

aaaflag = 0;
% Compress with AAA approximation (requires Chebfun aaa in the path)
if aaaflag
   [faaa,polaaa] = aaa(f(Z),Z,'mmax',N/2,'tol',0.1*tol,'lawson',0,'cleanup',0);
   if length(find(inpolygonc(polaaa,ww))) == 0    % AAA successful
      f = faaa; pol = polaaa;
      u = @(z) real(f(z));
   else                                           % AAA unsuccess.: pole in region
     %badpol = polaaa(find(inpolygonc(polaaa,ww)));% poles in polygon
      warning('STOKES AAA compression unsuccessful; returning uncompressed solution.')
   end
end

%% Finer mesh for a posteriori error check
Z2 = []; G2 = [];
for k = 1:nw
   newtt = mean([tt{k}(1:end-1) tt{k}(2:end)],2);
   newpts = pt{k}(newtt);
   Z2 = [Z2; newpts];
   G2 = [G2; g{k}(newpts)];
end
uu2 = u(Z2);
wt2 = build_wt(Z2,w,wc,scl,rel,mobflag);
err2 = norm(wt2.*(G2-uu2),inf);
maxerr = max(err,err2);   

%% Convergence curve plot
if plots
   ws = 'error'; if rel, ws = 'weighted error'; end
   if steps, figure, else clf, end, shg
   axes(PO,[.09 .65 .35 .26])
   semilogy(sqrt(Nvec),errvec,'.-k',LW,0.7,MS,10), grid on, hold on
   semilogy(sqrt(N),maxerr,'or',MS,7,LW,1), hold off
   errmin = .01*tol; axis([0 1.1*max(sqrt(Nvec)) 1e-14 100])
   set(gca,FS,fs-1), title('convergence',FS,fs,FW,NO)
   if arnoldi == 0, title('convergence - no Arnoldi',FS,fs,FW,NO), end
   xlabel('sqrt(DoF)',FS,fs), ylabel(ws,FS,fs)
   set(gca,'ytick',10.^(-16:4:0))
   ax2 = axis; x1 = ax2(1) + .05*diff(ax2(1:2));
   s = sprintf('solve time = %6.3f secs',tsolve);
   if ~steps, text(x1,4e-11,s,FS,fs), end
   z = randn(1000,1)+1i*randn(1000,1); z = z/10;
   tic, u(z); teval = 1e3*toc;
   s = sprintf('eval time = %4.1f microsecs per pt',teval);
   text(x1,4e-13,s,FS,fs)
end

%% Contour plot of solution
uz = u(Z);
u = @(z) u(z)+ubkg(z); psi = @(z) psi(z)+fbkg(z);   % add background solution
if plots
   uu = real(psi(zz)); uu(~inpolygonc(zz,ww)) = nan;
   uu = uu - min(uu(:));
   axes(PO,[.52 .34 .47 .56]), levels = linspace(min(uu(:)),max(uu(:)),21);
   
   umax = max(uu(:)); umin = min(uu(:));
   levels = [];
   if(umax>0), levels = [levels, linspace(0, umax, 5+16*(umax>-100*umin))]; end
   if(umin<0), levels = [levels, linspace(umin, 0, 5+16*(-umin>100*umax))]; end
   
   contour(sx,sy,uu,levels,LW,.7), colorbar, axis equal, hold on
   plot(ww,'-k',LW,1), plot(pol,'.r',MS,6)
   set(gca,FS,fs-1); axis(ax); plot(real(wc),imag(wc),'.k',MS,6);
   
   %plot(real(Z),imag(Z),'.k',MS,6);
   title(['dim(A) = ' int2str(M) ' x ' int2str(N) ' ', ...
       ' #poles = ' int2str(length(pol))],FS,fs,FW,NO), hold off
end

%% Error plot along boundary
if plots
   wc = mean(wr+1i*wi);
   axes(PO,[.09 .21 .35 .28])
   semilogy([-pi pi],maxerr*[1 1],'--b',LW,1), hold on
   semilogy(angle(Z2-wc),wt2.*abs(uu2-G2),'.r',MS,4)
   axis([-pi pi .0001*errmin 1]), grid on
   semilogy(angle(Z-wc),wt.*abs(uz-G),'.k',MS,4), hold off
   set(gca,'ytick',10.^(-16:4:0))
   set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
   set(gca,FS,fs-1), xlabel('angle on boundary wrt wc',FS,fs)
   title([ws ' on boundary'],FS,fs,FW,NO)
end

%% Put the logo in the corner
if plots
   axes(PO,[.82 .14 .12 .12]), lightninglogo
end


return
clf;
sx = 2*linspace(wr(1),wr(2),100); sy = linspace(wi(1),wi(2),100);
[xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
uu = real(fbkg(zz)); uu(xx<0 & yy<0) = nan;
surf(sx,sy,uu);

end   % end of main program

function [fZ] = fzeval(fld,Z,wc,cc,H,pol,d,arnoldi,mobflag,scl,n) 
    ZW = Z-wc; if mobflag, ZW = 1./ZW; end
    if arnoldi
       [Q0, Q1] = arnoldi_eval(ZW,H);
    else
       Q0 = (ZW/scl).^(0:n);
       Q1 = ((0:n)/scl).*(ZW/scl).^[0 0:n-1];
    end
    if mobflag, Q1 = Q1.*(-ZW.^2); end
    R0 = [Q0  d./(Z-pol)];
    R1 = [Q1 -d./(Z-pol).^2];
    i1 = 1:size(R0,2);
    i2 = i1(end)+(1:size(R0,2));
    if (fld==0)    % stream function
        fZ = -1i*(conj(Z-wc).*(R0*cc(i1)) +  R0*cc(i2));
    elseif(fld==1) % velocity u+1i*v
        fZ = conj(conj(Z-wc).*(R1*cc(i1)) - conj(R0*cc(i1)) + (R1*cc(i2)) );
    else % pressure
        fZ = conj(4*R1*cc(i1));
    end
end

function [A,H] = build_ls(n,Z,wc,pol,d,scl,T,II,arnoldi,mobflag,wt)
    H = zeros(n+1,n);    % Arnoldi Hessenberg matrix
    
    ZW = Z-wc; if mobflag, ZW = 1./ZW; end
    if arnoldi               
        [H, P0, P1] = arnoldi_fit(ZW,n); 
    else                                             % no-Arnoldi option
        P0 = (ZW/scl).^(0:n);                     % (for educational purposes)
        P1 = ((0:n)/scl).*(ZW/scl).^[0 0:n-1];
    end
    if mobflag, P1 = P1.*(-ZW.^2); end
    R0 = [P0  d./(Z-pol)];
    R1 = [P1 -d./(Z-pol).^2];
    
    ZZ = spdiags(conj(Z-wc),0,numel(Z),numel(Z));
    % slip boundary condition        [u; v]
    F1 = [ZZ*R1-R0, R1];
    F2 = [ZZ*R1+R0, R1];
    A = [real(F1) -imag(F1); imag(F2) real(F2)];
    % no-slip boundary condition     [u_T; psi]
    if any(II)
        ni = nnz(II);
        ZZ = spdiags(conj(Z(II)-wc),0,ni,ni);
        TT = spdiags(1i*T(II),0,ni,ni);
        F1 = [(TT*ZZ)*R1(II,:)+conj(TT)*R0(II,:), TT*R1(II,:)];
        F2 = [ZZ*R0(II,:), R0(II,:)];
        A([II,II],:) = [imag(F1) real(F1); imag(F2) real(F2)];
    else
        F2 = [ZZ*R0, R0];
        c = (wt'*F2)*(1/sum(wt));
        A = [A; imag(c) real(c)];
    end
    wp = (wt'*Z)/sum(wt); 
    zw = wp-wc; if mobflag, zw=1./zw; end
    if arnoldi
        [q0, q1] = arnoldi_eval(zw,H);
    else
        q1 = ((0:n)/scl).*(zw/scl).^[0 0:n-1];
    end
    
    if mobflag, q1=q1.*(-zw.^2); end
    p0 = [q1, -d./(wp-pol).^2, zeros(1,size(R0,2))];
    A = [A; real(p0) -imag(p0)];
end

function [wt,Kj] = build_wt(Z,w,wc,scl,rel,mobflag) % weights to measure error
    [~,Kj]=min(abs(Z-reshape(w,1,[])),[],2);
    if rel
        wt = abs(Z-w(Kj))/scl;
        if mobflag
            wm = abs((w(Kj)-wc)./(Z-wc));           
            rt = min(abs(diff(w([1:end,1]))))/scl;
            ii = wt>rt;
            wt(ii) = min(rt, rt*wm(ii));
            Kj(ii) = find(isinf(w),1);
        end
    else
        wt = ones(length(Z),1);
    end
end

function [H,P,D] = arnoldi_fit(z,n) % polyfit with Arnoldi
    M = length(z); 
    H = zeros(n+1,n); 
    P = ones(M,n+1); 
    D = zeros(size(P)); 
    for k = 1:n       
        p = z(:).*P(:,k);
        for j = 1:k
            H(j,k) = P(:,j)'*p;
            p = p - H(j,k)*P(:,j);
        end 
        for j = 1:k     % Gram-Schmidt twice
            DelH = P(:,j)'*p;
            p = p - DelH*P(:,j);
            H(j,k) = H(j,k)+DelH;
        end
        H(k+1,k) = norm(p);
        P(:,k+1) = p*(1/H(k+1,k));

        d=z(:).*D(:,k)+P(:,k);
        d=d-D(:,1:k)*H(1:k,k);
        D(:,k+1)=d*(1/H(k+1,k));
    end
end

function [P,D] = arnoldi_eval(z,H)  % polyval with Arnoldi
    M=numel(z); n=size(H,2); 
    P=ones(M,n+1);
    D=zeros(M,n+1);
    for k=1:n
        p=z(:).*P(:,k);
        p=p-P(:,1:k)*H(1:k,k);
        P(:,k+1)=p/H(k+1,k);

        d=z(:).*D(:,k)+P(:,k);
        d=d-D(:,1:k)*H(1:k,k);
        D(:,k+1)=d/H(k+1,k);
    end
end

function lightninglogo     % plot the lightning Laplace logo
    s = linspace(0,1,40)';
    v = exp(-.35i)*[0 1+2i .5+1.85i 1+3i 0+2.7i -.2+1.3i .1+1.4i];
    w = v(7) + (v(1)-v(7))*s;
    for k = 1:6; w = [w; v(k)+(v(k+1)-v(k))*s]; end
    w = w + .05*imag(w).^2;
    fill(real(w),imag(w),[1 1 .5]), axis equal, hold on
    plot(w,'-k','linewidth',.7)
    dots = .85*(v(3)+.01)*.72.^(0:5);
    dots = dots + .05*imag(dots).^2;
    for k = 1:6, plot(dots(k),'.r','markersize',13-2*k), end
    hold off, axis off
end

function [g, w, ww, pt, dw, tol, steps, plots, slow, ...
          rel, arnoldi, aaaflag, mobflag, ubkg, fbkg] = parseinputs(P,varargin)

%% Defaults
tol = 1e-6; steps = 0; plots = 1;
slow = 0; rel = 0; aaaflag = 0; arnoldi = 1;

%% First treat the domain, defined by P
demoflag = 0;
randomcirc = 0;
if ~iscell(P)                                       
   if isnumeric(P)
      if length(P) > 1, w = P;                    % vertices have been specified
      else
         if P < 0
            randomcirc = 1; P = -P;               % random circular arcs
         end
         w = exp(2i*pi*(1:P)/P).*(.1+rand(1,P));  % random vertices
      end
   else
      if strcmp(P,'L'), w = [2 2+1i 1+1i 1+2i 2i 0];
      elseif strcmp(P,'circleL'), P = {2 [2+1i -1] 1+2i 2i 0};
      elseif strcmp(P,'pent'), w = .7*exp(pi*2i*(1:5)/5);
      elseif strcmp(P,'snow'), P = exp(2i*pi*(1:12)/12);
                               w = P.*(1+.2*(-1).^(1:12)); w = w/1.4;
      elseif strcmp(P,'iso')
         w = [1+2i 1+3i 2i 1i+1 2+1i 2 3+1i 3+2i]-(1.5+1.5i); w = w/1.8;
      elseif strcmp(P,'ldc')
         w = [1+1i -1+1i -1-1i 1-1i]; demoflag=1;
      elseif strcmp(P,'step')
         w = [3+1i -3+1i -3 0 -1i 3-1i]; demoflag=2;
      elseif strcmp(P,'infchan')
         w = [3+1i 1i -3+1i inf -3 0 -1i 3-1i inf]; demoflag=3;
      end
   end
   if ~iscell(P), P = num2cell(w); end            % convert to cell array
   if randomcirc
      for k = 1:length(P)
         r = .6/rand;
         P{k} = [P{k} r*(-1)^double(randn>0)];
      end
   end
end

nw = length(P);
for k = 1:nw
    w(k) = P{k}(1);
end

w = w(:);
mobflag = any(isinf(w));
wf = w(~isinf(w));

nw = length(w);
pt = cell(nw,1);
dw = zeros(nw,1);
ww = [];          % bndry pts for plotting
for k = 1:nw
   kn = mod(k,nw)+1;   % index of next corner
   kp = mod(k-2,nw)+1;   % index of prev corner
   knn= mod(k+1,nw)+1;   % index of next next corner
   if(~isinf(w(k))), ww = [ww; w(k)]; end
   
   if isnumeric(P{k})
      
      if(isinf(w(k)))
        tmp = w(knn)-w(kn);
        dw(k) = abs(tmp);
        pt{k} = @(t) w(kn) - (dw(k)./(t)-1)*tmp;
        
      elseif(isinf(w(kn)))
        tmp = w(k)-w(kp);
        dw(k) = abs(tmp);
        pt{k} = @(t) w(k) + (dw(k)./(dw(k)-t)-1)*tmp;
        
      elseif length(P{k}) == 1                      % straight arc
         dw(k) = abs(w(kn)-w(k));                   % distance to next corner
         pt{k} = @(t) w(k) + t*(w(kn)-w(k))/dw(k);  % parametrization of arc
      else                                          %     circular arc
         r = P{k}(2);                               % radius of arc
         a = w(k); b = w(kn); ab = abs(b-a);        % endpoints of arc
         theta = asin(ab/(2*r));                    % half-angle of arc
         c = a + r*exp(1i*(pi/2-theta))*(b-a)/ab;   % center of arc
         dw(k) = 2*theta*r;                         % arc length of arc
         pt{k} = @(t) c - ...
             r*exp(1i*(pi/2+t/r-theta))*(b-a)/ab;   % parametrization of arc
         ww = [ww; pt{k}(linspace(0,dw(k),50)')];
      end
      k = k+1;
   else
      error('LAPLACE:parseinputs','general boundary arcs not yet implemented')
   end
end
ww = [ww; ww(1)]; 
Zplot = ww;

%% Next treat the boundary conditions
g = cell(nw,1);
[g{:}] = deal(@(z) nan(size(z))); %default
if(demoflag<=1)
    g{1} = @(z) -1+0*z;
elseif(demoflag==2)
    g{1} = @(z) 0*z; g{2} = @(z) 1-(2*imag(z)-1).^2; g{end} = @(z) (1-imag(z).^2)/2;
elseif(demoflag==3)
    [g{[1,2,3,end]}] = deal(@(z) nan(size(z))+1i/1.5); 
end

j = 1;
while j < nargin
   j = j+1;
   v = varargin{j-1};

   if ~ischar(v)                 % This block specifies Dirichlet bndry data g.
      if isa(v,'cell')           % if cell array, nothing to change
         g = v;
      elseif isa(v,'double')     % if vector, convert to cell array of fun. handles
         for k = 1:nw
            g{k} = @(z) v(k) + 0*z;
         end
      elseif isa(v,'function_handle')  % if fun. handle, convert to cell array
         for k = 1:nw
            g{k} = @(z) v(z);
         end
      else
         error('LAPLACE:parseinputs','boundary data g not in correct form')
      end
   elseif strcmp(v,'tol'), j = j+1; tol = varargin{j-1};
   elseif strcmp(v,'steps'), steps = 1; plots = 1;
   elseif strcmp(v,'noplots'), plots = 0;
   elseif strcmp(v,'slow'), slow = 1;
   elseif strcmp(v,'rel'), rel = 1;
   elseif strcmp(v,'noarnoldi'), arnoldi = 0;
   elseif strcmp(v,'aaa')
      if exist('aaa') == 2, aaaflag = 1;
      else error('LAPLACE:parseinputs','Chebfun aaa is not in the path'), end
   else error('LAPLACE:parseinputs','Unrecognized string input')
   end
end

[ubkg, fbkg, g]=parse_chan(w, g);

continuous = 1;         % check for disc. bndry data if 'rel' not specified
continuous = ~any(isinf(w)); 
for k = 1:nw
   j = mod(k-2,nw)+1;
   gkk = g{k}(w(k)); gjk = g{j}(w(k));
   if abs(gkk-gjk) > tol || isnan(gkk) || isnan(gjk)
      continuous = 0;   % continuity not enforced at Neumann corners
   end
end
if ~continuous
   rel = 1;
end

end   % end of parseinputs


function [ubkg, fbkg, g] = parse_chan(w,g)
    nw = length(w);
    winf = reshape(find(isinf(w)),1,[]);
    ni = length(winf);

    if(ni==0)
        fbkg = @(z) zeros(size(z));
        ubkg = @(z) zeros(size(z));
    else
        fbkg = @(z) 0*z;
        ubkg = @(z) 0*z;
        H = [1 -1 1 -1; 1 1 1 1; 0 1 -2 3; 0 1 2 3];
        for k = winf
            kpp= mod(k-3,nw)+1;
            kp = mod(k-2,nw)+1;
            kn = mod(k,nw)+1;
            w1 = w(kp); w2 = w(kn); w0 = (w1+w2)/2;
            g1 = g{kpp}; g2 = g{kn};
            sk = w1-w(kpp);
            w2 = w2 - real((w2-w1)*conj(sk))/conj(sk); 
            a0 = (w2-w1)/2i;
      
            s = a0/sk;
            coef = H\[imag(g1(0)); imag(g2(0)); 0; 0];
            fbkg = @(z) fbkg(z) + imag(goursat_cubic(coef,s,(z-w0)/a0));
            ubkg = @(z) ubkg(z) + (1/a0)*velocity_cubic(coef,s,(z-w0)/a0);
        end
    end
    ubkg = @(z) conj(ubkg(z));
    if(ni>0)
        for k = 1:length(g)
            if(isnan(g{k}(0)))
                g{k} = @(z) g{k}(z) - 1i*fbkg(z);
                g{k} = @(z) -ubkg(z);
            else
                g{k} = @(z) g{k}(z) - ubkg(z);
            end
        end
    end
end

function [f]=goursat_cubic(c, s, z)
    z2 = imag(z).*z;
    z3 = 0.25*(3*conj(z)-z).*z.^2;
    f = 1i*c(1)+z.*c(2)+z2.*c(3)+z3.*c(4);
    f = f.*(1+tanh(s*z))/2;
end

function [u]=velocity_cubic(c, s, z)
    z2 = imag(z).*z;
    z3 = 0.25*(3*conj(z)-z).*z.^2;
    f = 1i*c(1)+z.*c(2)+z2.*c(3)+z3.*c(4);
	th = (1+tanh(s*z))/2;
    u1 = c(2) + c(3)*(2*z-conj(z))/2i + c(4)*(0.75*z.*(2*conj(z)-z));
    u2 =        c(3)*conj(z)/2i + c(4)*(0.75*conj(z).^2);
    u = (0.5*s*sech(s*z).^2).*f + th.*u1 - conj(th).*u2;
end