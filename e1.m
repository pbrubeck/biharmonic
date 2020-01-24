function e1(kmax)
npoly=4; 
w=exp((2i*pi/npoly)*(1:npoly)')/(exp(2i*pi/npoly)-1); 
w0=0;
w=w-w0; 
w=w./max(max(imag(w)),max(real(w)));


res=zeros(kmax,1);
dofs=zeros(kmax,1);
for k=1:kmax
    N=k^2;
    n=ceil(N/2);
    [pol,zs,un]=newpoles(N,w);
    % Boundary data
    tol=1E-4;
    u0=(abs(un-1)<=tol);
    u1=zeros(size(zs));
    [u,res(k),dofs(k)]=bih(n,zs,un,u0,u1,pol);
end

% Plotting
nplt=64;
xx=linspace(min(real(w)),max(real(w)),nplt);
yy=linspace(min(imag(w)),max(imag(w)),nplt);
[xx,yy]=ndgrid(xx,yy); zz=xx+1i*yy;
[uu,du]=u(zz);
du([1,end],[1,end])=nan;
inpolygonc=@(z,w)inpolygon(real(z),imag(z),real(w),imag(w)); 
%uu(~inpolygonc(zz,w)) = NaN;

figure(1);
surf(real(zz),imag(zz),real(uu));
dx=max(real(zz(:)))-min(real(zz(:))); dy=max(imag(zz(:)))-min(imag(zz(:)));
dz=max(uu(:))-min(uu(:)); daspect([1,1,2*dz/min(dx,dy)]);
xlabel('x'); ylabel('y'); view(3); colorbar();
%print('-dpng','e1_soln');

figure(2); subplot(1,2,1);
surf(real(zz),imag(zz),real(du)); 
dz=max(real(du(:)))-min(real(du(:))); daspect([1,1,2*dz/min(dx,dy)]);
xlabel('x'); ylabel('y'); view(3); title('u_x');

figure(2); subplot(1,2,2);
surf(real(zz),imag(zz),imag(du)); 
dz=max(imag(du(:)))-min(imag(du(:))); daspect([1,1,2*dz/min(dx,dy)]);
xlabel('x'); ylabel('y'); view(3); title('u_y');
%print('-dpng','e1_grad');

figure(3); subplot(1,2,1);
semilogy(sqrt(dofs),res,'-ok');
grid on; set(gca,'xminorgrid','off','yminorgrid','off');
xlabel('sqrt(DoFs)'); ylabel('Weighted residual');

figure(3); subplot(1,2,2);
plot(w([1:end,1]),'-k'); hold on;
plot(zs,'.k'); plot(pol,'.r'); hold off; axis equal;
print('-depsc','e1_conv');
end