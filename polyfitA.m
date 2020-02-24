function [c,H] = polyfitA(x,y,n)   % polyfit with Arnoldi
M = length(x); V = ones(M,n+1); H = zeros(n+1,n); 
for k = 1:n       
    v = x(:).*V(:,k);
    for j = 1:k
        H(j,k) = V(:,j)'*v;
        v = v - H(j,k)*V(:,j);
    end
    
%     for j = 1:k     % Gram-Schmidt twice
%         DelH = V(:,j)'*v/M;
%         v = v - DelH*V(:,j);
%         H(j,k) = H(j,k)+DelH;
%     end
    H(k+1,k) = norm(v);
    V(:,k+1) = v/H(k+1,k);
end
c = V\y;
end
