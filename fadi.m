function [W, Y, res] = fadi(A1, A2, U, V, p, q, s, U0, V0)
%FADI Solve the Sylvester equation A1*X + X*A2' + U*V' = 0. 
% 
% The solution is returned in factored form X = W*Y'.

res = [];

m = size(U, 1);
n = size(V, 1);

if ~exist('U0', 'var')
    W = zeros(m, 0);
    Y = zeros(n, 0);
else
    W = U0; Y = V0;
    error('Not yet implemented');
end

if ~exist('s', 'var')
    s = length(p);
end

if length(p) < s
    l = ceil(s / length(p));
    p = kron(p(:), ones(l, 1));
    p = p(1:s);
    q = kron(q(:), ones(l, 1));
    q = q(1:s);
end

if isa(A1, 'hss')
    A1 = struct(... 
         'solve', @(nu, mu, x) shift(nu * A1, -mu) \ x , ...
         'multiply', @(rho, eta, x) rho * (A1 * x) - eta * x, ...
         'real', true, ...
         'issymmetric', true);
end

if ~isstruct(A1)
    Im = eye(m, 'like', A1);
    W = (A1 - q(1)*Im) \ U;
else
    W = A1.solve(1, q(1), U);
end

if isa(A2, 'hss')
    A2 = struct(... 
         'solve', @(nu, mu, x) shift(nu * A2, -mu) \ x , ...
         'multiply', @(rho, eta, x) rho * (A2 * x) - eta * x, ...
         'real', true, ...
         'issymmetric', true);
end

if ~isstruct(A2)
    In = eye(n, 'like', A2);
    Y = (A2 + p(1)*In) \ V;
else
    Y = A2.solve(1, -p(1), V);
end

k = size(U, 2);

for i = 2 : s
    if ~isstruct(A1)        
        W = [ W, (A1 - q(i)*Im) \ ((A1 - p(i-1)*Im) * W(:, end-k+1:end)) ];
    else
        W = [ W, A1.solve(1, q(i), A1.multiply(1, p(i-1), W(:, end-k+1:end))) ];
    end
    
    if ~isstruct(A2)
        Y = [ Y, (A2 + p(i)*In) \ ((A2 + q(i-1)*In) * Y(:, end-k+1:end)) ];
    else
        Y = [ Y, A2.solve(1, -p(i), A2.multiply(1, -q(i-1), Y(:, end-k+1:end))) ];
    end
end

D = kron( sqrt(-diag(q-p)), eye(k) );
W = -W * D;
Y = Y * D;

if nargout >= 3
    res = zeros(1, s);
    X = lyap(A1, A2, U*V');
    
    for i = 1 : s
        Xu = W(:, 1:k*i);
        Xv = Y(:, 1:k*i);
        
        res(i) = norm((A1*Xu)*Xv' + Xu*(Xv'*A2) + U*V', 'fro') / norm(X, 'fro');
    end    
end

end

