function [Xu, Xv] = rklyapsolve(A, B, U, V, p, q)
%
% Solve AX + XB' + U*V = 0, using the poles in p and q. 
%

bs = size(U, 2);

if isa(A, 'hss')
    A = struct(... 
         'solve', @(nu, mu, x) shift(nu * A, -mu) \ x , ...
         'multiply', @(rho, eta, x) rho * (A * x) - eta * x, ...
         'real', true, ...
         'issymmetric', true);
end

if isa(B, 'hss')
    B = struct(... 
         'solve', @(nu, mu, x) shift(nu * B, -mu) \ x , ...
         'multiply', @(rho, eta, x) rho * (B * x) - eta * x, ...
         'real', true, ...
         'issymmetric', true);
end

[VA, KA, HA, ~] = rat_krylov(A, U, [q,  inf]);
[VB, KB, HB, ~] = rat_krylov(B, V, [-p, inf]);

Cprojected = (VA(:,1:bs)' * U) * (VB(:,1:bs)'*V)';
Cprojected(bs * (1+length(p)), bs * (1+length(p))) = 0;

Ap = HA(1:end-bs, :) / KA(1:end-bs, :);
Bp = HB(1:end-bs, :) / KB(1:end-bs, :);

Y = lyap(Ap, Bp', Cprojected);

Xu = VA(:, 1:size(Y, 1)) * Y;
Xv = VB(:, 1:size(Y,1));

end

