function [R, nrm] = compute_residual(A, X, C)
%
% Compute the residual of the tensor equation
%
%  R := X x_1 A{1} + ... + X x_d A{d} - C


R = C;
for j = 1 : length(A)
    R = R - ttimes_dense(A{j}, X, j);
end

nrm = norm(R(:)) / norm(C(:));


end