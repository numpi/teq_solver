function [A, sA, X, C, FA] = create_example_fractional_space(n)
%FRACTIONAL SPACE DIFFERENTIAL PROBLEM, NONSYMMETRIC

d = length(n);

A = {};
sA = {};
FA = {};

alpha = 1.5;

for j = 1 : d
    [am, ap] = fractional_symbol(alpha, n(j));
    ap(n(j)) = 0;
    A{j} = hss('toeplitz', am, ap);
    A{j} = A{j} + A{j}';
    A{j} = symmetrize(A{j}, 'up');
    sA{j} = A{j};
    if nargout > 4
        FA{j} = toeplitz(am + ap);
    end
end

X = randn(n);
C = zeros(n);

for j = 1 : d
    C = C + ttimes_dense(A{j}, X, j);
end

end
