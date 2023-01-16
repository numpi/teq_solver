function [A, sA, X, C] = create_example_laplacian(n)

d = length(n);

A = cell(1, d);
sA = cell(1, d);

for j = 1 : d
    A{j} = hssgallery('laplacian', n(j));
    sA{j} = sparse(A{j});
end

C = randn(prod(n), 1);
C = reshape(C, n);
X = 0;
return;

C(n) = 0;

for j = 1 : d
    C = C + ttimes_dense(A{j}, X, j);
end

end
