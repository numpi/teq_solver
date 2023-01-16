function C = lyapnd(A, C)
%LYAPND Solve the equation \sum X \times_i A{i} + C = 0 by diagonalization.
%
% If C is a tensor with d + 1 indices, one more than the number of 
% coefficients A{i}, then the equation is solved for all slices of C. 

d = length(A);
% n = arrayfun(@(i) size(A{i}, 1), 1 : d);

% Diagonalize the coefficients
% t = tic;
Q = cell(1, d);
D = cell(1, d+1);
for k = 1 : d
    [Q{k}, D{k}] = eig(A{k}, 'vector');
    if issymmetric(A{k})
        C = ttimes_dense(Q{k}', C, k);
    else
        C = ttimes_dense(Q{k}, C, k, true);
    end
end
% fprintf("ttimes_dense_1: %f secs\n", toc(t));

% t = tic;
D{d+1} = zeros(size(C, d+1), 1);

M = D{1};
for k = 2 : d + 1
    M = M + D{k}.';
    M=M(:);
end

C = -C ./ reshape(M, size(C));
% fprintf("cauchy: %f secs\n", toc(t));

% t = tic;
for k = 1 : d
    C = ttimes_dense(Q{k}, C, k);
end
% fprintf("ttimes_dense_2: %f secs\n", toc(t));

end
