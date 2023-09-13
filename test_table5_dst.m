d = 2;
% Total timings:             [LYAPTIME, LYAPRES, ADITIME, ADIRES, RKTIME, RKRES, EKTIME, EKRES]
sizes = 2.^[11:15];
data = zeros(length(sizes), 2);

nmin = 2048;
tol = 1e-10;

for nj = 1 : length(sizes)
    nn = sizes(nj);
    n = ones(1, d) * nn;   
    % Create the matrices defining the equation
    [A, sA, ~, C] = create_example_laplacian(n);
    FA = cell(1, d);
    for i = 1 : d; FA{i} = full(A{i}); end

    tlap = tic;
    l = 2 + 2 * cos(pi ./ (n(1)+1) * (n(1):-1:1));
    X = C;
    X = dst(dst(X')');
    X = X ./ ( l + l' );
    X = idst(idst(X')');
    tlap = toc(tlap);
    
    [R, reslap] = compute_residual(FA, X, C);

	data(nj, :) = [ tlap, reslap ];

    fprintf('n = %d, tlap = %f, reslap = %e\n', ...
        nn, tlap, reslap);

    dlmwrite('test_table5_dst.dat', [ sizes(:), data ], '\t');
    
end
