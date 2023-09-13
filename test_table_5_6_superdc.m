d = 2;
% Total timings:             [LYAPTIME, LYAPRES, ADITIME, ADIRES, RKTIME, RKRES, EKTIME, EKRES]

addpath HSS/hss/
addpath FMM1D/fmm1d/src/
addpath FMM1D/fmm1d/utils/
addpath SuperDC/superdc_1.0.0/src/
addpath SuperDC/superdc_1.0.0/utils

sizes = 2.^[12:15];
data = zeros(length(sizes), 4);

nmin = 2048;
tol = 1e-10;

for nj = 1 : length(sizes)
    nn = sizes(nj);
    n = ones(1, d) * nn;
    
    % Create the matrices defining the equation
    [A, sA, ~, C] = create_example_laplacian(n);
    FA = cell(1, d);
    for i = 1 : d; FA{i} = full(A{i}); end

    [X, tlap] = lyap_superdc(FA, C, nmin, tol);
    [R, reslap] = compute_residual(FA, X, C);

    [A, sA, ~, C, FA] = create_example_fractional_space(n);

    [X, tlapfrac] = lyap_superdc(FA, C, nmin, tol);
    [R, reslapfrac] = compute_residual(FA, X, C);

    data(nj, :) = [ tlap, reslap, tlapfrac, reslapfrac ];

    fprintf('n = %d, tlap = %f, reslap = %e, tlapfrac = %f, reslapfrac = %e\n', ...
        nn, tlap, reslap, tlapfrac, reslapfrac);

    dlmwrite('test_table_5_6_superdc.dat', [ sizes(:), data ], '\t');
    
end
