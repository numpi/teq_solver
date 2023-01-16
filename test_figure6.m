d = 3;
% Total timings:             [LYAPTIME, LYAPRES, ADITIME, ADIRES] 
sizes = 2.^[6,6:10];
data = nan(length(sizes), 2);
hssoption('block-size', 32);
nadil = 5;
bs = [32, 32];

for nj = 1 : length(sizes)
    nn = sizes(nj);
    n = nn * ones(1,3);

    % Create the matrices defining the equation
    [A, sA, ~, C] = create_example_laplacian(n);
    fprintf('N = %d\n', nn);

    % Divide-and-conquer solvers
    tdac = tic;
    [XDAC, timedata] = dac_lyapnd(A, -C, 'nmin', bs, 'sA', sA, 'low_rank_solver', 'adi', 'nadil', nadil, 'tol', 1e-6, 'timedata', zeros(4, 1));
    tdac = toc(tdac);
    [~, resdac] = compute_residual(A, XDAC, C);      
    data(nj, 1) = tdac;
    data(nj, 2) = resdac;
    fprintf('    DAC, time = %f, res = %e, met = %s \n', tdac, resdac, 'adi');
    dlmwrite('test6_3D.dat', [ sizes(:), data ], '\t');
    timedata = timedata / tdac * 100;
    dlmwrite(sprintf('test6_3D_hist_n=%d.dat', n(1)), timedata, '\t');
end

generate_latex_histogram2('test6_3D_hist', [256, 512, 1024]);
