d = 2;
% Total timings:             [LYAPTIME, LYAPRES, ADITIME512, ADIRES512, ADITIME1024, ADIRES1024, ADITIME2048, ADIRES2048] 
sizes = 2.^[11:15];
data = nan(length(sizes), 8);
hssoption('block-size', 256);
bs = [2048, 4096, 8192];
for nj = 1 : length(sizes)
    nn = sizes(nj);
    n = ones(1, d) * nn;
    
    % Create the matrices defining the equation
    [A, sA, ~, C] = create_example_fractional_space(n);
    FA = cell(1, d);
    for i = 1 : d; FA{i} = full(A{i}); FA{i} = (FA{i} + FA{i}')/2; end

	% Dense solver
    tdense = tic;
    XD = lyapnd(FA, -C);
    tdense = toc(tdense);
    [~, resdense] = compute_residual(A, XD, C);
    fprintf('N = %d\n', nn);
    fprintf('  DENSE, time = %f, res = %e\n', tdense, resdense);    
    data(nj, 1) = tdense;
    data(nj, 2) = resdense;

	% Divide-and-conquer solvers
    for bsj = 1:length(bs)        
        tdac = tic;
        [XDAC, timedata] = dac_lyapnd(A, -C, 'nmin', bs(bsj), 'sA', sA, 'low_rank_solver', 'adi', 'tol', 1e-8, 'timedata', zeros(4, 1), 'spd_split', true);
        tdac = toc(tdac);
        [~, resdac] = compute_residual(A, XDAC, C);
        if bs(bsj) < size(A{1}, 1)
            data(nj, bsj*2+1) = tdac;
            data(nj, bsj*2+2) = resdac;
        end
        fprintf('    DAC, time = %f, res = %e, bs = %d \n', tdac, resdac, bs(bsj));       
    end

    dlmwrite('test_table4.dat', [ sizes(:), data ], '\t');
end
