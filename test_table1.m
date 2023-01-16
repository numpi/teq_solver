d = 2;
% Total timings:             [LYAPTIME, LYAPRES, ADITIME, ADIRES, RKTIME, RKRES, EKTIME, EKRES] 
sizes = 2.^[9:15];
data = zeros(length(sizes), 8);
hssoption('block-size', 256);
metods = {'adi', 'rk', 'ek'};
nadil = 5;
bs = 512;
for nj = 1 : length(sizes)
    nn = sizes(nj);
    n = ones(1, d) * nn;
    
    % Create the matrices defining the equation
    [A, sA, ~, C] = create_example_laplacian(n);
    FA = cell(1, d);
    for i = 1 : d; FA{i} = full(A{i}); end

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
    for met = 1:length(metods)         
        tdac = tic;
        %profile on
        [XDAC, timedata] = dac_lyapnd(A, -C, 'nmin', bs, 'sA', sA, 'low_rank_solver', metods{met}, 'nadil', nadil, 'tol', 1e-8 * (1e3 *(met==2)+1), 'timedata', zeros(4, 1), 'spd_split', false);
        %profile viewer
        tdac = toc(tdac);
        [~, resdac] = compute_residual(A, XDAC, C);      
        data(nj, met*2+1) = tdac;
        data(nj, met*2+2) = resdac;
        fprintf('    DAC, time = %f, res = %e, met = %s \n', tdac, resdac, metods{met});
        if strcmp(metods{met}, 'adi')
            timedata = timedata / tdac * 100;
            dlmwrite(sprintf('test1_2D_met=%s_hist_n=%d.dat', metods{met}, n(1)), timedata, '\t');
        end
    end

    dlmwrite('test1_2D.dat', [ sizes(:), data ], '\t');
    
end

data = generate_latex_histogram2('test1_2D_met=adi_hist', [4096, 8192, 16384]);
