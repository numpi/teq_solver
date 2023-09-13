d = 2;
% Total timings:             [LYAPTIME, LYAPRES, ADITIME, ADIRES, RKTIME, RKRES, EKTIME, EKRES] 
sizes = 2.^[12:15];
data = zeros(length(sizes), 8);
hssoption('block-size', 256);
metods = {'adi', 'rk', 'ek'};
nadil = 5;
bs = 2048;  
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
    for met = 1:length(metods) 
        %profile on        
        tdac = tic;
        [XDAC, timedata] = dac_lyapnd(A, -C, 'nmin', bs, 'sA', sA, 'low_rank_solver', metods{met}, 'nadil', nadil, 'tol', 1e-8 * (1e4 *(met==2)+1), 'timedata', zeros(4, 1), 'spd_split', true);
        tdac = toc(tdac);
        % profile viewer
        [~, resdac] = compute_residual(A, XDAC, C);      
        data(nj, met*2+1) = tdac;
        data(nj, met*2+2) = resdac;
        fprintf('    DAC, time = %f, res = %e, met = %s \n', tdac, resdac, metods{met});  
        if strcmp(metods{met}, 'adi')
            timedata = timedata / tdac * 100;
            dlmwrite(sprintf('test1_2D_bis_met=%s_hist_n=%d.dat', metods{met}, n(1)), timedata, '\t');
        end
    end

    dlmwrite('test_table2.dat', [ sizes(:), data ], '\t');
end

data = generate_latex_histogram2('test1_2D_bis_met=adi_hist', [4096, 8192, 16384]);


