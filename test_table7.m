if isempty(gcp('nocreate'))
    parpool('threads');
end

d = 3;
% Total timings:             [LYAPTIME, LYAPRES, ADITIME, ADIRES] 
sizes = 2.^[8:10];
data = zeros(length(sizes), 8);
hssoption('block-size', 128);

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
	for bsj = 1:3
		bs = [1024, 2^(6+bsj)];

		% Divide-and-conquer solvers
		tdac = tic;
		[XDAC, timedata] = dac_lyapnd(A, -C, 'nmin', bs, 'sA', sA, 'low_rank_solver', 'adi', 'tol', 1e-6, 'timedata', zeros(4, 1));
		tdac = toc(tdac);
		resdac = 0;
		[~, resdac] = compute_residual(A, XDAC, C);
		data(nj, 2*bsj+1) = tdac;
		data(nj, 2*bsj+2) = resdac;
		fprintf('    DAC, time = %f, res = %e, met = %s \n', tdac, resdac, 'adi');

		dlmwrite('test4_3D.dat', [ sizes(:), data ], '\t');
	end
end
