k = 1;
n = [ 800 800 ];

hssoption('block-size', 100);

eigmin = 1e-4;
lmin = eigmin;
shift = eigmin;
ntest = 10;


data = zeros(0, 3);

dd = linspace(0, 1, n(1)).^10;
A = cell(1, 2);

for ii = 1 : 32
    for jj = 1 : ntest
        AA = diag(-rand(n(1) - 1, 1), -1); AA = AA + AA';
        v = AA * ones(n(1), 1);
        fAA = AA - diag(v) + eigmin * eye(n(1));
        AA = hss(fAA);
        A0 = { AA, AA };
        
        [Q, l] = eig(fAA, 'vector');
        [~, p] = sort(l, 'descend');
        ip = p; ip(p)=[1:n];
        
        C = Q' * diag(dd(ip) .* rand(1, n(1))) * Q;
        for j = 1 : 2
            A{j} = A0{j} - (lmin - shift) * hss('eye', n(1));
        end

        cnd = cond(full(A{1}));
        sadi = @(tol, a1, b1, a2, b2) ceil( ...
            log(4 / tol) * ...
            log(16 * ((a1+b2) * (b1+a2)) / ((a1+a2)*(b2+b1)) ) / pi^2 );
        s = sadi(1e-6, 1, cnd, 1, cnd);    
		tempo = tic;
        [X] = dac_lyapnd(A, C, 'tol', 1e-6, 'nmin', 32, 'npoles', s);
        tempo = toc(tempo);
        [R, res] = compute_residual(A, X, -C);

        data = [ data ; cnd, res,  cnd * 1e-5 ];

        fprintf('jj = %d, cond = %e, res = %e, time = %.2f, npoles = %d\n', jj, data(end,1), data(end, 2), tempo, s);
        dlmwrite('test_figure2_right.dat', data, '\t');
    end

    shift = shift / sqrt(2);
end

close all

loglog(data(:, 1), sqrt(data(:, 1)) * 1e1 * data(1, 2) / sqrt(data(1, 1)), 'b--');
hold on;
scatter(data(:, 1), data(:, 2), 'rx');




