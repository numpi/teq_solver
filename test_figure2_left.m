k = 8;
n = [ 256 256 ];
Q = orth(triu(randn(n), -k), 0);
hssoption('block-size', 32);

data = zeros(0, 3);
ntest = 100;
dd = linspace(0, 1, n(1)).^10;

for i = 1 : 24
    p = 1 + (i - 1) * 0.05;
    d = 2 + 2 * cos(pi * (1 : n(1)) ./ (n(1) + 1));

    cc = max(d.^p) / min(d.^p);
    for j = 1 : ntest

        Q = orth(triu(randn(n), -k));
        A = hss(Q * diag(d.^p) * Q');

        C = zeros(n(1));
        for i = 1 : (256 / 32)
            for j = 1 : (256 / 32)
                AL = full(A((i-1) * 32 + 1 : i * 32, (i-1) * 32 + 1 : i * 32));
                AR = full(A((j-1)*32 + 1 : j*32, (j-1)*32 + 1 : j*32));
                [VL, DL] = eig(AL);
                [VR, DR] = eig(AR);
                [~, lpos] = min(diag(DL)); 
                [~, rpos] = min(diag(DR));
                C((i-1) * 32 + 1 : i * 32, (j-1)*32 + 1 : j*32) = VL(:, lpos) * VR(:, rpos)';
            end
        end
        tempo = tic;
        [X] = dac_lyapnd({ A, A }, C, 'tol', 1e-6, 'nmin', 32);
        tempo = toc(tempo);
        [R, res] = compute_residual({ A, A }, X, -C);
        data = [ data ; cc , res , cc * 1e-11 ];
        fprintf('j = %d, cond = %e, res = %e, time = %.2f\n', j, data(end,1), data(end, 2), tempo);
    end

    dlmwrite('test_figure2_left.dat', data, '\t');
end

figure

loglog(data(:, 1), data(:, 1) * 1e1 * data(1, 2) / data(1, 1), 'b--');
hold on;

loglog(data(:, 1), sqrt(data(:, 1)) * 1e1 * data(1, 2) / sqrt(data(1, 1)), 'b--');
scatter(data(:, 1), data(:, 2), 'rx');



