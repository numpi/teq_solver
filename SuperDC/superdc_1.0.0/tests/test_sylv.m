n = 8192;
w = 1;
nmin = 1024;
[tr, m] = npart(n, nmin); % HSS partition info

tol = 1e-10;

A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n);
B = A;

X = randn(n, n);
C = A*X + X*B;

[DA,UA,RA,BA] = band2hss(full(A),w,tr,m);
[DB,UB,RB,BB] = band2hss(full(B),w,tr,m);

tic;
[dA,QA] = superdc(DA,UA,BA,RA,tr,tol);
[dB,QB] = superdc(DB,UB,BB,RB,tr,tol);

% Transformed RHS
CC = superdcmv(QA, C,  1);
CC = superdcmv(QB, CC', 1)';

CC = CC ./ (dA + dB');

XX = superdcmv(QA, CC, 0);
XX = superdcmv(QB, XX', 0)';
toc

norm(X - XX, 'fro') / norm(XX, 'fro')
