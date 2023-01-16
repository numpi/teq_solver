function [XX, time] = lyap_superdc(A, C, nmin, tol)
%

d = length(A);

if d ~= 2
    error('ahi ahi ahi');
end

n = length(A{1});

[tr, m] = npart(n, nmin);

[DA,UA,RA,BA] = mat2hss(full(A{1}),tr,m);
[DB,UB,RB,BB] = mat2hss(full(A{2}),tr,m);

t1 = clock;

[dA,QA] = superdc(DA,UA,BA,RA,tr,tol);
[dB,QB] = superdc(DB,UB,BB,RB,tr,tol);

% Transformed RHS
CC = superdcmv(QA, C,  1);
CC = superdcmv(QB, CC', 1)';

CC = CC ./ (dA + dB');

XX = superdcmv(QA, CC, 0);
XX = superdcmv(QB, XX', 0)';

time = etime(clock, t1);

end