function nrm = froblr(U, V)
%FROBLR 

[~, RU] = qr(U, 0);
[~, RV] = qr(V, 0);

nrm = norm(RU * RV', 'fro');

end

