function T = ttimes_dense(A, T, ind, divide)
% Multiply the tensor T in the mode ind by the matrix A
	if size(T, ind) ~= size(A, 2)
		error('TTIMES::incompatible dimensions')
    end
    sz=size(T);
    sz(ind)=size(A, 1);
    if exist('divide', 'var') && divide
        T = fold(A \ unfold(T, ind), sz, ind);    
    else
        T = fold(A * unfold(T, ind), sz, ind);
    end
end
