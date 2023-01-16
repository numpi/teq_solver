function T = unfold(T, ind)
% Matricized the tensor T, positioning the index ind as row index
	n = size(T, ind);
    sz = size(T);
	d = ndims(T);
	p = [ind, 1:ind-1, ind+1:d];
	T = permute(T, p);
    
    if n ~= 0
        T = reshape(T, [n, numel(T)/n]);
    else
        T = zeros(0, prod(sz([ 1 : ind-1,ind+1:d ])));
    end
end
