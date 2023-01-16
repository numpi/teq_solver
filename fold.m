function T = fold(A, sz, ind)
% Tensorize the matrix A as a tensor of size sz, placing the index ind as first index (reverse operation w.r.t to unfold)
	d = length(sz);
	p = [ind, 1:ind-1, ind+1:d];
	ip(p) = [1:d];
	T = permute(reshape(A, sz(p)), ip);
end
