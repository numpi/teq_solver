function T = spectrum_tree(A, sA, low_rank_method)
    T = struct();
    if strcmp(low_rank_method, 'ek')
        T.sp = [1 2]; % Not needed for ek
    else
    	if issparse(sA)
    		[L, U, tp, tq] = lu(sA, 'vector');
    		ip = tp; ip(tp) = 1:length(tp);
    		iq = tq; ip(tq) = 1:length(tq);
    		Ainv = @(x) sparse_lu_solve(L, U, ip, iq, x);

    	else
			if issymmetric(sA)
				if isa(sA, 'hss')
					F = chol(sA);
					Ainv = @(x) chol_solve(F, x);
				else
					Ainv = @(x) A\x;
				end	
			else
				if isa(sA, 'hss')
					F = ulv(sA);
					Ainv = @(x) ulv_solve(F, x);
				else
					Ainv = @(x) A\x;
				end	
			end
		end	
        %T.sp = [ eigs(sA, 1, 'SM', struct('tol', 1e-2)), eigs(sA, 1, 'LM', struct('tol', 1e-2)) ];
        T.sp = [ eigs(Ainv, size(sA, 1), 1, 'SM', struct('tol', 1e-2)), eigs(sA, 1, 'LM', struct('tol', 1e-2)) ];
        T.sp = real(T.sp); % Possibile inghippo
    end
    
    if ~is_leafnode(A)        
        if isa(sA, 'hss')
            sA11 = A.A11;
            sA22 = A.A22;
            sA11.topnode = 1;
            sA22.topnode = 1;
        else
            n1 = size(A.A11, 2);
            sA11 = sA(1 : n1, 1 : n1);
            sA22 = sA(n1+1:end, n1+1:end);
        end
        
        T.A11 = spectrum_tree(A.A11, sA11, low_rank_method);
        T.A22 = spectrum_tree(A.A22, sA22, low_rank_method);
    end
end

function y = sparse_lu_solve(L, U, p, q, x)
	y = U\(L\(x(q,:))); 
	y = y(p, :);
end
