function [X, timedata] = dac_lyapnd(A, C, varargin)
% Possible argin:
%
% tol 				stopping criterion for Krylov methods
% sA				sparse representations of the matrix coefficients
% spectra			estimates of the spectral intervals obtained by calling spectrum_tree
% debug 			enables sme debugging prints
% low_rank_solver	The low-rank solver for the update equations. One in [ 'adi', 'rk', 'ek' ].
% spd_split 		exploit the spd structure to halve the rank of the right hand side when computing updates
% timedata  		vector for time measurings
%
%-----------------------------------------------------------------------------------------------------

d = length(A);

p = inputParser;
p.addOptional('tol', 1e-8, @(x) isscalar(x) && (x>=0));
p.addOptional('sA', A, @iscell);
p.addOptional('spectra', {}, @iscell);
p.addOptional('debug', false, @islogical);
p.addOptional('low_rank_solver', 'adi');
p.addOptional('spd_split', false, @islogical);
p.addOptional('npoles', 0, @isscalar);
p.addOptional('nadil', 4, @isscalar);
p.addOptional('nmin', hssoption('block-size'), @isvector);
p.addOptional('timedata', [0; 0; 0; 0], @isvector); % time dense, time low-rank, time forming rhs', time spectra estimate
p.parse(varargin{:});


sA = p.Results.sA;
spectra = p.Results.spectra;
timedata = p.Results.timedata;

opts.debug = p.Results.debug;
opts.low_rank_solver = p.Results.low_rank_solver;
opts.spd_split = p.Results.spd_split;
opts.npoles = p.Results.npoles;
opts.nadil = p.Results.nadil;
opts.tol = p.Results.tol;
opts.nmin = p.Results.nmin;

if length(opts.nmin) == 1
    opts.nmin = opts.nmin * ones(1, d-1);
end

if length(opts.nmin) < d - 1
    error('nmin should be either scalar or a vector with d - 1 entries');
end
   
if opts.npoles == 0
    if strcmp(opts.low_rank_solver, 'adi')
        sadi = @(tol, a1, b1, a2, b2) ceil( ...
            log(4 / tol) * ...
            log(16 * ((a1+b2) * (b1+a2)) / ((a1+a2)*(b2+b1)) ) / pi^2 );
    else
        sadi = @(tol, a1, b1, a2, b2) opts.nadil * floor( ...
            log(tol) / ...
            ( log(4) - opts.nadil*pi^2 / log(16 * ((a1+b2) * (b1+a2)) / ((a1+a2)*(b2+b1)) ) ));
    end
else
    sadi = @(tol, a1, b1, a2, b2) opts.npoles;
end

opts.sadi = sadi;

if d > 2
    t1 = clock;
    X0 = solve_baserec(A, C, opts.nmin(d-1));
    t2 = clock;
    timedata(1) = etime(t2, t1);
else
    X0 = [];
end

[X, timedata] = dac_lyapnd_ric(A, C, {ones(1,d), size(C)}, sA, spectra, opts, timedata, X0);

end

function [X, timedata] = dac_lyapnd_ric(A, C, idx, sA, spectra, opts, timedata, X0)

d = length(A);
n = arrayfun(@(i) size(A{i}, 1), 1 : d);

% Decide on which dimensions we need to split
split_modes = [];
nmax = max(n);

for j = 1 : d
    if n(j) > nmax / 2 && ~is_leafnode(A{j}) && n(j) > opts.nmin(d-1)
        split_modes = [ split_modes, j ];
    end
end

found = ~isempty(split_modes);

if ~found    
	t1 = clock;
    AA = cell(1, d);
    for j = 1 : d
        AA{j} = full(sA{j});
        if opts.spd_split || issymmetric(sA{j})
            AA{j} = (AA{j} + AA{j}') / 2;
        end
    end

    ind = cell(1, d);
    for j = 1 : d
        ind{j} = idx{1}(j) : idx{2}(j);
    end
    
    if d > 2
        X = X0(ind{:});
        return;
    end
    
    if issymmetric(AA{1})
        X = lyapnd(AA, C(ind{:}));
    else

        if d == 2 
            X = lyap(AA{1}, AA{2}', C(ind{:}));
        else
            X = laplace_recursive(AA, -C(ind{:}), 16);
        end
    end
    t2 = clock;
    timedata(1) = timedata(1) + etime(t2, t1);
    % X = lyap(full(A{1}), full(A{2})', C);
    
    % In case we need to check the residual
    %R = -C;
    %for ii = 1 : d
    %    R = R - ttimes_dense(AA{ii}, X, ii);
    %end
    %norm(R(:)) / norm(C(:))
    
else
    if isempty(spectra)
    	t1 = clock;
        spectra = cell(1, length(A));
        for i = 1 : length(A)
            spectra{i} = spectrum_tree(A{i}, sA{i}, opts.low_rank_solver);
        end
        t2 = clock;
        timedata(4) = timedata(4) + etime(t2, t1);
    end
    
    if d == 2 && ( ~is_leafnode(A{1}) && ~is_leafnode(A{2}) ) && ( max(n) < 2 * min(n) )
        m1 = size(A{1}.A11, 1); n1 = size(A{2}.A11, 1);
	    if opts.spd_split
		    Adiag = A;
		    sAdiag = sA;
		    [U1, V1] = offdiag(A{1}, 'lower');
		    [U2, V2] = offdiag(A{2}, 'lower');
		    %Adiag{1}.A11 =  Adiag{1}.A11 + hss('low-rank', V1, V1, 'cluster', cluster(Adiag{1}.A11));
		    %Adiag{1}.A22 =  Adiag{1}.A22 + hss('low-rank', U1, U1);	
		    U1 = [V1; -U1];
		    V1 = -U1;
		    U2 = [V2; -U2];	
		    V2 = -U2;
		    Adiag{1} = symmetrize(Adiag{1} + hss('low-rank', U1, U1, 'cluster', cluster(Adiag{1})));
		    Adiag{2} = symmetrize(Adiag{2} + hss('low-rank', U2, U2, 'cluster', cluster(Adiag{2})));
		    if issparse(sAdiag{1})
			    sAdiag{1} = sAdiag{1} + sparse(U1) * sparse(U1)';
		    else
			    sAdiag{1} = Adiag{1};
		    end
		    if issparse(sAdiag{2})
			    sAdiag{2} = sAdiag{2} + sparse(U2) * sparse(U2)';
		    else
			    sAdiag{2} = Adiag{2};
		    end
	    else
		    Adiag = A; 
		    sAdiag = sA;
		    [U1, V1] = offdiag(A{1}, 'all');
		    [V2, U2] = offdiag(A{2}, 'all'); % We extract the offdiagonal blocks for dA' instead of dA
	    end	
         
        sAA = sAdiag;
        
        % Block X11
        AA = { Adiag{1}.A11, Adiag{2}.A11 };
        AA{1}.topnode = 1; AA{2}.topnode = 1;    
        if isa(sAA{j}, 'hss')
            sAA{1} = AA{1};
            sAA{2} = AA{2};
        else
            sAA{1} = sAdiag{1}(1 : m1, 1 : m1);
            sAA{2} = sAdiag{2}(1 : n1, 1 : n1);
        end
        
        [X11, timedata] = dac_lyapnd_ric(AA, C, { idx{1}, idx{1} + [m1, n1] - 1}, sAA, { spectra{1}.A11, spectra{2}.A11 }, opts, timedata, X0);
        
        % Block X12
        AA = { Adiag{1}.A11, Adiag{2}.A22 };
        AA{1}.topnode = 1; AA{2}.topnode = 1;    
        if isa(sAA{j}, 'hss')
            sAA{1} = AA{1};
            sAA{2} = AA{2};
        else
            sAA{1} = sAdiag{1}(1 : m1, 1 : m1);
            sAA{2} = sAdiag{2}(n1+1:end, n1+1:end);
        end
        
        [X12, timedata] = dac_lyapnd_ric(AA, C, { idx{1} + [0, n1], [ idx{1}(1) + m1-1, idx{2}(2) ]}, sAA, { spectra{1}.A11, spectra{2}.A22 }, opts, timedata, X0);
        
        % Block X21
        AA = { Adiag{1}.A22, Adiag{2}.A11 };
        AA{1}.topnode = 1; AA{2}.topnode = 1;    
        if isa(sAA{j}, 'hss')
            sAA{1} = AA{1};
            sAA{2} = AA{2};
        else
            sAA{1} = sAdiag{1}(m1+1:end, m1+1:end);
            sAA{2} = sAdiag{2}(1 : n1, 1 : n1);
        end
        
        [X21, timedata] = dac_lyapnd_ric(AA, C, { [idx{1}(1) + m1, idx{1}(2) ], [idx{2}(1), idx{1}(2) + n1 - 1] }, sAA, { spectra{1}.A22, spectra{2}.A11 }, opts, timedata, X0);
        
        % Block X22
        AA = { Adiag{1}.A22, Adiag{2}.A22 };
        AA{1}.topnode = 1; AA{2}.topnode = 1;    
        if isa(sAA{j}, 'hss')
            sAA{1} = AA{1};
            sAA{2} = AA{2};
        else
            sAA{1} = sAdiag{1}(m1+1:end, m1+1:end);
            sAA{2} = sAdiag{2}(n1+1:end, n1+1:end);
        end
        [X22, timedata] = dac_lyapnd_ric(AA, C, { idx{1} + [m1, n1], idx{2} }, sAA, { spectra{1}.A22, spectra{2}.A22 }, opts, timedata, X0);
        
        % Build the rhs for the low-rank equation
		t1 = clock;
        X = [ X11, X12 ; X21 X22 ];
        [U, V] = compress_factors([ U1, X * U2 ], [ X' * V1, V2 ], opts.tol); 
        t2 = clock;
        timedata(3) = timedata(3) + etime(t2, t1);   
        nrmrhs = froblr(U, V);
        
        a1 = spectra{1}.sp(1); b1 = spectra{1}.sp(2);
        a2 = spectra{2}.sp(1); b2 = spectra{2}.sp(2);    
        
        t1 = clock;
        switch opts.low_rank_solver
            case 'adi'
                s = opts.sadi(opts.tol, a1, b1, a2, b2);
                [p, q] = zolotarev_poles(s, a1, b1, a2, b2);
                [Xu, Xv] = fadi(sA{1}, sA{2}, U, V, p, q);
            case 'adi_l'
                s = opts.sadi(opts.tol, a1, b1, a2, b2);
                [p, q] = zolotarev_poles(opts.nadil, a1, b1, a2, b2);
                if isa(sA{1}, 'hss')
                	A1S = create_struct_poles(sA{1}, q);
                elseif issparse(sA{1})
                	A1S = create_struct_poles_sparse(sA{1}, q);
                else
                	error('Not supported format for sA{1}')
                end	
                if isa(sA{1}, 'hss')
                	A2S = create_struct_poles(sA{2}, -p);
                elseif issparse(sA{1})
                	A2S = create_struct_poles_sparse(sA{2}, -p);
                else
                	error('Not supported format for sA{2}')
                end	
                [Xu, Xv] = fadi(A1S, A2S, U, V, p, q, s);
            case 'rk'
                s = opts.sadi(opts.tol, a1, b1, a2, b2);
                if (s+2) * size(U, 2) > min(size(U, 1), size(V, 1))
                    s = floor(min(size(U, 1), size(V, 1)) / size(U, 2)) - 2;
                    warning('The estimated rational Krylov dimension is too large: please increase the block size');
                end            
                [p, q] = zolotarev_poles(s, a1, b1, a2, b2);
                [Xu, Xv] = rklyapsolve(sA{1}, sA{2}', U, V, p, q);
            case 'ek'
                [Xu, Xv] = ek_sylv(sA{1}, sA{2}', U, V, inf, @(r,n) r < opts.tol * nrmrhs, false, 'fro', true);
        end
        t2 = clock;
        timedata(2) = timedata(2) + etime(t2, t1);
        
 		t1 = clock;
        X = X + Xu * Xv';
        t2 = clock;
        timedata(3) = timedata(3) + etime(t2, t1);
    else
        nsplit = length(split_modes);
        U = cell(1, nsplit);
        V = cell(1, nsplit);
        
        % Make a copy of the coefficients that we can use in the recursive
        % calls
        Adiag = A;
        sAdiag = sA;
        
        X = zeros(n);
        
        newspectra11 = spectra;
        newspectra22 = spectra;
        
        for t = 1 : nsplit
            j = split_modes(t);
        
            % Split over dimension j
            if opts.spd_split
                [U{t}, V{t}] = offdiag(A{j}, 'lower');
                U{t} = [V{t}; -U{t}];
                V{t} = -U{t};
    
                Adiag{j} = symmetrize(Adiag{j} + hss('low-rank', U{t}, U{t}, 'cluster', cluster(Adiag{j})));
    
                if issparse(sAdiag{j})
                    sAdiag{j} = sAdiag{j} + sparse(U{t}) * sparse(U{t})';
                else
                    sAdiag{j} = Adiag{j};
                end
            else
                [U{t}, V{t}] = offdiag(A{j}, 'all');
            end
            
            % Compute the spectra of the diagonal blocks
            newspectra11{j} = spectra{j}.A11;
            newspectra22{j} = spectra{j}.A22;
        end
        
        % Solve the 2^nsplit equations by recursion
        for tt = 0 : 2^nsplit-1
            child_selector = dec2bin(tt, nsplit);
            
            AA = Adiag;
            sAA = sAdiag;
            
            newspectra = spectra;
            
            t1 = clock;
            
            % Indices of the part of the solution that we will compute; these
            % will be restricted in the next for loop
            indices = cell(1, d);
            for j = 1 : d
                indices{j} = 1 : n(j);
            end
            
            idx_rec = idx;

            for t = 1 : nsplit
                j = split_modes(t);
                
                if child_selector(t) == '0'
                    % Select A11 in mode j
                    AA{j} = Adiag{j}.A11;
                    
                    if isa(sAA{j}, 'hss')
                        sAA{j} = AA{j};
                        sAA{j}.topnode = 1;
                    else
                        sAA{j} = sAdiag{j}(1 : size(AA{j}, 1), 1 : size(AA{j}, 2));
                    end
                    idx_rec{2}(j) = idx_rec{1}(j) + size(AA{j}, 2) - 1;
                    indices{j} = 1 : size(AA{j}, 2);
                    newspectra{j} = newspectra11{j};
                else
                    % Select A22 in mode j
                    AA{j} = Adiag{j}.A22;
                    if isa(sAA{j}, 'hss')
                        sAA{j} = AA{j};
                        sAA{j}.topnode = 1;
                    else
                        sAA{j} = sAdiag{j}(size(AA{j}, 1) + 1 : end, size(AA{j}, 2) + 1 : end);
                    end
                    idx_rec{1}(j) = idx_rec{1}(j) + size(AA{j}, 2);
                    indices{j} = size(AA{j}, 2)+1:n(j);
                    newspectra{j} = newspectra22{j};
                end
                
                AA{j}.topnode = 1;
                % newspectra{j} = [ eigs(sAA{j}, 1, 'SM', 'Tolerance', 1e-2), eigs(sAA{j}, 1, 'LM', 'Tolerance', 1e-2) ];
            end
            t2 = clock;
            timedata(3) = timedata(3) + etime(t2, t1);
            
            % Solve the equation, overwriting the indices in the solution.            
            [X(indices{:}), timedata] = dac_lyapnd_ric(AA, C, idx_rec, sAA, newspectra, opts, timedata, X0);
        end
        
        % We now solve the update equations, in all modes that we have split;
        % we create a separate dX variable, since the X solving the equation at
        % the lower level is used to compute the RHS.

        Xu = cell(1, nsplit);
        Xv = Xu;
        
  		t1 = clock;      
        for t = 1:nsplit
        	j = split_modes(t);
            

             % Form the RHS using the computed solution X of the block-diagonal
             % system

             if exist('tensorprod', 'builtin')
                 V{t} = reshape(tensorprod(V{t}', X, 2, j), size(V{t}, 2), prod(size(X, [1:j-1,j+1:d])));
             else
                 V{t} = V{t}' * unfold(X, j);
             end
             V{t} = V{t}';
        end
        t2 = clock;
        timedata(3) = timedata(3) + etime(t2, t1);

        t1 = clock;
        
        parfor t = 1 : nsplit
             
             j = split_modes(t);
             [U{t}, V{t}] = compress_factors(U{t}, V{t}, opts.tol);
             
             AA = A;
             sAA = sA;
             AA(j) = [];
             sAA(j) = [];
             
    
    		 t1 = clock;
             if length(AA) == 1
                 B = sAA{1};
    
                 a2 = spectra{1}.sp(1);
                 b2 = spectra{1}.sp(2);
             else
                 B = struct(... 
                     'solve', @(nu, mu, x) op_solve(AA, nu, mu, x, opts, spectra([1:j-1,j+1:d]), sAA), ...
                     'multiply', @(rho, eta, x) rho * op_multiply(sAA, x) - eta * x, ...
                     'real', true, ...
                     'issymmetric', true);
    
                 a2 = 0; b2 = 0;
                 for jj = 1 : d
                     if jj ~= j
                         a2 = a2 + spectra{jj}.sp(1);
                         b2 = b2 + spectra{jj}.sp(2);
                     end
                 end    
             end
             % t2 = clock;
			 % timedata(2) = timedata(2) + etime(t2, t1);    
    
             a1 = spectra{j}.sp(1);
             b1 = spectra{j}.sp(2);    
    
             nrmrhs = froblr(U{t}, V{t});
    		 t1 = clock;
             switch opts.low_rank_solver
                 case 'adi'
                     s = opts.sadi(opts.tol, a1, b1, a2, b2);
                     [p, q] = zolotarev_poles(s, a1, b1, a2, b2);
                     [Xu{t}, Xv{t}] = fadi(sA{j}, B, U{t}, V{t}, p, q);
                 case 'rk'
                     % [Xu, Xv] = rk_sylv(@(index,i) poles_handle(a,b,c,d,i,index), sA{j}, B, U, V, inf, @(r, n) r < nrmrhs * opts.tol, opts.debug);
                     s = opts.sadi(opts.tol, a1, b1, a2, b2);
                     if (s+2) * size(U{t}, 2) > min(size(U{t}, 1), size(V{t}, 1))
                         s = floor(min(size(U{t}, 1), size(V{t}, 1)) / size(U{t}, 2)) - 2;
                         warning('The estimated rational Krylov dimension is too large: please increase the block size');
                     end
                     [p, q] = zolotarev_poles(s, a1, b1, a2, b2);
                     [Xu{t}, Xv{t}] = rklyapsolve(sA{j}, B, U{t}, V{t}, p, q);            
                 case 'ek'
                     [Xu{t}, Xv{t}] = ek_sylv(sA{j}, B, U{t}, V{t}, inf, @(r, n) r < nrmrhs * opts.tol, opts.debug);
             end
             % t2 = clock;
             % timedata(2) = timedata(2) + etime(t2, t1);             
        end

        t2 = clock;
        timedata(2) = timedata(2) + etime(t2, t1);
        
        for t = 1:nsplit
             t1 = clock;
             j = split_modes(t);
             X = X + fold(Xu{t} * Xv{t}', n, j);
             t2 = clock;
             timedata(3) = timedata(3) + etime(t2, t1);
        end
    end
end

end
%----------------------------------------------------
function nrm = normoperator(A)
    nrm = 0;
    for k = 1 : length(A)
        nrm = nrm + normest(A{k}, 1e-2);
    end
end

function w = op_multiply(A, v)
    d = length(A);     
    n = arrayfun(@(i) size(A{i}, 1), 1 : d);
    w = zeros(size(v));
    
    for i = 1 : size(v, 2)
        vi = v(:,i);
        vi = reshape(vi, n);
        wi = zeros(size(vi));
        
        for j = 1 : d
            wi = wi + ttimes_dense(A{j}, vi, j);
        end

        w(:,i) = wi(:);
    end
end
%----------------------------------------------------
function w = op_solve(A, nu, mu, v, opts, spectra, sA)
% w = (nu*A - mu*I) \ v

d = length(A);
n = arrayfun(@(i) size(A{i}, 1), 1 : d);

if nu == 0
    w = (-1/mu) * v;
else    
    s = sign(nu);
    
    nu = nu / s;
    mu = mu / s;
    
    for j = 1 : d
        % A{j} = nu * A{j} - mu/d * eye(size(A{j}), 'like', A{j});
        
        if isa(sA{j}, 'hss')
            sA{j} = shift(nu * sA{j}, -mu/d);
        else
            sA{j} = nu * sA{j} - mu/d * eye(size(sA{j}), 'like', sA{j});
        end
        
        A{j} = shift(nu * A{j}, -mu/d); 
        spectra{j}.sp = nu * spectra{j}.sp - mu/d;
    end
    
    w = zeros(size(v));

    if ( (d == 2) && (max(n) <= 512)) &&  issymmetric(A{1}) && issymmetric(A{2})
        w = lyapnd({ full(A{1}), full(A{2}) }, -reshape(v, [n, size(v, 2)]));
        w = reshape(w, size(v));
    else
        for i = 1 : size(v, 2)
            vi = reshape(v(:, i), n);    
            % wi = dac_lyapnd(A, -vi, 'spd_split', spd_split, 'low_rank_solver', low_rank_solver, 'tol', tol, 'spectra', spectra, 'sA', sA);
            if d > 2
                error('ahi ahi ahi');
            end
            wi = dac_lyapnd_ric(A, -vi, {ones(1, length(n)), size(vi)}, sA, spectra, opts, zeros(4,1), []);
            w(:,i) = wi(:);
        end
    end
    
    w = w * s;
    
end

end
%----------------------------------------------------
function p = poles_handle(a, b, c, d, j, index)
    [p, q] = zolotarev_eds(a, b, c, d, j);
    
    if index == 2
        p = -q;
    end
end
%-----------------------------------------------------
function X = solve_baserec(A, C, nmin)
    n = size(C);
    d = length(A);

    eq = find_equations(A, C, nmin);

    parallel = true;
    if length(eq) == 1
        parallel = false;
    end
    
    if parallel
       parfor k = 1 : length(eq)
         XX{k} = lyapnd(eq{k}.A, eq{k}.C);
       end
    else
       for k = 1 : length(eq)
         XX{k} = lyapnd(eq{k}.A, eq{k}.C);
       end
    end

    
    X = zeros(n);
    for k = 1 : length(eq)
        indices = cell(1, d);
        for j = 1 : d
            indices{j} = eq{k}.start(j) : eq{k}.end(j);
        end
        X(indices{:}) = XX{k};
    end
end

function eq = find_equations(A, C, nmin)
	d = length(A);
	sz = zeros(1, d);
	for j = 1 : d
		g{j} = get_cluster(A{j}, nmin);
		sz(j) = length(g{j});
	end

    I = cell(1, d);    
    eq = cell(1, prod(sz));
	for k = 1 : prod(sz)
        [I{:}] = ind2sub(sz, k);
        eq{k} = struct;
        
        for j = 1 : d
        	eq{k}.A{j} = g{j}(I{j}).A;
        	eq{k}.start(j) = g{j}(I{j}).start;
        	eq{k}.end(j) = g{j}(I{j}).end;
        end
        
        Cindices = cell(1, d);
        for j = 1 : d
            Cindices{j} = eq{k}.start(j) : eq{k}.end(j);
        end
        eq{k}.C = C(Cindices{:});
	end
end

function cl = get_cluster(A, nmin)
	if is_leafnode(A) || size(A, 1) <= nmin
		% cl = [1; size(A, 1)];
		cl = struct('A', full(A), 'start', 1, 'end', size(A, 1));
	else
		cl1 = get_cluster(A.A11, nmin);
		cl2 = get_cluster(A.A22, nmin);
        for j = 1 : length(cl2)
            cl2(j).start = cl2(j).start + size(A.A11, 1);
            cl2(j).end = cl2(j).end + size(A.A11, 1);
        end
		cl = [cl1, cl2];
	end
end

