function varargout = solver_FDCH_direct(tspan,y0,params)
  %This function uses finite difference for solving Cahn-Hilliard equation
  %Due to copyright issues and hence inability to share my modified version of ode15s that uses Krylov method. This version uses the direct method. Exact jacobian for  an approximate jacobian 
  %input:
  %tspan: time span for simulation
  %y0: initial condition. Must provide full initial condition if additional degrees of freedom is used. Otherwise, if empty, will use random field with an average of 0.5, if a vector of length Nspecies, will use use a random field with an average given by y0
  %params.N, L, dx: number of grid points, length, grid spacing in each direction, at least provide 2
  %params.Nspecies: number of species. must provide for multicomponent Cahn-Hilliard
  %params.D0: mobility tensor (by default 1)
  %params.kappa: the gradient penalty in free energy must be set
  %params.kappa can either be a scalar, in which case the kappa tensor be kappa*I
  %params.kappa can also be a matrix like D
  %params.mu: the homogeneous part of the chemical potential. This is a struct that functional handle and parameters (see below). If params.mu is nonexistent, uses regular solution model using params.omega (by default 3)
  %the user must use params.mu to return the chemical potential, params.mu.func should return an array of size [params.N,Nspecies], the last dimension is the i in mu_i, params.mu.grad should return an array of size [params.N,Nspecies,Nspecies], the last two dimensions are the i and j in dmu_i / dc_j
  %params.diagonalD: set to true if D is diagonal. If multicomponent and mobility tensor is a constant, params.D0 should be a vector (row or column), if not homgeneous, params.D0 or output of params.D should be an array of size [params.N,Nspecies]. dD/dc should be an array of [prams.N,Nspecies,Nspecies], where the last two dimension i and j and d(Dii)/dc_j
  %params.reaction: add reaction to all species except for nutrient
  %this mode is on if user provides params.R, a function handle that takes mu and c as input (params.R(mu,c)), returns the reaction rate for all species (just like mu) and dR/dmu, dR/dc if requested, if dR/dmu or dR/dc is empty, then they are zero. dR/dmu and dR/dc follows params.mu usage. If single component,
  %user can also provide a time-dependent reaction rate, so the input will be params.R(mu,c,t).
  %params.message: turn on message during solution
  %params.add: add additional degrees of freedom to the equations.
  %additional y is placed at the end
  %now we support when additional degrees of freedom changes the chemical potential. params.add.RHS (defined below) will take over params.mu in terms of its output (mu and dmu). dmu must still return dmu/dy.
  %params.add.N (mandatory) number of additional variables.
  %params.add.RHS (mandatory) is a function handle that takes in (t,y,yadd,params), where y is the input when there is no additional variables (same below), yadd is the additional variables in a single column, and params is the entire params struct, and returns {dyadd,mu} where dyadd is the RHS that corresponds to the additional degrees of freedom and must be a column vector. 
  %params.add.jacobian (mandatory)  the input is (t,y,yadd,params)
  %in add mode, user must provide y0. Default y0 will result in an error
  %output:
  %if the length of tspan is 1, then return {yp0,params} where yp0 is dy/dt at tspan
  %otherwise, return {tout,y,params} where tout is the time points and y is the solution at those time points

  if ~isfield(params,'dx') && isfield(params,'N') && isfield(params,'L')
    params.dx = params.L ./ params.N;
  elseif ~isfield(params,'N') && isfield(params,'dx') && isfield(params,'L')
    params.N = params.L ./ params.dx;
  elseif ~isfield(params,'L') && isfield(params,'N') && isfield(params,'dx')
    params.L = params.N .* params.dx;
  elseif ~all(isfield(params,{'dx','N','L'}))
    params.N = [128,128];
    params.L = [5,5];
    params.dx = params.L./params.N;
  end
  n = prod(params.N);
  params.n = n;
  N = params.N;
  d = length(N);
  params.d = d;
  dx = params.dx;
  if ~isfield(params,'D0') || isempty(params.D0)
    params.D0 = 1;
  end
  if ~isfield(params,'mu') || isempty(params.mu)
    if ~isfield(params,'omega') || isempty(params.omega)
      params.omega = 3;
    end
    params.mu.func = @(x,coeff) log(x./(1-x)) + params.omega*(1-2*x);
    params.mu.grad = @(x,coeff) 1./x + 1./(1-x) - 2*params.omega;
    params.mu.params = [];
  end
  params.multicomponent = isfield(params,'Nspecies') && params.Nspecies>1;
  if ~isfield(params,'diagonalD')
    params.diagonalD = false;
  end
  if params.multicomponent
    if params.diagonalD
      params.D0_jac = kron(diag(params.D0),sparse(1:n,1:n,ones(1,n)));
      params.D0 = permute(params.D0(:),[1+(1:d),1]);
    else
      %spatially constant diffusivity only
      %form a matrix form of D0 that includes all degrees of freedom (in the form of jacobian)
      params.D0_jac = kron(params.D0,sparse(1:n,1:n,ones(1,n)));
      %for D0 tensor product, let the species (as well as higher order tensor index of species such as diffusivity tensor) index be the last index, shuffle to the back
      species_perm_ind = [2+(1:d),2,1]; %2 then 1 because D*mu is sum_j(D_ij * mu_j) = sum_j(D'_ji * mu_j), although D should be symmetric
      params.D0 = permute(params.D0,species_perm_ind);
    end
    Nspecies = params.Nspecies;
    if ~isscalar(params.kappa)
      %similar to D0
      %form a matrix form of kappa that includes all degrees of freedom
      params.kappa_jac = kron(params.kappa(1:Nspecies,1:Nspecies),sparse(1:n,1:n,ones(1,n)));
      species_perm_ind = [2+(1:d),2,1];
      params.kappa = permute(params.kappa,species_perm_ind);
    else
      params.kappa_jac = params.kappa;
    end
  else
    Nspecies = 1;
    params.Nspecies = Nspecies;
    if ~isscalar(params.kappa)
      %not multicomponent, but sometimes with substrate, kappa can be multicomponent
      params.kappa_jac = params.kappa(1); %jacobian only needs the first component
      species_perm_ind = [2+(1:d),2,1];
      params.kappa = permute(params.kappa,species_perm_ind);
    else
      params.kappa_jac = params.kappa;
    end
  end
  dof = n*Nspecies;
  params.hasReaction = isfield(params,'R') && isa(params.R,'function_handle');
  if Nspecies==1 && d==1
    N = [N,1];
    params.N = N;
  end

  %create operators
  params.L2 = sparse(n,n);
  for i = 1:d
    Ni = N(i); %number of points in the ith direction
    %operator list in this dimension
    op = [];
    E = sparse(1:Ni-1,2:Ni,1,Ni,Ni);
    op.L2 = E + E' - 2*speye(Ni); %L2 in one dimension
    op.L2(1,Ni) = 1;
    op.L2(Ni,1) = 1;
    op.L2 = op.L2 / dx(i)^2;
    for j = 1:d
      if j==i
        O = op.L2;
      else
        O = speye(N(j));
      end
      if j==1
        Oi = O;
      else
        Oi = kron(O,Oi);
      end
    end
    params.L2 = params.L2 + Oi;
  end
  if params.multicomponent
    params.L2 = kron(sparse(1:Nspecies,1:Nspecies,ones(1,Nspecies)),params.L2);
  end

  if params.multicomponent
    %i and j for sparse to create a jacobian (including all dof) for dmu/dc from dmu returned by params.mu, which has the unit of [params.N,params.Nspecies,params.Nspecies]
    params.dmu2jac_ind.i = reshape(repmat(1:dof,1,Nspecies),[],1);
    params.dmu2jac_ind.j = reshape(repmat(permute(reshape(1:dof,n,Nspecies),[1,3,2]),[1,Nspecies,1]),[],1);
  end

  options = odeset;
  params.message = isfield(params,'message') && params.message;
  if params.message
    options = odeset(options,'OutputFcn',@(t,y,flag) simOutput(t,y,flag,true));
  end
  if isfield(params,'options') && ~isempty(params.options)
    options = odeset(options,params.options);
  end

  %set up intialization
  if isempty(y0) || length(y0)==Nspecies
    %initialization
    if isempty(y0)
      n0 = 0.5;
    else
      n0 = y0;
    end
    if params.multicomponent
      n0 = n0(:);
      n0 = permute(n0,[1+(1:d),1]);
    end
    sigma = 0.01;
    rng(1);
    y0 = n0 + sigma*randn([N,Nspecies]);
    y0 = y0(:);
  end
  if isempty(tspan)
    tspan = linspace(0,0.5,100);
  end
  t0 = tspan(1);
  tf = tspan(end);
  params.t0 = t0;
  params.tf = tf;
  yp0 = RHS(t0,y0,params);
  if length(tspan)==1
    varargout = {yp0,params};
    return;
  end
  options = odeset(options,'InitialSlope',yp0);

  odeFcn = @(t,y) RHS(t,y,params);
  options = odeset(options,'Jacobian',@(t,y) jacobian(t,y,params));
  [tout,y] = ode15s(odeFcn,tspan,y0,options);
  varargout = {tout,y,params};

end

function dy = RHS(t,y,params)
  N = params.N;
  n = params.n;
  dx = params.dx;
  d = params.d;
  if isfield(params,'add')
    yadd = y((end-params.add.N+1):end);
    y = y(1:(end-params.add.N));
  end
  if params.multicomponent
    N = [N, params.Nspecies];
    n = n*params.Nspecies;
  end
  y = reshape(y,N);
  if isfield(params,'add')
    [dyadd,mu] = params.add.RHS(t,y,yadd,params);
  else
    mu = customizeFunGrad(params,'mu','fun',y);
  end
  d2y = zeros(size(y));
  for i = 1:d
    d2y = d2y + (circshift(y,1,i)+circshift(y,-1,i)-2*y)/dx(i)^2;
  end
  if isscalar(params.kappa)
    mu = mu - params.kappa * d2y;
  else
    mu = mu - squeeze(sum(params.kappa.*d2y,d+1));
  end
  D = params.D0;
  dy = zeros(N);
  %diffusion
  for i = 1:d
    dy = dy + (circshift(mu,1,i)+circshift(mu,-1,i)-2*mu)/dx(i)^2;
  end
  if ~params.multicomponent || params.diagonalD
    dy = D.*dy;
  else
    dy = squeeze(sum(D.*dy,d+1));
  end

  if params.hasReaction
    if nargin(params.R)==2
      dy = dy + params.R(mu,y);
    elseif nargin(params.R)==3
      dy = dy + params.R(mu,y,t);
    end
  end

  dy = dy(:);

  if isfield(params,'add')
    dy = [dy; dyadd];
  end

end

function dfdy = jacobian(t,y,params)
  %L2 = Lxx + Lyy + ... (Laplacian)
  %L is a cell, L{i} is Li (central differencing in the ith direction), needed if mobility is not constant
  N = params.N;
  n = params.n;
  d = params.d;
  dx = params.dx;
  L2 = params.L2;
  if isfield(params,'add')
    yadd = y((end-params.add.N+1):end);
    y = y(1:(end-params.add.N));
  end
  if params.multicomponent
    N = [N, params.Nspecies];
    n = n*params.Nspecies;
  end
  if params.multicomponent
    y = reshape(y,N);
  end

  [mu,dmu] = customizeFunGrad(params,'mu','fungrad',y);
  if params.multicomponent
    % dmu = params.mu.grad_multi(y,params.mu.params);
    dmusz = size(dmu);
    %size of dmu excluding spatial, just number of components in mu and c
    dmusz = dmusz(d+1:end);
    if dmusz(1)~=dmusz(2)
      %asymmetry dmu
      dmu2jac_ind.i = reshape(repmat(1:params.n*dmusz(1),1,dmusz(2)),[],1);
      dmu2jac_ind.j = reshape(repmat(permute(reshape(1:params.n*dmusz(2),params.n,dmusz(2)),[1,3,2]),[1,dmusz(1),1]),[],1);
      dmu = sparse(dmu2jac_ind.i,dmu2jac_ind.j,dmu(:));
    else
      dmu = sparse(params.dmu2jac_ind.i,params.dmu2jac_ind.j,dmu(:),n,n);
    end
    dmu = dmu - params.kappa_jac * L2;
  else
    dmu = sparse(1:n,1:n,dmu(:)) - params.kappa_jac * L2;
  end

  if params.multicomponent
    dfdy = params.D0_jac * L2 * dmu;
  else
    dfdy = params.D0 * L2 * dmu;
  end

  if params.hasReaction
    if nargin(params.R)==2
      [~, dRdmu, dRdc] = params.R(mu,y);
    else
      [~, dRdmu, dRdc] = params.R(mu,y,t);
    end
    if ~isempty(dRdmu)
      if params.multicomponent
        dRdmu = sparse(params.dmu2jac_ind.i,params.dmu2jac_ind.j,dRdmu(:),n,n);
      else
        dRdmu = sparse(1:n,1:n,dRdmu);
      end
      dfdy = dfdy + dRdmu * dmu;
    end
    if ~isempty(dRdc)
      if params.multicomponent
        dRdc = sparse(params.dmu2jac_ind.i,params.dmu2jac_ind.j,dRdc(:),n,n);
      else
        dRdc = sparse(1:n,1:n,dRdc);
      end
      dfdy = dfdy + dRdc;
    end
  end

  %as an approximation of the actual jacobian, here we include only the block diagonal terms that correspond to y (the concentration field) and yadd (additional degrees of freedom), while ignoring the off diagonal terms
  if isfield(params,'add') && isfield(params.add,'jacobian')
    jac_add = params.add.jacobian(t,y,yadd,params);
    dfdy = blkdiag(dfdy, jac_add);
  end

end