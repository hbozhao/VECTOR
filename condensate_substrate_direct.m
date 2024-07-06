function params = condensate_substrate_direct(params)
  %Due to copyright issues and hence inability to share my modified version of ode15s that uses Krylov method. This version uses the direct method, where the jacobian is an estimate of the full jacobian, where the off block diagonal betweeen the condensate concentration field and the substrate degrees of freedom is set to 0.
  %set up params to be used in solver_FDCH for simulating condensate-substrate interaction
  %in params, provide:
  %params.substrate.location [N,ndims], each row is the coordinates of the substrate center. Mandatory always. Typically it is the initial condition of the coordinates of the center are additional degrees of freedom, this function will set params.add.N and params.substrate.N.
  %params.substrate.U0, interaction strength, a vector whose length is the same as the number of species. Positive for repulsive. Negative for attractive.   %all substrates are identical and have the same interaction kernel with all species, use U0 to specify interaction strength with each species
  %params.substrate.kernel.func, a function handle which is a function of r2 and params.substrate.kernel.params (NB. it it is r squared!) it should return G(r2), dG/dr2, and d^2G/d(r2)^2, if requested
  %params.substrate.kernel.cutoff, set of cutoff length, will only compute kernel in a box in x_i-cutoff<= x<=x_i+cutoff, where x_i is the center substrate location.
  %note that we also have a predefined kernel which can be accessed by setting kernel.func to the name.
  %Gaussian: exp(-r2/(2*sigma^2)), user must provide sigma in params.substrate.kernel.params. Will automatically turn it to exp(k*r2), will automatically set params.substrate.kernel.cutoff = 3*sigma (user can also customize after using this function)
  %params.substrate.dynamics:
  %dashpot: dx/dt = mobility * F, mobility = 1/friction (alternative defintion, if friction is provided, it will also be converted to mobility, overwriting older mobility). Provide mobility or friction as a scalar or vector (each substrate can can different friction coefficient) in params.substrate
  %Kevin-Voigt, dx/dt = mobility * (F - k (x - x0) ), mobility same usage as above. User need to provide k and x0 in params.substrate.spring.stiffness and params.substrate.spring.origin, which must have the dimension of [N,ndims] (same as params.substrate.location). Stiffness, similar to mobility, can either be a scalar or an array with different stiffness. Anisotropic spring not yet supported.
  %SFM: standard fluid model (Jeffreys). dx/dt = mobility * (F - k(x - x0)), and tau * dx0/dt = (x-x0), mobility, k, and x0 same usage as above. provide tau in params.substrate.tau. tau can be a scalar or an array just like friction and stiffnesss
  %For all models, condensate-substrate interaction is periodic.
  %Rouse: Rouse model. See details below. User can tune params.substrate.Rouse.kappa_to_eta. Currently only support identical mobility for all beads on all chains.
  %params.grid. coordinates of the grid, expressed as a cell, {x,y,...}
  %params.mu provide the homogeneous mu struct

  d = length(params.N);
  if isfield(params,'Nspecies') && params.Nspecies>1
    params.substrate.U0 = permute(params.substrate.U0(:),[1+(1:d),1]);
  else
    params.Nspecies = 1;
  end
  %number of substrate pointss
  Ns = size(params.substrate.location,1);
  params.substrate.N = Ns;
  %original degrees of freedom (not including additional varibles)
  params.dof = prod([params.N,params.Nspecies]);
  if isequal(params.substrate.kernel.func,'Gaussian')
    params.substrate.kernel.func = @Gaussian;
    sigma = params.substrate.kernel.params;
    params.substrate.kernel.params = -1/(2*params.substrate.kernel.params^2);
    if ~isfield(params.substrate.kernel,'cutoff')
      params.substrate.kernel.cutoff = 3*sigma;
    end
  end
  params.add.RHS = @RHS;
  params.add.jacobian = @jacobian;
  Nadd = d*Ns;
  if isequal(params.substrate.dynamics,'SFM')
    Nadd = Nadd*2;
  end
  %further specify whether to compute off-diagonal terms of dFdr (dFx/dry, and dFy/drx)
  if ~isfield(params.substrate,'dFdr_diag')
    params.substrate.dFdr_diag = true;
  end
  if isfield(params.substrate,'friction')
    params.substrate.mobility = 1./params.substrate.friction;
  end
  params.substrate.mobility = params.substrate.mobility(:);
  if isfield(params.substrate,'spring') && isfield(params.substrate.spring,'stiffness')
    params.substrate.spring.stiffness = params.substrate.spring.stiffness(:);
    stiff = params.substrate.spring.stiffness;
    if isscalar(stiff)
      stiff = stiff*ones(Ns,1);
    end
    if params.substrate.dFdr_diag
      springdFdr = -stiff;
    else
      springdFdr = zeros([d,d,Ns]);
      for i = 1:Ns
        for j = 1:d
          springdFdr(j,j,i) = -stiff(i);
        end
      end
    end
    params.substrate.spring.dFdr = springdFdr;
  end
  if isfield(params.substrate,'tau')
    params.substrate.tau = params.substrate.tau(:);
  end
  if isequal(params.substrate.dynamics,'Rouse')
    %with Rouse model, model half of the chain, with the middle one being the bead that interacts with the condensate and have symmetric boundary condition on the chain. The last bead has Dirichlet boundary condition
    params.substrate.Rouse.N = 100; %the number of beads in each chain in addition to the middle bead
    Nb = params.substrate.Rouse.N+1; %number of beads per chain
    params.substrate.Rouse.Ns = Ns; %number of chains
    if ~isfield(params.substrate.Rouse,'kappa_to_eta')
      kappa_to_eta = 1; %kappa is the coefficient in front of the nabla, eta is the coefficient on the LHS in front of dx/dt (see droppullet_theory). We can always set this to 1.
    else
      kappa_to_eta = params.substrate.Rouse.kappa_to_eta;
    end
    %corresponding to kappa_to_eta = 1 (which is when time is normalized to process time), let's set dx = 0.1, so the total length is 10 on either side. The factor in front of second order differencing is hence kappa / dx^2. Note that eta is 1/mobility
    Rouse_dx = 0.1;
    params.substrate.Rouse.factor = kappa_to_eta / params.substrate.mobility / Rouse_dx^2;
    Nadd = d * Ns * Nb;
    %keep the initial condition, need it for Dirichlet boundary condition
    params.substrate.Rouse.init = params.substrate.location;
    %expand location to all beads. user can use it for initial condition.
    %the ordering is that beads on the same chain goes first, then the next chain
    params.substrate.location = repelem(params.substrate.location,Nb,1);
    %now compute the jacobian for the Rouse dynamics
    %the degrees of freedom are ordered as: beads on a chain, chain, dimensionality. Get the operator for the beads on the same chain, then do kron.
    jac_Rouse = spdiags(ones(Nb,1)*[1,-2,1],-1:1,Nb,Nb);
    jac_Rouse(1,1) = -1;
    jac_Rouse = kron(speye(Ns),jac_Rouse);
    jac_Rouse = kron(speye(d),jac_Rouse);
    params.substrate.Rouse.jac = jac_Rouse * params.substrate.Rouse.factor * params.substrate.mobility;
  end
  params.add.N = Nadd;

end

function [dx,mu] = RHS(t,y,yadd,params)
  d = length(params.N);
  if isequal(params.substrate.dynamics,'SFM')
    nx = params.add.N/2;
    x0 = yadd((nx+1):2*nx);
    x0 = reshape(x0,[],d);
    yadd = yadd(1:nx);
  end
  x = reshape(yadd,[],d);
  if isequal(params.substrate.dynamics,'Rouse')
    %use the first bead to represent the bead in the middle of the Rouse chain, which is the only bead needed for the capillary interaction
    xx = reshape(x,params.substrate.Rouse.N+1,[],d);
    x = squeeze(xx(1,:,:));
  end
  Nx = size(x,1);
  dV = prod(params.dx);
  mu = params.mu.func(y,params.mu.params);
  ysum = sum(y .* params.substrate.U0,d+1);
  F = zeros(size(x));
  U = zeros(params.N);
  for i = 1:Nx
    roi = true(params.N);
    for j = 1:d
      dr = params.grid{j} - x(i,j);
      dr = dr - params.L(j).*round(dr./params.L(j));
      roi = roi & abs(dr)<params.substrate.kernel.cutoff;
    end
    r2 = 0;
    for j = 1:d
      r2 = r2 + (params.grid{j}(roi)-x(i,j)).^2;
    end
    [G,dG] = params.substrate.kernel.func(r2,params.substrate.kernel.params);
    U(roi) = U(roi) + G;
    for j = 1:d
      F(i,j) = sum(2*(params.grid{j}(roi) - x(i,j)) .* dG .* ysum(roi));
    end
  end
  %turn F from summation to integration
  F = F * dV;
  mu = mu + U .* params.substrate.U0;
  switch params.substrate.dynamics
  case 'Kevin-Voigt'
    F = F - params.substrate.spring.stiffness.*(x - params.substrate.spring.origin);
  case 'SFM'
    F = F - params.substrate.spring.stiffness.*(x - x0);
  case 'Rouse'
    %F is the force on the middle bead
    F0 = F;
    F = zeros(size(xx));
    %the first bead has symmetric boundary condition
    F(1,:,:) = xx(2,:,:) - xx(1,:,:);
    %the last bead has Dirichlet boundary condition
    F(end,:,:) = xx(end-1,:,:) + permute(params.substrate.Rouse.init,[3,1,2]) - 2*xx(end,:,:);
    F(2:end-1,:,:) = xx(1:end-2,:,:) + xx(3:end,:,:) - 2*xx(2:end-1,:,:);
    F = F * params.substrate.Rouse.factor;
    F(1,:,:) = F(1,:,:) + permute(F0,[3,1,2]);
  end
  dx = F .* params.substrate.mobility;
  dx = dx(:);
  if isequal(params.substrate.dynamics,'SFM')
    dx0 = (x-x0) ./ params.substrate.tau;
    dx = [dx; dx0(:)];
  end
end

function jac = jacobian(t,y,yadd,params)
  d = length(params.N);
  if isequal(params.substrate.dynamics,'SFM')
    nx = params.add.N/2;
    x0 = yadd((nx+1):2*nx);
    x0 = reshape(x0,[],d);
    yadd = yadd(1:nx);
  end
  x = reshape(yadd,[],d);
  if isequal(params.substrate.dynamics,'Rouse')
    %use the first bead to represent the bead in the middle of the Rouse chain, which is the only bead needed for the capillary interaction
    xx = reshape(x,params.substrate.Rouse.N+1,[],d);
    x = squeeze(xx(1,:,:));
    Nb = params.substrate.Rouse.N+1; %number of beads per chain
  end
  Nx = size(x,1);
  dV = prod(params.dx);
  dx = params.dx;
  ysum = sum(y .* params.substrate.U0,d+1);
  dFdr = zeros([d,d,Nx]);
  for i = 1:Nx
    roi = true(params.N);
    for j = 1:d
      dr = params.grid{j} - x(i,j);
      dr = dr - params.L(j).*round(dr./params.L(j));
      roi = roi & abs(dr)<params.substrate.kernel.cutoff;
    end
    r2 = 0;
    for j = 1:d
      r2 = r2 + (params.grid{j}(roi)-x(i,j)).^2;
    end
    [G,dG,d2G] = params.substrate.kernel.func(r2,params.substrate.kernel.params);
    dFdr_eye = -sum(2 * dG .* ysum(roi));
    for j = 1:d
      dFdr(j,j,i) = dFdr_eye;
      for k = 1:d
        dFdr(j,k,i) = dFdr(j,k,i) - 4*sum( (params.grid{j}(roi) - x(i,j)) .* (params.grid{k}(roi) - x(i,k)) .* d2G .* ysum(roi));
      end
    end
  end
  dFdr = dFdr * dV;
  switch params.substrate.dynamics
  case {'Kevin-Voigt','SFM'}
    dFdr = dFdr + params.substrate.spring.dFdr;
    tau = params.substrate.tau;
    if isscalar(tau)
      tau = tau*ones(Nx,1);
    end
  end
  dFdr = dFdr .* params.substrate.mobility;
  %now convert to a matrix
  jac_N = params.add.N;
  jac = sparse(jac_N,jac_N);
  for p = 1:d
    for q = 1:d
      if isequal(params.substrate.dynamics,'Rouse')
        jac = jac + sparse(((1:Nx)-1)*Nb+1+(p-1)*Nb*Nx, ((1:Nx)-1)*Nb+1+(q-1)*Nb*Nx, squeeze(dFdr(p,q,:)), jac_N, jac_N);
      else
        jac = jac + sparse((1:Nx)+(p-1)*Nx, (1:Nx)+(q-1)*Nx, squeeze(dFdr(p,q,:)), jac_N, jac_N);
        if p==q && ismember(params.substrate.dynamics,{'SFM','Kevin-Voigt'})
          jac = jac + sparse(d*Nx+(1:Nx)+(p-1)*Nx, d*Nx+(1:Nx)+(q-1)*Nx, -1./tau, jac_N, jac_N);
        end
      end
    end
  end
  if isequal(params.substrate.dynamics,'Rouse')
    %add the Rouse operator
    jac = jac + params.substrate.Rouse.jac;
  end
end


function [G,dG,d2G] = Gaussian(r2,k)
  G = exp(k*r2);
  if nargout>1
    dG = G .* k;
  end
  if nargout>2
    d2G = dG .* k;
  end
end