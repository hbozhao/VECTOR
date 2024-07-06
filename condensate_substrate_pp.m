function [F,mu] = condensate_substrate_pp(params,y,force,varargin)
  %post-processing of condensate substrate interaction
  %compute the force on the substrate F, or if requested the chemical potential mu
  %provide: params, y (solution in the format of solver_FDCH output, Ntime*dof), and force
  %where force can be 'net', which is the net force on the substrate except (drag force is the opposite of net force), including capillary (chemical) and viscoelastic (mechanical) forces, depending params.substrate.dynamics
  %or it can be 'capillary', which is only the capillary (chemical) force due to condensate_substrate interaction
  %params must be returned by condensate_substrate
  %the output F has the dimension of Ntime * Npoints * d, where Npoints is the number of substrate center, and d is the dimensionality
  %the optional output mu has the dimension of Ntime * params.N * params.Nspecies
  ps = inputParser;
  addParameter(ps,'tspan',[]);
  ps.CaseSensitive = false;
  parse(ps,varargin{:});
  ps = ps.Results;


  dynamics = params.substrate.dynamics;
  Nt = size(y,1);
  if isequal(dynamics,'fixed')
    params.substrate.dynamics = 'dashpot';
    params.mu = params.muh;
    params = condensate_substrate(params);
    x0 = params.substrate.location(:).';
    y = [y, repmat(x0,Nt,1)];
  end
  if ~isfield(params,'dx')
    params.dx = params.L ./ params.N;
  end
  d = length(params.N);
  Ns = params.substrate.N;
  if isequal(force,'capillary')
    params.substrate.dynamics = 'dashpot';
    if isequal(dynamics,'Rouse')
      Ns = params.substrate.Rouse.Ns;
      params.substrate.N = Ns;
      real_add_N = params.add.N;
      params.add.N = Ns*d;
    end
    Nadd = Ns*d;
  else
    Nadd = params.add.N;
  end
  F = zeros(Nt,Ns,d);
  if isequal(dynamics,'Rouse') && isequal(force,'net')
    Nchain = params.substrate.Rouse.Ns;
    F = zeros(Nt,Nchain,d);
  end
  if nargout>1
    mu = zeros(Nt,prod(params.N)*params.Nspecies);
  end

  if ~isscalar(params.kappa)
    %not multicomponent, but sometimes with substrate, kappa can be multicomponent
    species_perm_ind = [2+(1:d),2,1];
    params.kappa = permute(params.kappa,species_perm_ind);
  end

  for i = 1:Nt
    thisyadd = y(i,(end-params.add.N)+(1:Nadd));
    thisyadd = thisyadd(:);
    thisy = y(i,1:(end-params.add.N));
    if isequal(dynamics,'Rouse') && isequal(force,'capillary')
      %just get the capillary force on the middle bead
      thisyadd = y(i,(end-real_add_N)+(1:real_add_N));
      thisyadd = reshape(thisyadd,[],Ns,d);
      thisyadd = squeeze(thisyadd(1,:,:));
      thisyadd = thisyadd(:);
      thisy = y(i,1:(end-real_add_N));
    end
    thisy = reshape(thisy,[params.N,params.Nspecies]);
    if isempty(ps.tspan)
      thist = [];
    else
      thist = ps.tspan(i);
    end
    [dx,thismu] = params.add.RHS(thist,thisy,thisyadd,params);
    dx = dx(1:(Ns*d));
    dx = reshape(dx,Ns,d);
    dx = dx ./ params.substrate.mobility;
    if isequal(dynamics,'Rouse') && isequal(force,'net')
      %just get the force on the middle bead
      dx = reshape(dx,[],Nchain,d);
      dx = squeeze(dx(1,:,:));
    end
    F(i,:,:) = dx;
    if nargout>1
      mu(i,:) = thismu(:);
    end
  end

  if nargout>1
    mu = reshape(mu,[Nt,params.N,params.Nspecies]);
  end
