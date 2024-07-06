function [y,dy,Iy] = ChemPotential_Flory_Huggins(phi,params,dyflag)
  %general Flory-Huggins model
  %params.r, size of all species, including solvent
  %set params.r(i) = inf to exclude entropy term for species i
  %params.chi, the chi matrix of all species in the same order as params.r
  %params may contain phi_coeff, which is a vector for determining the reference species volume fraction, which is 1 - phi * phi_coeff(:), where each column of phi is a non-reference species volume fraction. normally phi_coeff is all ones, but in some special circumstances user may want to override that (e.g. filament_hex_local_func dimerization). if provided, the length of phi_coeff should be params.Nspecies
  %phi is the volume fraction
  %the last dimension of phi is component. Note that solvent volume fraction is not included
  %y is the diffusional chemical potential, that is, the solvent is a dependent species, defined as dg/dphi_i - dg/dphi_s
  %y has the same dimension as phi
  %by default dyflag = 'full', dy is the full Jacobian
  %user can set dyflag = 'info', dy will be the info to be used by Jacobian_multiplier
  %Iy is the free energy
  %this function also supports params.mu_add_on, which means user can add additional terms to the chemical potential / free energy. The usage is [y,dy,Iy] = params.mu_add_on(phi), where phi includes the solvent species and is linearized for each species, that is, each column is all the phi of ones species, but the returned y, and dy should contain those except the solvent species
  %user can also set params.mu_add_on_params will be passed to params.mu_add_on as additional parameters, [y,dy,Iy] = params.mu_add_on(phi,params.mu_add_on_params{:})
  r = params.r(:).';
  chi = params.chi;
  sz = size(phi);
  ns = length(r); %number of species
  phi = reshape(phi,[],ns-1);
  %add solvent fraction
  if isfield(params,'phi_coeff')
    phi = [phi, 1-phi*params.phi_coeff(:)];
  else
    phi = [phi, 1-sum(phi,2)];
  end
  %partial derivatives of free energy without normalization constraint
  dg = phi * chi';
  entr = ~isinf(r);
  if any(entr)
    dg(:,entr) = dg(:,entr) + (log(phi(:,entr))+1) ./ r(entr);
  end
  y = dg(:,1:(ns-1)) - dg(:,ns);
  y = reshape(y,sz);

  if nargout>1
    if nargin>2 && isequal(dyflag,'info')
      if any(entr)
        dy = 1./phi(:,entr) ./ r(entr);
      else
        dy = [];
      end
    else
      %diagonal element of the solvent and chi matrix
      dy = zeros(size(phi,1),(ns-1)^2);
      if entr(ns)
        if isfield(params,'phi_coeff')
          dy = dy + 1./phi(:,ns)/r(ns) .* repelem(params.phi_coeff(:).',ns-1);
        else
          dy = dy + 1./phi(:,ns)/r(ns);
        end
      end
      dy = dy + params.chi_reduced(:).';
      %diagonal element of all other species
      %linear index of digonal elements
      ind = ns*(1:(ns-1)) - (ns-1);
      sub = 1:(ns-1);
      entr = ~isinf(r(sub));
      ind = ind(entr);
      sub = sub(entr); %nonreference species with entropy
      if ~isempty(ind)
        dy(:,ind) = dy(:,ind) + 1./phi(:,sub)./r(sub);
      end
      dy = reshape(dy,[sz,ns-1]);
    end

    if nargout>2
      logphi = log(phi);
      logphi(phi==0) = 0;
      Iy = sum(phi.*logphi./r,2) + 1/2 * sum(phi .* (phi * chi'),2);
    end
  end

  if isfield(params,'mu_add_on')
    if isfield(params,'mu_add_on_params')
      mu_add_on_params = {params.mu_add_on_params};
    else
      mu_add_on_params = {};
    end
    switch nargout
    case 1
      yadd = params.mu_add_on(phi,mu_add_on_params{:});
    case 2
      [yadd,dyadd] = params.mu_add_on(phi,mu_add_on_params{:});
    case 3
      [yadd,dyadd,Iyadd] = params.mu_add_on(phi,mu_add_on_params{:});
    end
    yadd = reshape(yadd,sz);
    y = y + yadd;
    if nargout>1
      dyadd = reshape(dyadd,[sz,ns-1]);
      dy = dy + dyadd;
    end
    if nargout>2
      Iy = Iy + Iyadd;
    end
  end

end
