function params = Flory_Huggins_setup(params,flag)
  %constructor function for params for Flory Huggins model simulation
  %fields that are not handled in this construct and that user should supply themselves:
  %fields that will be constructed in this function
  %1. chi: the full chi interaction matrix involving all species
  %2. chi_reduced: the reduced chi matrix excluding the solvent species (obeying normalization constraint)
  %3. mu: the function handle of the chemical potential model
  %the params input to ths function contain:
  %1. species. A cell that contain list of species names. y in the solver will be in this order. The last species MUST be the solvent/reference/dependent species
  %2. r. A list of sizes in the same order as params.species
  %3. chi. This chi can be the full interaction matrix, or a cell that has the format of {'A','B',chi_AB,...}, for all of the pairs, pairs with zero interaction strength can be omitted
  %optionally, params may contain lambda, in which case kappa will be set to (lambda^2*(-chi_reduced))
  %params may contain phi_coeff, which is a vector for determining the reference species volume fraction, which is 1 - phi * phi_coeff(:), where each column of phi is a non-reference species volume fraction. normally phi_coeff is all ones, but in some special circumstances user may want to override that (e.g. filament_hex_local_func dimerization). if provided, the length of phi_coeff should be params.Nspecies. Note that this is only when some species have identical dynamics. The chemical potential of non reference species does not change.
  %set flag = 'multiplier' to use dmu and kappa multiplier (by default not used)

  Ns = length(params.species);
  params.Nspecies = Ns-1;
  if iscell(params.chi)
    params.chi = reshape(params.chi,3,[]);
    chi = zeros(Ns);
    for i = 1:size(params.chi,2)
      A = ismember(params.species,params.chi(1,i));
      A = find(A);
      B = ismember(params.species,params.chi(2,i));
      B = find(B);
      thischi = params.chi{3,i};
      if A==B
        chi(A,B) = 2*thischi;
      else
        chi(A,B) = thischi;
        chi(B,A) = thischi;
      end
    end
    params.chi = chi;
  end
  chi = params.chi;
  if isfield(params,'phi_coeff')
    coeff = params.phi_coeff(:);
    params.chi_reduced = chi(1:(Ns-1),1:(Ns-1)) - ones(Ns-1,1)*chi(Ns,1:(Ns-1)) - chi(1:(Ns-1),Ns)*coeff' + ones(Ns-1,1)*chi(Ns,Ns)*coeff';
  else
    params.chi_reduced = chi(1:(Ns-1),1:(Ns-1)) - ones(Ns-1,1)*chi(Ns,1:(Ns-1)) - chi(1:(Ns-1),Ns)*ones(1,Ns-1) + chi(Ns,Ns);
  end
  names = {'r','chi','chi_reduced','mu_add_on','mu_add_on_params','phi_coeff'};
  %add the relevant fields to muparams
  for i = 1:length(names)
    if isfield(params,names{i})
      muparams.(names{i}) = params.(names{i});
    end
  end
  params.mu = struct('func',@ChemPotential_Flory_Huggins,'params',muparams);
  if isfield(params,'lambda')
    params.kappa = -params.lambda^2 * params.chi_reduced;
  end

end
