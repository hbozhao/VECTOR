function [phi,gamma,psi,flag] = phaseEqm(params,x0,psi,options,grad)
  %find equilibrium phases given the average composition x0 (see LLPS.pdf -> numerical methods -> energy minimization)
  %provide params which contains customize function struct 'mu' or 'ChemicalPotential'. The third
  %provide avg x0, which is a (M-1) vector, where M is the number of components
  %provide initial guess (psi), a (N-1)*M array, where N is the number of phases, M is the number of components, guess_ij is the fraction of phase i that component j is in (probability of being in phase i given it's component j)
  %options: for fminunc
  %grad: by default true. Indicate whether the free energy model params.mu.func also outputs chemial potential in the first output
  %output:
  %phi: fraction of component j in phase i (N*M)
  %gamma: fraction of each phase (N*1)
  %psi: fraction of phase i that component j is in (N*M)
  %use can also provide options for fminunc, we will add Display='none' on top of user provided option
  x0 = x0(:).';
  x0 = [x0,1-sum(x0)];
  if isfield(params,'mu')
    funstr = params.mu;
  elseif isfield(params,'ChemicalPotential')
    funstr = params.ChemicalPotential;
  end
  psi = phi2z(psi')';
  if nargin>3 && ~isempty(options)
    options = optimoptions(options,'Display','none');
  else
    options = optimoptions('fminunc','Display','none');
  end
  if nargin<5 || isempty(grad)
    grad = true;
  end
  options = optimoptions(options,'SpecifyObjectiveGradient',grad);
  if grad
    fun = @(x) eqn_with_grad(funstr,x0,x);
  else
    fun = @(x) eqn(funstr,x0,x);
  end
  [psi,~,flag] = fminunc(fun,psi(:),options);
  if flag>0
    psi = reshape(psi,[],length(x0));
    psi = z2phi(psi',0)'; %N*M
    gamma = sum(x0.*psi,2);
    phi = x0.*psi ./ gamma;
  else
    phi = NaN;
    gamma = NaN;
    psi = NaN;
  end

end

function y = eqn(funstr,x0,psi)
  psi = reshape(psi,[],length(x0));
  psi = z2phi(psi',0)'; %N*M
  gamma = sum(x0.*psi,2);
  phi = x0.*psi ./ gamma;
  phi = phi(:,1:end-1); %N*(M-1)
  [~,~,g] = funstr.func(phi,funstr.params);
  y = sum(gamma.*g);
end


function [y,dy] = eqn_with_grad(funstr,x0,psi)
  M = length(x0);
  psi = reshape(psi,[],length(x0));
  Npsi = numel(psi);
  N = size(psi,1)+1;
  [psi,dpsi] = z2phi(psi',0);
  psi = psi'; %N*M
  dpsi = full(dpsi);
  dpsi = reshape(dpsi,M,N-1,M,N-1); %dpsi here based on psi' which has the size of M*(N-1), so reshape into M*(N-1)*M*(N-1), where the first two indices is output of z2phi, the last two indices are the input of z2phi
  dpsi = cat(2, dpsi, -sum(dpsi,2)); %get M*N*M*(N-1)
  dpsi = permute(dpsi,[2,1,4,3]); %N*M*(N-1)*M, now first two indices is output of z2phi transposed, and the second indices are the input of z2phi transposed
  dpsi = reshape(dpsi,N,M,Npsi); %N*M*((N-1)*M), the last index is the input of eqn_with_grad
  gamma = sum(x0.*psi,2);
  dgamma = sum(x0.*dpsi,2);
  phi = x0.*psi ./ gamma;
  dphi = x0.*dpsi ./ gamma - x0.*psi ./ gamma.^2 .* dgamma;
  phi = phi(:,1:end-1); %N*(M-1)
  dphi = dphi(:,1:end-1,:);
  [mu,~,g] = funstr.func(phi,funstr.params);
  dg = sum(mu .* dphi, 2);
  dg = squeeze(dg);
  y = sum(gamma.*g);
  dy = sum(gamma.*dg + squeeze(dgamma).*g,1);
  if any(~isfinite(dy))
    %if dy returns NaN, roots in linesearch will throw error
    y = NaN;
  end
end
