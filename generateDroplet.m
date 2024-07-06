function varargout = generateDroplet(x,y,center,d,gamma,phi,phim,periodic,weight)
  %generate the concentration (or phase) field of one or more droplets and one or more phases that has a target phase fraction gamma
  %the length of gamma n is the total number of phases - 1 (the background phase is not included)
  %the order of gamma is from inner to outer layer of the droplet
  %the output g will contain n channels [Nx,Ny,n]
  %optional: phi is the phase composition [n+1,m], where m is the number of components
  %if phi is provided, we output [c,g], where c is the composition of all species [Nx,Ny,m], the order of rows of phi is from inner to outer layer of the droplet
  %if phi is not provided, we output g
  %if gamma is not provided, please provide phim and phi (if gamma provided, phim is ignored), which is the mean composition. In this case, we require that the total number of species is n+1, phi should be [n+1,n] (otherwise no sln), the last mean concentration of the last species is added automatically to satisfy normalization, and similarly phim should be [1,n] or [n,1], we will solve for gamma. The output will be [c,g], size of c is [Nx,Ny,n] (the last species is not included), size of g is still [Nx,Ny,n]
  %c = g * phi (matrix product) c_i = g_j * phi_ji (here g is expanded to [Nx,Ny,n+1] to include the last phase)
  %for n=1, g is defined to be:
  % sum_i{ 1 / (1+exp(|r-center(i,:)|-r0)/d) }, where r is the location, center(i,:) is the origin of the ith droplet and d is the width of the droplet. g is ~1 inside the droplet and 0 outside
  %user must make sure that the droplet must be the minority phase.
  %r0 is the unknown to be solved for. Droplets have the same size
  %user can either input the major and minor phase concentration (one component only) c and average composition cm, in which case varargin = {c,cm}, varargout = {c,g}
  %or input the fraction of the minor phase gamma, in which case varagin = {gamma}. varargout = {g}
  %the equation to be solved is the target minor phase fraction: gamma = mean(g(:))
  %the concentration field is
  % c = (c2-c1) * g + c1, or better: c2*g + c1*(1-g)
  % where c1 is the major phase concentration, c2 is the minor phase concentration
  %input center is an array [N*2], each row is a center, [x,y]
  %on 10/5/2022, I generalized this to multiple phases, in which case we only support varargin being gamma, which is an array of the fraction of each phase
  %each phase is represented by a field, the output g is a [Nx,Ny,n] array, where n is the number of phases. g(x,y,i)=1 indicates the presence of phase i at (x,y)
  %multiple droplet is still supported
  %the inner most layer is
  % 1 / (1+exp((|r-center|-r_1)/d))
  %the nth layer is
  % 1 / (1+exp((|r-center|-r_n)/d)) - 1 / (1+exp((|r-center|-r_(n-1))/d))
  %r_1, ..., r_n will be solved
  %this definition can it very easy to code, because we just need to generate 1 / (1+exp((|r-center|-r_i)/d)) and solve for r_i such that the average is equal sum(gamma(1:i)) for i = 1:n, separately. And then substrate to get each field
  %on 10/6/2022 center is generalized to allow for different center positions for the different phases (doesn't have to be concentric, but user is responsible for making sure the inner phase is completely contained in outer phase)
  %to enable this make center a [N*2*n] array, where N is the number of droplets, n is the number of phases in each droplet
  %on 10/30/2022, I added an option periodic. If provided and set to true, (false by default), we will use periodic boundary for determining r
  %weight: if provided, the output c = (c2*g + c1*(1-g))*weight. The output g remains the same as before, not multiplied by weight. Without weight, the definition of gamma is the mean of g, with weight, the definition of gamma is the weighted mean of g, or sum(g.*weight)./sum(weight).
  %if provided, weight must have the same size as x
  %on 7/11/2023, I extended this to nD. To generate droplet in arbitrary dimension, provide {x,y,...} a cell of the coordinates in place of x, when x is a cell, y will be ignored
  %on 7/12/2023, in varargout, append another output r0, the solution to the equation so varargout is {c,g,r0} or {g,r0}

  no_gamma = isempty(gamma);
  if no_gamma
    Nphase = size(phi,1)-1;
    phi_extend = [phi,1-sum(phi,2)];
    if length(phim)~=Nphase
      error('phim length wrong');
    end
    gamma = phi_extend \ [phim(:);1-sum(phim)];
    gamma = gamma(1:end-1);
  end
  if ~iscell(x)
    x = {x,y};
  end
  dim = length(x);
  L = cellfun(@(x) max(x(:))-min(x(:)),x);
  if nargin>7
    if isscalar(periodic)
      if periodic
        %turn periodic to the domain size
        periodic = L;
      else
        periodic = [];
      end
    end
  else
    periodic = [];
  end
  if nargin<9
    weight = [];
  end
  N = size(center,1);
  Nphase = length(gamma);
  Ngrid = size(x{1});
  %guess for r_i, presumably it is rectangular grid but doesn't have to be
  area = prod(L);
  gamma = cumsum(gamma);
  switch dim
  case 2
    r0 = sqrt(area*gamma/N/pi);
  case 3
    r0 = (area*gamma/N*3/4/pi).^(1/3);
  end
  g = zeros([prod(Ngrid),Nphase]);
  prevg = zeros(Ngrid);

  opt = optimoptions('fsolve','Display','off');
  for i = 1:Nphase
    if ndims(center)==3
      thiscenter = center(:,:,i);
    else
      thiscenter = center;
    end
    [r0(i),~,flag] = fsolve(@(r0) eqn(x,thiscenter,d,r0,gamma(i),periodic,weight), r0(i), opt);
    if flag<=0
      thisg = NaN(Ngrid);
    else
      thisg = gfunc(x,thiscenter,d,r0(i),periodic);
    end
    g(:,i) = thisg(:) - prevg(:);
    prevg = thisg;
  end
  g_extended = [g,1-sum(g,2)];
  g = reshape(g,[Ngrid,Nphase]);
  if (nargin>5 && ~isempty(phi)) || no_gamma
    c = g_extended * phi;
    if ~isempty(weight)
      c = c.*weight(:);
    end
    c = reshape(c,[Ngrid,size(c,2)]);
    varargout = {c,g,r0};
  else
    varargout = {g,r0};
  end

end

function g = gfunc(x,center,d,r0,periodic)
  g = 0;
  for i = 1:size(center,1)
    r = 0;
    for j = 1:length(x)
      dx = x{j}-center(i,j);
      if ~isempty(periodic)
        dx = dx - periodic(j).*round(dx./periodic(j));
      end
      r = r + dx.^2;
    end
    r = sqrt(r);
    g = g + 1./(1+exp((r-r0)/d));
  end
end

function y = eqn(x,center,d,r0,gamma,periodic,weight)
  g = gfunc(x,center,d,r0,periodic);
  if isempty(weight)
    y = mean(g(:)) - gamma;
  else
    y = sum(g(:).*weight(:)) ./ sum(weight(:)) - gamma;
  end
end
