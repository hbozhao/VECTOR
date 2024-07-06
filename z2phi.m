function [phi,dfdy] = z2phi(z,flag)
  %a multicomponent generalization of conversion from z to fraction phi
  %each column of z and phi is a species
  %flag = -1 : output phi has N-1 components (default), the last component phi_N is omitted
  %flag = 0  : output phi has N components
  %dfdy is the jacobian dphi/dz (not including phi of the Nth component, phi are linearized by column major)
  if nargin<2
    flag = -1;
  end
  z = 1./(1+exp(-z));
  %within diff: 0, phi_1, phi_1+phi_2, ..., phi_1+...+phi_(N-1)
  Nc = size(z,1);
  phi = diff([zeros(Nc,1), cumprod(z,2,'reverse')],[],2);
  if flag == 0
    phi = [phi, 1-sum(phi,2)];
  end

  if nargout>1
    %compute the jacobian, the variables are linearized (column major)
    %jacobian does not include the phi of the Nth component
    %xi = 1/(1+exp(-z)), dzdz below is dxi/dz
    dzdz = z.*(1-z);
    ind = 1:numel(phi);
    ind = reshape(ind,size(phi));
    N = size(z,2)+1;
    n = size(phi,1)*(N-1);
    dfdy = sparse(n,n);
    %see LLPS.pdf Numerical methods
    %phi_i = (1-z_(i-1)) * z_i * ... * z_(N-1) if i>1
    %phi_1 = z_1 * ... * z_(N-1)
    for i = (N-1):-1:1
      for j = (N-1):-1:i
        dfdy = dfdy + sparse(ind(:,i),ind(:,j),phi(:,i)./z(:,j),n,n);
      end
      if i>1
        dfdy = dfdy + sparse(ind(:,i),ind(:,i-1),-phi(:,i)./(1-z(:,i-1)),n,n);
      end
    end
    dfdy = dfdy .* dzdz(:).';
  end
end
