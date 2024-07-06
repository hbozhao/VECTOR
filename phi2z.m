function [z,dfdy] = phi2z(phi,flag)
  %a multicomponent generalization of conversion from fraction phi to z which are unbounded
  %each column of phi is a species, sum of phi is 1
  %flag = -1 : phi has N-1 components (default)
  %flag = 0  : phi has N components
  %dfdy is the jacobian dz/dphi (not including phi of the Nth component, phi are linearized by column major)
  if nargin<2
    flag = -1;
  end
  if flag==-1
    phi = [phi, 1-sum(phi,2)];
  end
  N = size(phi,2);
  z = log(cumsum(phi(:,1:N-1),2) ./ phi(:,2:N));
  if nargout>1
    %compute the jacobian, the variables are linearized (column major)
    %jacobian does not include the phi of the Nth component
    ind = 1:numel(phi);
    ind = reshape(ind,size(phi));
    csphi = cumsum(phi(:,1:N-1),2);
    n = size(phi,1)*(N-1);
    dfdy = sparse(n,n);
    for i = 1:N-1
      for j = 1:i
        dfdy = dfdy + sparse(ind(:,i),ind(:,j),1./csphi(:,i),n,n);
      end
      if i<N-1
        dfdy = dfdy + sparse(ind(:,i),ind(:,i+1),-1./phi(:,i+1),n,n);
      else
        for j = 1:N-1
          dfdy = dfdy + sparse(ind(:,i),ind(:,j),1./phi(:,i),n,n);
        end
      end
    end
  end
end
