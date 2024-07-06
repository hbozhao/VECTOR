function [x,y,z] = makegrid(L,N,pos)
  %create a uniform grid given the domain size L and number of grid points in each direction N
  %The domain is divided into N(1)*N(2)*... adjacent cells, x and y are the center of the cells
  %it's very important to not let the boundary point overlap in the periodic domain, otherwise particles get attracted to edges
  %L: [Ly,Lx]
  %N: [Ny,Nx], image coordinate
  %pos: positioning of the grid, 'center' (default) or 'corner'
  %center means the grid is centered at the origin, corner means the corner of the grid is at the origin
  if nargin < 3
    pos = 'center';
  end
  d = length(L);
  for i = 1:d
    thisgrid = linspace(0,L(i),N(i)+1);
    thisgrid = (thisgrid(1:end-1)+thisgrid(2:end))/2;
    if strcmp(pos,'center')
      thisgrid = thisgrid - L(i)/2;
    end
    grid{i} = thisgrid;
  end
  if d==2
    [x,y] = ndgrid(grid{1},grid{2});
  elseif d==3
    [x,y,z] = ndgrid(grid{1},grid{2},grid{3});
  end

end
