function yy = solverParser(y,params,channel,channelrgb,varargin)
  %y is the output of the ode solver
  %the second dimension is ordered by points in space then species
  %the output is a struct, each field (with the name of the species) contains the concentration of all species (including the solvent) which has the dimension of [params.N,Ntime]
  %user can also not provide species, in which case the y will be reshaped into [params.N,Nspecies,Ntime] and put in yy.y. rgb will not work without species. metaspecies still works.
  %use channel to specify three components that you want to output an RGB image of, that is, if rgb = {'A','B','C'}, the output will contain a field that has the dimension of [params.N, 3, Ntime], the color of the three channels will be red, green and blue
  %in addition to computing the concentration of species in params.species, user can also provide params.metaspecies to compute linear combinations of species, such as the total concentration of monomers
  %for example: if params.species = {'A','B','AB'}, params.metaspecies.At = [1,0,ra/rab] computes At = [A] + ra/rab*[AB]
  %rgb can also request metaspecies, like {'At','C','S'}
  %by default rgb represents three channels: red, green and blue, but user can also specify their own colors for each of the channels as a matrix channelrgb, each row is the rgb value of that channel. By defining their own channelrgb, user can also then specify more than three channel
  %additional degrees of freedom will be kept in yy.add. If params.substrate exists, yy.add will be in the format of [Ntime, Ns, d], where Ns is the number of centers, otherwise yy.add is in the format of [Ntime,Nadd]
  %if params.add does not exist, all extra degrees of freedom will be put in yadd, otherwise only the first params.add.N extra degrees of freedom will be put in yadd

  %Additional input through varargin (field-parameter pair)
  %clim: user can set the clim for each channel, each row of clim is a two-element array that corresponds to each channel. Values beyond clim are saturate. clim by default is [0,1] for all channels.

  ps = inputParser;
  addParameter(ps,'clim',[]);
  ps.CaseSensitive = false;
  parse(ps,varargin{:});
  ps = ps.Results;


  Ntime = size(y,1);
  if ~isfield(params,'Nspecies')
    params.Nspecies = 1;
  end
  dof = prod(params.N) * params.Nspecies;
  if isfield(params,'add') || size(y,2)>dof
    if ~isfield(params,'add')
      params.add.N = size(y,2)-dof;
    end
    yadd = y(:,dof+(1:params.add.N));
    y = y(:,1:dof);
    if isfield(params,'substrate')
      d = length(params.N);
      yadd = yadd(:,1:(params.substrate.N*d));
      yadd = reshape(yadd,Ntime,params.substrate.N,d);
    end
    yy.add = yadd;
  end
  y = reshape(y,[],params.Nspecies);
  y = y.';
  nospecies = ~isfield(params,'species');
  if nospecies
    params.species = {};
  end
  for i = 1:length(params.species)
    if i==length(params.species)
      thisy = 1-sum(y,1);
    else
      thisy = y(i,:);
    end
    thisy = reshape(thisy,Ntime,[]).';
    thisy = reshape(thisy,[params.N,Ntime]);
    yy.(params.species{i}) = thisy;
  end
  if isfield(params,'metaspecies')
    names = fieldnames(params.metaspecies);
    for i = 1:length(names)
      thisy = params.metaspecies.(names{i}) * y;
      thisy = reshape(thisy,Ntime,[]).';
      thisy = reshape(thisy,[params.N,Ntime]);
      yy.(names{i}) = thisy;
    end
  end
  if nospecies
    y = reshape(y.',Ntime,[]).';
    yy.y = reshape(y,[params.N,params.Nspecies,Ntime]);
  end
  if nargin>2 && ~isempty(channel) && ~nospecies
    if nargin<4
      channelrgb = [];
    end
    yy = c2rgb(yy,params,channel,channelrgb,ps.clim);
  end
