function y = c2rgb(y,params,channel,channelrgb,clim)
  %convert concentration fields to rgb
  %params contains the dimensionality info
  %y should be a struct that contains concentration fields, rgb will be added to y
  %concentration field has the dimension [params.N,Ntime]
  %use channel to specify three components that you want to output an RGB image of, that is, if rgb = {'A','B','C'}, the output will contain a field that has the dimension of [params.N, 3, Ntime], the color of the three channels will be red, green and blue. If number of components provided is less than 2 and channelrgb, they will use channels in the order of red, green and then blue
  %by default rgb represents three channels: red, green and blue, but user can also specify their own colors for each of the channels as a matrix channelrgb, each row is the rgb value of that channel. By defining their own channelrgb, user can also then specify more than three channel
  %user can set the clim for each channel, each row of clim is a two-element array that corresponds to each channel. Values beyond clim are saturate. clim by default is [0,1] for all channels.
  %like said above, channelrgb can be [# of channels, rgb]. But it can also be [Ntime, # of channels, rgb] so that each time can have a different color
  %the formula for rgb is
  %rgb = sum_i { scale(y.(channel{i})) * channel(i,:)}
  %where scale is (y-clim(i,1))/(clim(i,2)-clim(i,1)) and cutting off values beyond 0 and 1

  Nchannel = length(channel);
  if nargin<5 || isempty(clim)
    clim = repmat([0,1],Nchannel,1);
  end
  nd = length(params.N);
  y.rgb = [];
  for i = 1:length(channel)
    val = y.(channel{i});
    val = (val - clim(i,1))/(clim(i,2)-clim(i,1));
    val(val<0) = 0;
    val(val>1) = 1;
    y.rgb = cat(nd+2,y.rgb,val);
  end
  if nargin<4 || isempty(channelrgb)
    %fill to 3 channels
    if Nchannel<3
      multiplier = zeros(3,1);
      multiplier(1:Nchannel) = 1;
      multiplier = permute(multiplier,[1+(1:nd),1]);
      y.rgb = y.rgb .* multiplier;
    end
  else
    channelrgb = cast(channelrgb,class(y.rgb));
    if length(size(channelrgb))<3
      y.rgb = squeeze(sum(y.rgb .* permute(channelrgb,[2+(1:(nd+1)),1,2]), nd+2));
    else
      y.rgb = squeeze(sum(y.rgb .* permute(channelrgb,[3+(1:nd),1,2,3]), nd+2));
    end
  end
  if ndims(y.rgb)>nd+1
    %ndims(y.rgb)=nd+1 if Ntime=1, then no need to permute
    y.rgb = permute(y.rgb,[1:nd,nd+2,nd+1]);
  end

end
