function VECTOR_package(tspan,yraw,snapshots,td,lambda,params,clim,t_on,filepath,study,xy)
  %a wrapper function to visualize the simulations results and save the post-processing figures
  %snapshots are index of key frames to display
  %clim is for limits of tel and core channels used by solverParser, in the order of At and C by order (corelet,telomere)
  %t_on is the duration of light-on
  %generate a package of figures for the droppullet
  %filepath gives the directory of files to be saved
  %study is the study name
  %files that will be saved: (for every pdf there is also a svg)
  %(study).mp4 the entire movie
  %(pp_study_1).pdf, the distance-time plot (version 1, no highlighted time points, distance normalized with respect to 2*R0) and kymograph
  %(pp_study_2).pdf, the distance-time plot (version 2, with keyframe time points highlighted, distance normalized with respect to initial distance) and kymograph
  %(pp_study_3).pdf, similar to 2, with kymograph from C, AB channels separately and together
  %(study_snapshots_i).png, key frames, C and AB channels, in color
  %(study_snapshots_C_i).png, key frames, C channel, in black and white
  %(study_snapshots_AB_i).png, key frames, C channel, in black and white
  %notice that if study contains the string 'asym', then instead of outputing the distance between the two, output the position of the two telomeres (relative to the initial center of mass)
  %on 7/25/2022 I changed the 5th input Rcore to lambda, which is sqrt(params.kappa)
  %xy: 'old' (default), old x y order, 'new', new x y order (x is the first dimension, y is the second)
  if nargin < 11
    xy = 'old';
  end

  linecolor = [175,51,185]/256;
  patchcolor = [38,169,224]/256;
  telcolor = [1,0,1];
  corecolor = [0,1,0];

  asym = contains(study,'asym');

  if iscell(yraw)
    yy = solverParser(cat(1,yraw{:}),params,{'At','C'},[corecolor;telcolor],'clim',clim);
  else
    yy = solverParser(yraw,params,{'At','C'},[corecolor;telcolor],'clim',clim);
  end
  if ~isfield(yy,'add')
    %this is for coilin, use peaks along the mid-plane to determine the position distance
    C = yy.C(round(params.N(1)/2),:,:);
    C = squeeze(C);
    x = linspace(-params.L(2)/2,params.L(2)/2,params.N(2));
    x = x(:);
    dx = params.L(2)/params.N(2);
    TF = islocalmax(C,1,'MaxNumExtrema',4,'MinSeparation',5,'MinProminence',0.01);
    nmax = sum(TF,1);
    Ntime = size(C,2);
    % distance = sum(C.*x,1) ./ sum(C,1);
    distance = zeros(Ntime,1);
    for i = 1:Ntime
      if nmax(i) > 1
        ind = find(TF(:,i));
        ind = ind(1);
        ind = ind+(-10:10);
        distance(i) = sum(C(ind,i).*x(ind)) ./ sum(C(ind,i)) * 2;
      end
    end
    distance(nmax>1) = smoothdata(distance(nmax>1),'movmean',10);
    distance = abs(distance);
    single = false;
  else
    single = size(yy.add,2)==1;
    if ~single
      %if there is only a single telomere. Don't plot distance-time.
      if isequal(params.substrate.dynamics,'Rouse')
        thisyadd = reshape(yy.add,[],1+params.substrate.Rouse.N,2,2);
        thisyadd = squeeze(thisyadd(:,1,:,1));
        distance = diff(thisyadd,[],2);
      else
        distance = abs(diff(yy.add(:,:,1),[],2));
      end
    end
    if asym
      distance = yy.add(:,:,1);
    end
  end
  t_on = [0,t_on];
  channel = {'At','C'};

  for i = 1:3
    fig = figure;
    if i==3 || single
      axpos = {{'row_spacing',0.05},...
              {'spacing',0.15,'axes',[0.2,1],'spacing',0.05,'alignment','top'},...
              {'row_spacing',0.1},...
              {'spacing',0.15,'axes',[0.2,1],'spacing',0.05,'alignment','top'},...
              {'row_spacing',0.1},...
              {'spacing',0.15,'axes',[0.2,1],'spacing',0.05,'alignment','top'},...
              {'row_spacing',0.1},...
              {'spacing',0.15,'axes',[0.7,1],'spacing',0.05,'alignment','top'},...
              {'row_spacing',0.15}};
    else
      axpos = {{'row_spacing',0.05},...
              {'spacing',0.15,'axes',[0.2,1],'spacing',0.05,'alignment','top'},...
              {'row_spacing',0.1},...
              {'spacing',0.15,'axes',[0.7,1],'spacing',0.05,'alignment','top'},...
              {'row_spacing',0.15}};
    end
    if single
      axpos(end-1:end) = [];
    end
    [ax,figsize] = specifyaxes(axpos);
    fig.Position = [100,100,figsize/figsize(1)*500];
    thisax = ax(1);
    axes(thisax);
    switch xy
    case 'old'
      im = permute(yy.rgb(round(params.N(1)/2),:,:,:),[2,4,3,1]);
    case 'new'
      im = permute(yy.rgb(:,round(params.N(1)/2),:,:),[1,4,3,2]);
    end
    if single && yy.add(1,1,1)<0
      %if single telomere, let the telomere sit at the top
      im = flip(im);
    end
    image(tspan/td,params.L(2)/sqrt(params.kappa),im);
    if ~single
      thisax.XTick = [];
    end
    thisax.YTick = [];
    if i==3 || single
      for k = 1:length(channel)
        thisax = ax(k+1);
        axes(thisax);
        switch xy
        case 'old'
          im = permute(yy.(channel{k})(round(params.N(1)/2),:,:),[2,3,1]);
        case 'new'
          im = permute(yy.(channel{k})(:,round(params.N(1)/2),:),[1,3,2]);
        end
        thisclim = clim(k,:);
        im = (im-thisclim(1))/diff(thisclim);
        im(im<0) = 0;
        im(im>1) = 1;
        im = im.*permute(ones(3,1),[2,3,1]);
        if single && yy.add(1,1,1)<0
          im = flip(im);
        end
        image(tspan/td,params.L(2)/sqrt(params.kappa),im);
        if ~single
          thisax.XTick = [];
        end
        thisax.YTick = [];
      end
    end
    if single
      break;
    end
    thisax = ax(end);
    axes(thisax);
    if i==1
      dnorm = lambda;
      if asym
        thisax.YLim = distance(1,:)/dnorm*1.2;
      else
        thisax.YLim = [0,distance(1)/dnorm*1.2];
      end
    else
      if asym
        thisax.YLim = [-1,1];
      else
        thisax.YLim = [0,1.2];
      end
      dnorm = abs(distance(1,:));
    end
    if diff(t_on)>0
      patch([t_on,flip(t_on)]/td,repelem(thisax.YLim,2),patchcolor,'FaceColor',patchcolor,'EdgeColor','none','FaceAlpha',0.3);
    end
    hold on;
    plot(tspan/td,distance./dnorm,'LineWidth',2,'Color',linecolor);
    if i~=1
      hold on;
      plot(tspan(snapshots)/td,distance(snapshots,:)./dnorm,'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
    end
    thisax.XLim = tspan([1,end])/td;
    xlabel('Time ($t/t_d$)','Interpreter','latex');
    if i==1
      ylabel('Distance ($d/\lambda$)','Interpreter','latex');
    else
      ylabel('Distance ($d/d_0$)','Interpreter','latex');
    end
    thisax.TickDir = 'out';
  end
  for i = 1:length(snapshots)
    im = yy.rgb(:,:,:,snapshots(i));
    if isequal(xy,'new')
      im = permute(im,[2,1,3]);
    end
    imwrite(im,[filepath,study,'_snapshots_',num2str(i),'.png']);
    for k = 1:length(channel)
      im = yy.(channel{k})(:,:,snapshots(i));
      if isequal(xy,'new')
        im = permute(im,[2,1]);
      end
      thisclim = clim(k,:);
      im = (im-thisclim(1))/diff(thisclim);
      im(im<0) = 0;
      im(im>1) = 1;
      imwrite(im,[filepath,study,'_snapshots_',channel{k},'_',num2str(i),'.png']);
    end
  end

  vid = VideoWriter([filepath,study],'MPEG-4');
  vid.FrameRate = 30;
  open(vid);
  rgb = yy.rgb;
  rgb(rgb<0) = 0;
  if isequal(xy,'new')
    rgb = permute(rgb,[2,1,3,4]);
  end
  writeVideo(vid,rgb);
  close(vid);

  if isequal(params.substrate.dynamics,'Rouse')
    %if rouse. make another movie that visualize the dynamics of the chain

    yadd = reshape(yy.add,[],1+params.substrate.Rouse.N,2,2); %[Ntime * Npoints * 2 chains * 2 dimensions]
    %only take the first dimension
    yadd = yadd(:,:,:,1);
    %complement the other half of the chain with symmetry condition
    yadd = cat(2,flip(yadd(:,2:end,:,:),2),yadd);
    im_xrange = [min(params(1).grid{1}(:)),max(params(1).grid{1}(:))];
    im_yrange = [min(params(1).grid{2}(:)),max(params(1).grid{2}(:))];
    %the second dimension is constant and make it linear and span the entire y direction * 2, so zoom in 2x
    yadd_y = linspace(im_yrange(1),im_yrange(2),size(yadd,2)) * 5;
    yadd_y = repmat(yadd_y,[size(yadd,1),1,size(yadd,3)]);
    yadd = cat(4,yadd,yadd_y);
    fig = figure;
    axpos = {{'axes',[1,1],'alignment','top'}};
    [ax,figsize] = specifyaxes(axpos);
    fig.Position = [100,100,figsize/figsize(1)*300];
    vid = VideoWriter([filepath,study,'_chain'],'MPEG-4');
    vid.FrameRate = 30;
    open(vid);
    for i = 1:size(yadd,1)
      thisyadd = squeeze(yadd(i,:,:,:));
      image(im_xrange,im_yrange,rgb(:,:,:,i),'Interpolation','bilinear');
      ax.Visible = 'off';
      hold on
      for j = 1:size(thisyadd,3)
        plot(thisyadd(:,j,1),thisyadd(:,j,2),'LineWidth',2,'Color',ones(1,3));
      end
      hold off
      set(gcf,'InvertHardCopy','off')
      writeVideo(vid,print('-RGBImage','-r600'));
    end
    close(vid);
  end