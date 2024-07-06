function [ax,figsize] = specifyaxes(pos)
  %create subplots by specify the positions of axes with a cell pos like table in latex
  %going from top to bottom, each entry is a also a cell that corresponds to a row
  %{'row_spacing',width} is used to specify a spacing between rows with the specified width
  %{'spacing',width,'axes',[height,width],'spacing',width,...,'alignment',alignment}, to specify from left to right the axes and spacing on this row, 'alignment' specifies the type of alignment ('top','middle','bottom'), alignment can be place anywhere except for the first one (recommended to put at the end)
  %all dimensions are relative
  %axes position refers to the limits defined in the axes.Position property
  %returns the axes handles in the order from top to bottom, left to right
  %figsize is the size of the figures [width,height] (consistent with the convention in figure positions)
  %the position of the axes in the original length scale is kept in ax.UserData

  %compute row positions and total width
  y = 0; %current vertical position
  %save the top and bottom of positions of the row
  rowpos = cell(1,length(pos));
  %save the alignment of the row
  rowalign = cell(1,length(pos));
  figwidth = 0;
  figheight = 0;
  for i = 1:length(pos)
    if isequal(pos{i}{1},'row_spacing')
      thisheight = pos{i}{2};
      y = y + thisheight;
    else
      j = 1;
      thiswidth = 0;
      thisheight = 0;
      while j<=length(pos{i})
        if isequal(pos{i}{j},'alignment')
          align = pos{i}{j+1};
        elseif isequal(pos{i}{j},'spacing')
          thiswidth = thiswidth + pos{i}{j+1};
        else
          x = pos{i}{j+1};
          thiswidth = thiswidth + x(2);
          thisheight = max(thisheight,x(1));
        end
        j = j+2;
      end
      rowpos{i}(1) = y;
      rowpos{i}(2) = y + thisheight;
      rowalign{i} = align;
      y = y + thisheight;
      figwidth = max(figwidth,thiswidth);
    end
    figheight = figheight + thisheight;
  end
  figsize = [figwidth,figheight];

  y = 0;
  counter = 0; %counter for axes index
  scale = [figwidth,figheight,figwidth,figheight];
  for i = 1:length(pos)
    if isequal(pos{i}{1},'row_spacing')
      continue;
    end
    j = 1;
    x = 0; %current horizontal position
    while j<=length(pos{i})
      if isequal(pos{i}{j},'spacing')
        x = x + pos{i}{j+1};
      elseif isequal(pos{i}{j},'axes')
        sz = pos{i}{j+1};
        %get the position of the bottom of the axes
        switch rowalign{i}
        case 'top'
          bottom = rowpos{i}(1) + sz(1);
        case 'middle'
          bottom = mean(rowpos{i}) + sz(1)/2;
        case 'bottom'
          bottom = rowpos{i}(2);
        end
        axpos = [x,figheight-bottom,flip(sz)];
        counter = counter + 1;
        ax(counter) = axes('Position',axpos./scale,'UserData',axpos);
        x = x + sz(2);
      end
      j = j+2;
    end
  end
