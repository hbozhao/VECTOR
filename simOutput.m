function status = simOutput(t,y,flag,progress)
%show the progress of the simulation, time, and ETA (no saving in this version)
persistent t0 tf 
if progress
  switch flag
  case 'init'
    t0 = t(1);
    tf = t(end);
    tic;
  case []
    if length(t)>1
      for i=1:(length(t)-1)
        disp(['Current physical time: ',num2str(t(i))]);
      end
    end
    tnow = t(end);
    runtime=toc;
    disp(['Current physical time: ',num2str(tnow),' Elapsed time (min): ',num2str(runtime/60),' ETA (min): ',num2str(runtime/60*(tf-tnow)/(tnow-t0))]);
  end
end

status=0;

