function varargout = customizeFunGrad(customize,name,request,varargin)
  %allow func to output both fun and grad eval
  %request can be 'fun','grad',or 'fungrad', and the output will be the requested evaluation in the corresponding order.
  %evaluate fun and grad from func function at the same time if grad doesn't exist
  %Another form of struct accepted: fungrad, in which case the function handle accepts varargin{:},params and request
  funstr = customize.(name);
  params = funstr.params;
  if isfield(funstr,'fungrad')
    [fun,grad] = funstr.fungrad(varargin{:},params,request);
  elseif (~isfield(funstr,'grad'))
    if isequal(request,'fun')
      fun = funstr.func(varargin{:},params);
    else
      [fun,grad] = funstr.func(varargin{:},params);
    end
  else
    if isequal(request,'fun') || isequal(request,'fungrad')
      fun = funstr.func(varargin{:},params);
    end
    if isequal(request,'grad') || isequal(request,'fungrad')
      grad = funstr.grad(varargin{:},params);
    end
  end
  switch request
  case 'fun'
    varargout = {fun};
  case 'grad'
    varargout = {grad};
  case 'fungrad'
    varargout = {fun,grad};
  end
end
