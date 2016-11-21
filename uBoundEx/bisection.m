
function [ur,ur_rmse]=bisection(opt,func,cond,ul,ur,thresh)
  if(~exist('thresh','var')) thresh=eps; end
  ur_rmse=0; ul_rmse=0;
  while(ur-ul>1e-5*ur)
    fprintf('%10g(%g) <-> %10g(%g)\n',ul,ul_rmse,ur,ur_rmse);
    opt.u=(ur+ul)/2;
    out=func(opt);
    if(isstruct(out))
      out=out.alpha;
    end
    rmse=cond(out);
    if(rmse<=thresh)
      ur=opt.u; ur_rmse=rmse;
    else
      ul=opt.u; ul_rmse=rmse;
    end
  end
  fprintf('u=%g rmse=%g\n',opt.u,ur_rmse);
end

