
function [ur,ur_rmse]=bisection(func,cond,ul,ur,thresh)
  if(~exist('thresh','var')) thresh=eps; end
  ur_rmse=0; ul_rmse=0;
  while(ur-ul>1e-5*ur)
    fprintf('%10g(%g) <-> %10g(%g)\n',ul,ul_rmse,ur,ur_rmse);
    u=(ur+ul)/2;
    out=func(u);
    if(isstruct(out))
      out=out.x;
    end
    rmse=cond(out);
    if(rmse<=thresh)
      ur=u; ur_rmse=rmse;
    else
      ul=u; ul_rmse=rmse;
    end
  end
  fprintf('u=%g rmse=%g\n',ur,ur_rmse);
end

