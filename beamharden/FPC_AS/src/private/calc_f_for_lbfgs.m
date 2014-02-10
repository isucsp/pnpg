function val = calc_f_for_lbfgs(xT)
    global MT AT cT bb nfe_sub;
    persistent resi xtmp;

    nfe_sub = nfe_sub + 1;
    if isstruct(AT); xtmp=zeros(AT.n,1); xtmp(AT.col)=xT; resi = AT.A*xtmp - bb; else resi = AT*xT - bb; end
    if isempty(MT)
        val = cT*xT + 0.5*norm(resi)^2;
    else
        val = cT*xT + 0.5*(resi'*(MT*resi)); 
    end
end