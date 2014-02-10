function grad = calc_g_for_lbfgs(xT)
    global MT AT cT bb nge_sub;
    persistent resi xtmp cor

    nge_sub = nge_sub + 1;
    % there is a way to save this calculate if resi is already calculated in calc_f_for_lbfgs
    if isstruct(AT);
        xtmp=zeros(AT.n,1); xtmp(AT.col)=xT; resi = AT.A*xtmp - bb;
        if isempty(MT)
            cor = (AT.A')*resi;
        else
            cor = (AT.A')*(MT*resi);
        end
        grad = cT' + cor(AT.col);
    else
        resi = AT*xT - bb;
        if isempty(MT)
            grad = (cT + resi'*AT)';
        else
            grad = (cT + resi'*MT*AT)';
        end
    end
end