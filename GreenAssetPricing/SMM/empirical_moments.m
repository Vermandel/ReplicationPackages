function ICemp = empirical_moments(data,z,varname)
    % compute dates
    IC = 0.95;
    fd = fieldnames(data);
    if any(strcmp(varname,fd))
        YY=eval(['data.' varname]);
        nT = length(YY);
        switch z
            case 1
            % mean
            % 95% interval
            th = tinv([1-(1-IC)/2,(1-IC)/2],nT-1);
	        XX=ones(size(YY));
	        MU=inv(XX'*XX)*XX'*YY;
	        EE=var(YY-MU);
            ICemp = [MU (MU+th(2)*sqrt(diag(EE*inv(XX'*XX)))) (MU+th(1)*sqrt(diag(EE*inv(XX'*XX))))];
            case 2
            % std
            % 95% interval
            ch = chi2inv([1-(1-IC)/2,(1-IC)/2],nT-1);
            ICemp = [ std(YY) sqrt(((nT-1)*std(YY)^2)/ch(1)),sqrt(((nT-1)*std(YY)^2)/ch(2))];
            case 3
	        [thecorr,~,corr_lb,corr_ub] = corrcoef(YY(2:end),YY(1:end-1),'Alpha',1-IC);
            ICemp = [ thecorr(2,1) corr_lb(2,1) corr_ub(2,1)];
        end
    else
        ICemp = 0;
    end

end

