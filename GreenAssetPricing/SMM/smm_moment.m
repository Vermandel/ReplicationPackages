function output = smm_moment(param_file,M_,oo_,options_,wcc,shocks,pp)

    load([param_file '_smm_temp.mat']);
    
    for ix=1:length(params)
        M_.params(strmatch(deblank(ve_names{ix}),M_.param_names,'exact')) = params(ix);
%        eval([deblank(ve_names{ix}) '= '  num2str(params(ix))])
    end

    if nargin > 6 && ~isempty(pp)
        for ix=1:2:size(pp,2)
            M_.params(strmatch(pp{ix},M_.param_names,'exact'))=pp{ix+1};
        end
     end


    [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, {});


    nt = size(smm.random_shocks{1},1);
    nd = size(smm.random_shocks,2);
    y0 = oo_.dr.ys;
    if isfield(M_,'histval_dseries')
      for ix=find(M_.histval_dseries('0Y').data~=0)
	       idx = strmatch(M_.histval_dseries.name{ix},M_.endo_names,'exact');
	       y0(idx) = M_.histval_dseries('0Y').data(ix);
      end
    end    
    y_sims   = nan(M_.endo_nbr,nt+1,nd);
    y_moms   = nan(M_.endo_nbr,3,nd);
    y_betas  = nan(M_.endo_nbr,M_.endo_nbr,nd);
    for ix=1:nd
         y_sims(:,:,ix) = simult_(M_, options_, y0, oo_.dr, smm.random_shocks{ix}, options_.order);
         y_moms(:,1,ix)  = mean(y_sims(:,(smm.nb_drops+1):end,ix),2);
         y_moms(:,2,ix)  = std(y_sims(:,(smm.nb_drops+1):end,ix),[],2);
         for iy = 1:M_.endo_nbr
            themat =  corrcoef(y_sims(iy,(smm.nb_drops+2):end,ix),y_sims(iy,(smm.nb_drops+1):end-1,ix));
            y_moms(iy,3,ix) =  themat(1,2);
            % inv(XX'*XX)*XX'*YY
            warning off;
            for iyx = 1:M_.endo_nbr
                YY = y_sims(iy,(smm.nb_drops+1):end,ix)';
                XX = [ ones(nt+1-smm.nb_drops,1) y_sims(iyx,(smm.nb_drops+1):end,ix)'];
                betas = inv(XX'*XX)*XX'*YY;
                if isfinite(betas(2))
                y_betas(iyx,iy,ix) = betas(2);
                end
            end
         end
    end

    output.mom = mean(y_moms,3);
    output.var_names = M_.endo_names;
    output.ve_names = ve_names;
    output.param = params;
    output.y_beta = mean(y_betas,3);

    if exist('x_SD') == 0
        x_SD = nan(size(params));
    end
    output.x_SD   = x_SD;
    output.smm    = smm;
    output.y_sims = y_sims;
    output.y_det =  simult_(M_, options_, y0, oo_.dr, smm.random_shocks{ix}*0, options_.order);
    
    if nargin > 4 && ~isempty(wcc)
		tolerance = 1e-06;
	    options = optimset('Display','on','TolFun',tolerance,'TolX',tolerance,'TolFun',tolerance,'DiffMinChange',tolerance,'Algorithm','Interior-Point','MaxFunEvals',2500,'MaxIter',2500);
        u_stoch = output.mom(strmatch(wcc{1},M_.endo_names,'exact'),1);
        output.wc = (1-fminsearch('get_welfcost',0.996,[],M_,oo_,options_,smm,wcc{2},wcc{1},u_stoch,y0))*100;
        ln_u=nan(size(output.y_sims,3),1);
        for ix=1:size(output.y_sims,3)
            ln_u(ix) = mean(output.y_sims(strmatch(wcc{1},M_.endo_names,'exact'),(smm.nb_drops+1):end,ix)./output.y_det(strmatch(wcc{1},M_.endo_names,'exact'),(smm.nb_drops+1):end));
        end
        output.ln_u = mean(ln_u);
    end
    
    if nargin > 5 && ~isempty(shocks)
        std_y = M_.params(strmatch('std_y',M_.param_names,'exact'));
        they0 = mean(y_sims(:,(1+smm.nb_drops),:),3);
        output.y_simul3 =  simult_(M_, options_, they0, oo_.dr, shocks/std_y, options_.order);
        % order 2
        options_.order = 2;
        [oo_.dr,infos,M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
        for ix=1:nd
             y_sims(:,:,ix) = simult_(M_, options_, y0, oo_.dr, smm.random_shocks{ix}, options_.order);
        end
        they0 = mean(y_sims(:,(1+smm.nb_drops),:),3);
        output.y_simul2 =  simult_(M_, options_, they0, oo_.dr, shocks/std_y, options_.order);
        % order 1
        options_.order = 1;
        [oo_.dr,infos,M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
        for ix=1:nd
             y_sims(:,:,ix) = simult_(M_, options_, y0, oo_.dr, smm.random_shocks{ix}, options_.order);
        end
        they0 = mean(y_sims(:,(1+smm.nb_drops),:),3);
        output.y_simul1 =  simult_(M_, options_, they0, oo_.dr, shocks/std_y, options_.order);
    end


end
