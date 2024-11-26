function [pw] = build_pw_struct(oo_,M_,options_)	 
	 
	global M_;

    
	init_filename = M_.fname;

	%%% check whether system has good steady state
	myres=resid(1);
	if sum(isnan(myres)) > 0 
		disp(myres);
		error('Wrong steady state');
	end  
	pw.M_ = M_;
	pw.oo_ = oo_;
	pw.options_ = options_;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET PQ PER DYNARE FILE
	%% PARSING DYNARE FILE
	% look for constraint definitions
	i_bind = []; % line number
	num_bind = []; % constraint id
	regstr='@#define\s?c[1-9]_bind\s?=\s?[01]';
	% find model equations
	i_n = 0;
	% Read txt into cell file00
	fid = fopen([init_filename '.mod'],'r');
	i = 1;
	tline = fgetl(fid);
	file00{i} = tline;
	while ischar(tline)
		if regexp(tline,regstr) == 1
			i_bind=[i_bind;i];
			x = 0;
			splstr=strsplit(tline,'_bind');
			num_bind = [num_bind;str2num(splstr{1}(end:end))];
		end
		% if in dynamics equations: look for 'end;'
		if i_n == 1
			if regexp(tline,'\s?end\s?\;') == 1
				% if found the end, keep line number
				i_n = i;
			end
		end	
		% activate search of end 'model;' if we are in dynamic equations
		if ~isempty(regexp(tline,'\s?model\s?\;')) || ~isempty(regexp(tline,'\s?model(linear)\s?\;'))
			i_n = 1;
		end
		i = i+1;
		tline = fgetl(fid);
		file00{i} = tline;
	end
	fclose(fid);

	pw.Nr 				= length(i_bind);
	pw.combinations 	= dec2bin(0:2^pw.Nr-1) - '0'; % the - '0' turns the dec2bin result into a double

	% write the dynare file we are gonna use
	filexx = file00; % make a copy
	% start writing
	fid = fopen([init_filename '_Occbin.mod'], 'w');
	for i = 1:i_n % stop at end of dynamic equations
		if filexx{i+1} == -1
			fprintf(fid,'%s', file00{i});
			break
		elseif sum(i==i_bind) == 0 % if not constraint definition, then just write
			fprintf(fid,'%s\n', file00{i});
		end
	end
	fclose(fid);	


	%%% LOAD EQUATIONS OF EACH SUB-MODEL
	for i1=1:size(pw.combinations,1)
		
		fname = sprintf('%i',pw.combinations(i1,:));
		
		if i1 > 1 % diet simulation
			command_to_paste = ['dynare ' init_filename '_Occbin.mod console nostrict noclearall minimal_workspace nowarn nolog'];
		else % normal simulation
			command_to_paste = ['dynare ' init_filename '_Occbin.mod nostrict noclearall '];	
		end
		
		% add constraints definition inside commands
		for j=1:length(num_bind)
			command_to_paste = [ command_to_paste sprintf(' -Dc%i_bind=%i',num_bind(j),pw.combinations(i1,j)) ];
        end	
		eval(command_to_paste);
        M_.lead_lag_incidence
        copyfile(['+' init_filename '_Occbin'],['+' init_filename '_' fname])
            
       if i1 == 1

			%% Since state space (t+1,t,t-1) implies to have 3*n_variables
			%% the steady state must be replicated 3 times for each time horizon
			myss = nan(max(max(pw.M_.lead_lag_incidence)),1);
			for i2=1:size(myss,1)
				[~, col] = find(pw.M_.lead_lag_incidence == i2);
				myss(i2) = pw.oo_.dr.ys(col);
			end


			% GET POLICY FUNCTION
			nvar 	= pw.M_.endo_nbr;
			nstate 	= length(pw.oo_.dr.state_var);
			ghx  	= pw.oo_.dr.ghx(:,:);
			ghu  	= pw.oo_.dr.ghu(:,:);
			pw.QQ = ghu;
			% build PP
			jj = 0;
			for i1 = 1:nvar
				% if state variable, pick in ghx
				if sum( pw.oo_.dr.order_var(i1) ==  pw.oo_.dr.state_var)
					jj = jj + 1;
					PPnew = ghx(:,jj);
				else
					PPnew = zeros(nvar,1);
				end
				if i1 == 1
					pw.PP = PPnew;
				else
					pw.PP = [pw.PP PPnew];
				end
			end
		end
		
		eval(['[pw.FF' fname ',pw.GG' fname ',pw.HH' fname ',pw.MM' fname ',pw.CC' fname ',pw.Jac' fname '] = read_pq(myss,''' init_filename '_Occbin'',pw.oo_,pw.M_,pw.options_);'])
	end

	M_ = pw.M_;
	oo_ = pw.oo_;
	options_ = pw.options_;

	M_.fname = init_filename;


    % check if inversion is working
    fname = sprintf('%i',pw.combinations(1,:));
	if sum(sum(abs(eval(['-inv(pw.FF' fname '*pw.PP+pw.GG' fname ')*pw.HH' fname '-pw.PP'])))) > .001 && sum(sum(abs(eval(['-inv(pw.FF' fname '*pw.PP+pw.GG' fname ')*pw.MM' fname '-pw.QQ'])))) > 0.001
		error('Inversion Failed')
	end
	
	try
		delete([init_filename '_Occbin_dynamic.m'])
		delete([init_filename '_Occbin_results.mat'])
		delete([init_filename '_Occbin_set_auxiliary_variables.m'])
		delete([init_filename '_Occbin_static.m'])
		delete([init_filename '_Occbin.m'])
		delete([init_filename '_Occbin.mod'])
		delete([init_filename '_Occbin.log'])
	catch
	end

end

