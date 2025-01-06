%addpath('../core_files/')

% get data
data =  readtable('DataFinance.xlsx','Sheet','compact');


% CRRA LAISSEZ-FAIRE
eval(['dynare call_CARA  -DOPTIMAL=0 -DEZ=0 -DDIET=0 ;'])
CRRA_LF  = smm_moment('data_CARA',M_,oo_,options_);
% CRRA OPTIMAL POLICY
eval(['dynare call_CARA  -DOPTIMAL=1 -DEZ=0 -DDIET=0 ;'])
CRRA_OPT = smm_moment('data_CARA',M_,oo_,options_);


% EZWH LAISSEZ-FAIRE
eval(['dynare call_EZ  -DOPTIMAL=0 -DEZ=1 -DDIET=0 ;'])
EZWH_LF  = smm_moment('data_ez',M_,oo_,options_);
% EZWH OPTIMAL POLICY 
eval(['dynare call_EZ  -DOPTIMAL=1 -DEZ=1 -DDIET=0 ;'])
EZWH_OPT = smm_moment('data_ez',M_,oo_,options_);

% variables to display:
to_show = {'AVG-gy', 'STD-gy', 'STD-gc', 'STD-giT', 'AVG-cyrat',  'AVG-kHkYrat', 'AVG-rFquat', 'STD-rFquat','AVG-realEPquat','STD-realEPquat', 'AVG-mu','AVG-X','AVG-SCC','STD-SCC','AVG-U','STD-U','AVG-gd','AVG-bpann','BETA-ln_SCC-ln_c'};


% show outcome
struct_list = {'EZWH_LF','EZWH_OPT','CRRA_LF','CRRA_OPT'};
fprintf('\n') 
fprintf('TABLE: ESTIMATED PARAMETERS\n')
fprintf('\n') 
fprintf('\t\t\t')
pnames = cell(0);ix = 1;
for iy = 1:2:size(struct_list,2)
    thestruc=eval(struct_list{iy});
    for iyy=1:length(thestruc.param)
        pnames{ix} = thestruc.ve_names{iyy};
        ix = ix+1;
    end
    fprintf(['\t %7s %7s'],struct_list{iy},'(SD)');
end
pnames=unique(pnames);
fprintf('\n')
for ix = 1:size(pnames,2)
    fprintf('%15s',pnames{ix});
    for iy = 1:2:size(struct_list,2)
        thestruc=eval(struct_list{iy});
       
        idx = strmatch(pnames{ix},thestruc.ve_names,'exact');
        if ~isempty(idx)
            %if  output.x_SD
             fprintf('\t %.4f (%.4f)',thestruc.param(idx),thestruc.x_SD(idx));
        else
            fprintf('\t %6s %6s',' - ','(   -  )');
        end
        
    end
    fprintf('\n');
end
fprintf('\n')
fprintf('\n')
fprintf('TABLE: MOMENTS FOR ALTERNATIVE MODELS\n')
fprintf('\n')
fprintf('%15s  \t \t ','Variables');
for iy = 1:size(struct_list,2)
     fprintf(['\t %.8s '],struct_list{iy});
end
fprintf('\n')
readstruc  = @(x,y,z) x.mom(strmatch(y,x.var_names,'exact'),z);
readstrucb = @(x,y,z) x.y_beta(strmatch(y,x.var_names,'exact'),strmatch(z,x.var_names,'exact'));
for ix=1:size(to_show,2)
    rez = split(to_show{ix},'-');
    if strcmp(rez{1},'AVG')
        z = 1;
        strz = 'E';
    elseif strcmp(rez{1},'STD')
        z = 2;
        strz = 'SD';
    elseif strcmp(rez{1},'BETA')
        z = 4;
        strz = 'B';
    else 
        z = 3; strz = 'COR';
    end
    % display name
    if z < 4
        display_name = rez{2};
    else
        display_name = [rez{2} ',' rez{3}];
    end
    fprintf('%15s \t',[strz '(' display_name ')']);

    ICemp = 0;


    if ~ICemp
        fprintf(' \t ','      ','      ','      ');
    else
        is = 2;
        if abs(ICemp(1)) > 100
            is = is -2;
        elseif abs(ICemp(1)) > 10
            is = is -1;
        end
        fprintf([' \t %.' int2str(is) 'f [%.' int2str(is) 'f;%.' int2str(is) 'f]'],ICemp);
    end
    for iy = 1:size(struct_list,2)
        if z < 4
            yshow = readstruc(eval(struct_list{iy}),rez{2},z);
        else
            yshow = readstrucb(eval(struct_list{iy}),rez{3},rez{2});
        end
        % check if estimated
        dummy_estimated_param = 0 ;
        is = 2;
        if isempty(yshow)
            fprintf(['\t %' int2str(is+1) 's  '],' -    ');
        else
            if abs(yshow) > 100
                is = is -2;
            elseif abs(yshow) > 10
                is = is -1;
            end
            if yshow<0
                is=is-1;
            end
            if dummy_estimated_param
                estim_str = '*   ';
            else
                estim_str='    ';
            end
            fprintf(['\t %.' int2str(is) 'f%1s '],yshow,estim_str);
        end
    end
    fprintf('\n',strz,rez{2})
end

if exist('wc') > 0
     fprintf('%15s \t',['E(100ln(' wc{1} '/uss))']);
     fprintf(' \t \t \t \t \t ','      ','      ','      ');
     for iy = 1:size(struct_list,2) 
         is =  2;
        thestruc = eval(struct_list{iy});
        fprintf(['\t %.' int2str(is) 'f%1s '],100*log(thestruc.ln_u),' ');
     end
     fprintf('\n');

     fprintf('%15s \t',['welfcost']);
     fprintf(' \t \t \t \t \t ','      ','      ','      ');
     for iy = 1:size(struct_list,2) 
         is =  2;
        thestruc = eval(struct_list{iy});
        fprintf(['\t %.' int2str(is) 'f%1s '],thestruc.wc,' ');
     end
     fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



