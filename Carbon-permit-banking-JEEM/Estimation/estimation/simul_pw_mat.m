function [y_pw,y_lin,rbool] = simul_pw_mat(pw,zz)
%SIMUL_OCB Summary of this function goes here
%   Detailed explanation goes here

Rg = pw.Nr;
M_ = pw.M_;
oo_ = pw.oo_;
Tsample=size(zz,2);

if size(zz,1)~=M_.exo_nbr
	error('Wrong number of shocks');
end

rbool   = zeros(Rg,Tsample+1);
sbool   = '[';
for i1 = 1:Rg
    if i1 > 1
        sbool   = [sbool ';'];
    end
    sbool   = [sbool eval(['read_bind(pw.c' num2str(i1) '_bind,M_)']) ];
end
sbool   = [sbool ']'];

qMAX = 500; % max number of forward periods for guess and try

y_lin = zeros(M_.endo_nbr,Tsample+1);	% state variables
y_pw  = zeros(M_.endo_nbr,Tsample+1);	% state variables



vv = [];
% initialize with the normal regime:
reg_name=sprintf([repmat('%i',[1 Rg])],pw.combinations(1,:));
vv.tree=char(reg_name);
idx     = 1;
vv.SS(:,:,idx) =   eval(['-inv(pw.FF' reg_name '*pw.PP+pw.GG' reg_name ')']);
vv.PP(:,:,idx) =   eval(['vv.SS(:,:,idx)*pw.HH' reg_name ]);
vv.QQ(:,:,idx) =   eval(['vv.SS(:,:,idx)*pw.MM' reg_name ]);
vv.RR(:,:,idx) =   eval(['vv.SS(:,:,idx)*(pw.CC' reg_name '+pw.FF' reg_name '*zeros(pw.M_.endo_nbr,1))']);

    
for t = 2:(Tsample+1)
    
	% linear case
	y_lin(:,t)	= pw.PP*y_lin(:,t-1)+pw.QQ*zz(:,t-1);

	% piece wise
	ynew 		= pw.PP*y_pw(:,t-1)+pw.QQ*zz(:,t-1);
	rbool(:,t)  = eval(sbool);
	
	if sum(rbool(:,t))>0 % then some regimes bind
		%save pw_mat; error('STOP');
 				
		%% step 1 :
		% get linear duration (inital guess)
		b = rbool(:,t);
		ib = 0; % counter
		tbool = rbool(:,t);
		while sum(b) > 0
				ib = ib + 1;
				ynew=pw.PP*(ynew);
				tbool = [ tbool   eval(sbool) ];
				b =  sum(tbool(:,end));
		end
		q=ib;
		tbool = tbool(:,1:ib); % take away last empty row
		
		
		
		% step 2: guess and try
		violation = 1;
		iter = 0;
		while violation == 1
            iter=iter+1;
			if iter>1000
				save gg;
				error('lol')
			end
%			[t q]
			% i1 = 0 --> no constraint binds anymore
			% i1 > 0 --> binds
			% creates matrices
			for i1 = 1:q
				% name of the history tree based on P(t+1)
				vname = sprintf(['.' repmat('%i',[1 Rg])],flip(tbool(:,(end+1-i1):end),2));
				vname = vname(2:end);
				
				idx = strmatch(vname,vv.tree,'exact');
				
				% Compute current matrices
				if isempty(idx) 
						
					% initiate the future matrix of P(+1) and R(+1)
					if i1==1
						% next period is normal times
						idnext = 1;
					else
						% find id within the tree
						prevname = sprintf(['.' repmat('%i',[1 Rg])],flip(tbool(:,(end+2-i1):end),2));
						idnext = strmatch(prevname(2:end),vv.tree,'exact');
					end
					% get next matrix PP and RR
					PPnext = vv.PP(:,:,idnext);
					RRnext = vv.RR(:,:,idnext);
					
					% name of the current system to be used (ie F, G, H, M, C)
					reg_name=sprintf([repmat('%i',[1 Rg])],tbool(:,end+1-i1));
					% load current system matrices
					FF	    = 	eval(['pw.FF' reg_name ]);
					GG	    = 	eval(['pw.GG' reg_name ]);
					HH	    = 	eval(['pw.HH' reg_name ]);
					MM	    = 	eval(['pw.MM' reg_name ]);
					CC	    = 	eval(['pw.CC' reg_name ]);

					% matrices not yet created:
					% S(t) = -inv(F*P(t+1)+G);
					% P(t) = S*H;
					% Q(t) = S*M;
					% R(t) = S*(C+F*R(t+1));
					
					% get new row id
					ndim = size(vv.SS,3)+1;
					vv.tree = char(vv.tree,vname);
					% compute and save new dr
					vv.SS(:,:,ndim) =  -inv(FF*PPnext+GG);
					vv.PP(:,:,ndim) =   vv.SS(:,:,ndim)*HH;
					vv.QQ(:,:,ndim) =   vv.SS(:,:,ndim)*MM;
					vv.RR(:,:,ndim) =   vv.SS(:,:,ndim)*(CC + FF*RRnext );
					
				end
							
			end % for loop ends
			% from here we have the value of P(t), Q(t) and R(t)
			
			if q == 0
                error('errrrr')
            end
			% now let's see if PQ is consistent
			% keep prev history and use state dependent policy rule
			y_test = [y_pw(:,t-1) zeros(M_.endo_nbr,q)];
			for i1 = 1:q
				% update name (in reverse order here!)
				% we use the longest PQ
				vname = sprintf(['.' repmat('%i',[1 Rg])],flip(tbool(:,i1:end),2));
				idx = strmatch(vname(2:end),vv.tree,'exact');
				% propagation PP(t)*x(-1)
				ynew = vv.PP(:,:,idx)*y_test(:,i1);
				% constant term RR(t)
				ynew = ynew + vv.RR(:,:,idx);
				if i1 ==1 % then add the current shock
					ynew = ynew + vv.QQ(:,:,idx)*zz(:,t-1);
				end
				%save result
				y_test(:,i1+1) = ynew;
				% check whether the path is consistent with guess
				if i1 > 1 && sum(eval(sbool)-tbool(:,i1))~=0
					%booldiff = eval(sbool)-tbool(:,i1);
					
					%error('err..');
					% if not consistent
					% we are not sure of the outcome in the next periods
					% so delete what's next
					q = i1-1;
					y_test = y_test(:,1:(i1));
					tbool = tbool(:,1:(i1-1));
					%[t q q]
					% stop current loop here
					% the path is good until now
					break; 
				end
			end % for loop ends
			
			% check whether condition is violated at next period
			ynew = pw.PP*y_test(:,end);
			tempbool = eval(sbool);
			if 1 == 0
				for i1 = 1:M_.endo_nbr
					fprintf('%5s = %f \n', M_.endo_names(oo_.dr.order_var(i1),:),ynew(i1));
				end
			end

			if sum(tempbool) == 0 || q==qMAX
				% termination condition met
				%q
				% leave the loop plz
				violation = 0;
			else
				
				% not finished
				q = q +1;
				% update the booleans
				tbool = [tbool tempbool];
				
				% if still binding now
				% check if we can improve guess
				if sum(tempbool) > 0
					% guess and try again to find new normal regime
					while sum(tempbool) ~= 0
						% check whether condition is violated at next period
						ynew = pw.PP*ynew;
						tempbool = eval(sbool);
						% update the booleans
						if sum(tempbool) > 0
							tbool = [tbool tempbool];
							q=q+1;
						end
					end
				end
			end
			
        end

		% get last iteration
		ynew=y_test(:,2);

		% save new values
		y_pw(:,t) = ynew;
		
	else	
		y_pw(:,t) = ynew;
	end	

end
% add steady state
y_lin = y_lin(oo_.dr.inv_order_var,:)' + repmat(pw.oo_.dr.ys',[Tsample+1 1]);
y_pw = y_pw(oo_.dr.inv_order_var,:)' + repmat(pw.oo_.dr.ys',[Tsample+1 1]);

end

