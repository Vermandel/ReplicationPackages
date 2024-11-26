function [vv,P,Q,R] = update_dr(vv,q,tbool,pw)
%UPDATE_DR 
% create state dependent decision rules given a duration q
    Rg=pw.Nr;
	% Compute state-dependent DR from duration q
 	
    % tip for the loop:
    % i1 = 0 --> no constraint binds anymore
	% i1 > 0 --> binds
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
    
    vname = sprintf(['.' repmat('%i',[1 Rg])],flip(tbool(:,(end+1-q):end),2));
    idx = strmatch(vname(2:end),vv.tree,'exact');
    P=vv.PP(:,:,idx);
    Q=vv.QQ(:,:,idx);
    R=vv.RR(:,:,idx);
end

