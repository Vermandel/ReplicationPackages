% this function convert variables, parameters and steady states into dynare 
% variables to increase the speed of computation

function constraint1 = read_bind(constraint,M_)

% create a list of delimiters that can separate parameters and endogenoous
% variables in the string that expresses the constraint
delimiters = char(',',';','(',')','+','-','^','*','/',' ','>','<','=');

% split the string that holds the constraint into tokens
tokens = tokenize(constraint,delimiters);

ntokens = length(tokens);

rnd_dec = 5; % number of decimals when rounding

% search for tokens that match the list of endogenous variables
for i=1:ntokens
	token = deblank(tokens{i});

    if size(token,1) == 1 &&  strcmp(token,'<')
		tokens{i} = ['- 1.0e-' num2str(rnd_dec) ' <'];
    end
    if size(token,1) == 1 &&  strcmp(token,'>')
		tokens{i} = ['+ 1.0e-' num2str(rnd_dec) ' >'];
    end
	% read token
	if size(token,2) > 3 &&  strcmp(token(end-2:end),'_ss') 
		% then it's a ss
		% look for ss
		idx = strmatch(token(1:end-3),M_.endo_names,'exact');
		tokens{i} = ['round(pw.oo_.dr.ys(' num2str(idx) '),8)'];
	else
		idv = strmatch(token,M_.endo_names,'exact');
		idp = strmatch(token,M_.param_names,'exact');
		if ~isempty(idv)
			% its a endogenous variable
			%tokens{i} = ['round((ynew(pw.oo_.dr.inv_order_var(' num2str(idv) '))+pw.oo_.dr.ys(' num2str(idv) ')),' num2str(rnd_dec) ')-1e-' num2str(rnd_dec)  ];
			tokens{i} = ['(ynew(pw.oo_.dr.inv_order_var(' num2str(idv) '))+pw.oo_.dr.ys(' num2str(idv) '))'  ];
		elseif ~isempty(idp)
			% its a parameter
			tokens{i} = ['pw.M_.params(' num2str(idp) ')'];		
		end
	end
end

% reassemble the tokens to create a string that expresses the constraint
constraint1 = strmerge(tokens);