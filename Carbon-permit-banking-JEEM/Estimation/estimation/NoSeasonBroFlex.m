function [dt] = NoSeasonBroFlex(Y,time_scale)


global options_
global oo_
	
	if nargin < 2
		time_scale = 4; % 12 = Monthly
	end
	
	
	%Y = mydef;
	
	Y_bk = Y;
	% no nan
	noNaN = find(~isnan(Y));
	
	Y = Y(noNaN);
	T = length(Y);
		
	
	%%%%%%%%%%%%%% Step 2. Detrend the data using a 13-term moving average
	sW13 = [1/24;repmat(1/12,time_scale-1,1);1/24]*12/time_scale;
	
	Ys = conv(Y,sW13,'same');
	Ys(1:(time_scale/2)) = Ys((time_scale/2+1));
	Ys((T-time_scale/2+1):T) = Ys(T-time_scale/2);


	xt = Y./Ys;

	%h = plot(date,Y,'b',date,Ys,'r','LineWidth',2);
	%legend(h,'13-Term Moving Average')
	%hold off




	%%%%%%%%%%%%%% Step 3. Create seasonal indices
	s = time_scale;
	for i = 1:s
	 sidx{i,1} = i:s:T;
	end

	%%%%%%%%%%%%%% Step 4. Apply an S3×3 filter
	% S3x3 seasonal filter
	% Symmetric weights
	sW3 = [1/9;2/9;1/3;2/9;1/9];
	% Asymmetric weights for end of series
	aW3 = [.259 .407;.37 .407;.259 .185;.111 0];

	% Apply filter to each month
	shat = NaN*Y;
	for i = 1:s
		Ns = length(sidx{i});
		first = 1:4;
		last = Ns-3:Ns;
		dat = xt(sidx{i});
		
		sd = conv(dat,sW3,'same');
		sd(1:2) = conv2(dat(first),1,rot90(aW3,2),'valid');
		sd(Ns-1:Ns) = conv2(dat(last),1,aW3,'valid');
		shat(sidx{i}) = sd;
	end

	% 13-term moving average of filtered series
	%sW13 = [1/24;repmat(1/12,11,1);1/24];
	sW13 = [1/24;repmat(1/12,time_scale-1,1);1/24]*12/time_scale;
	sb = conv(shat,sW13,'same');
	%sb(1:6) = sb(s+1:s+6); 
	sb(1:(time_scale/2)) = sb(s+1:s+time_scale/2); 
	%sb(T-5:T) = sb(T-s-5:T-s);
	sb((T-time_scale/2+1):T) = sb((T-time_scale/2+1-s):T-s);
	
	% Center to get final estimate
	s33 = shat./sb;

	%figure
	%plot(date,s33)


	%%%%%%%%%%%%% Step 5. Apply a 13-term Henderson filter
	% Deseasonalize series
	dt = Y./s33;
	
	if time_scale == 12
		% Henderson filter weights
		Sym   = [ -0.01935	-0.02787	0	0.06549	0.14736	0.21434	0.24006	0.21434	0.14736	0.06549	0	-0.02787	-0.01935 ]';
		Asym6 =	[ 0	-0.01643	-0.02577	0.00127	0.06594	0.14698	0.21314	0.23803	0.21149	0.14368	0.06099	-0.00533	-0.03401]';
		Asym5 =	[ 0	0	-0.011	-0.02204	0.0033	0.06626	0.14559	0.21004	0.23324	0.20498	0.13547	0.05108	-0.01695]';
		Asym4 =	[ 0	0	0	-0.00814	-0.0202	0.00413	0.06608	0.14441	0.20784	0.23002	0.20076	0.13024	0.04483]';
		Asym3 =	[ 0	0	0	0	-0.01604	-0.02487	0.00267	0.06784	0.14939	0.21605	0.24144	0.2154	0.1481]';
		Asym2 =	[ 0	0	0	0	0	-0.04271	-0.03864	0.00182	0.0799	0.17436	0.25392	0.29223	0.2791]';
		Asym1 =	[ 0	0	0	0	0	0	-0.09187	-0.05812	0.01202	0.11977	0.2439	0.35315	0.4211]';
		% Henderson filter weights
		%sWH = [-0.019;-0.028;0;.066;.147;.214;.24;.214;.147;.066;0;-0.028;-0.019];
		sWH = Sym;

		aWH = [Asym6(end:-1:2) Asym5(end:-1:2) Asym4(end:-1:2) Asym3(end:-1:2) Asym2(end:-1:2) Asym1(end:-1:2)];

		% Asymmetric weights for end of series
		%aWH = [-.034  -.017   .045   .148   .279   .421;
		%	   -.005   .051   .130   .215   .292   .353;
		%		.061   .135   .201   .241   .254   .244;
		%		.144   .205   .230   .216   .174   .120;
		%		.211   .233   .208   .149   .080   .012;
		%		.238   .210   .144   .068   .002  -.058;
		%		.213   .146   .066   .003  -.039  -.092;
		%		.147   .066   .004  -.025  -.042  0    ;
		%		.066   .003  -.020  -.016  0      0    ;
		%		.001  -.022  -.008  0      0      0    ;
		%	   -.026  -.011   0     0      0      0    ;
		%	   -.016   0      0     0      0      0    ];
	end
	if time_scale == 4
		% weights
		sWH		= [ -0.07343	0.29371	0.55944	0.29371	-0.07343 ]';
		Asym2	= [ 0	-0.04419	0.29121	0.52522	0.22776 ]';
		Asym1	= [ 0	0	-0.13181	0.36713	0.76467 ]';

		aWH = [Asym2(end:-1:2) Asym1(end:-1:2)];
	end
	
	
	
	
	
	
	% Apply 13-term Henderson filter
	first = 1:time_scale;
	last = T-(time_scale-1):T;
	h13 = conv(dt,sWH,'same');
	h13(T-time_scale/2+1:end) = conv2(dt(last),1,aWH,'valid');
	h13(1:time_scale/2) = conv2(dt(first),1,rot90(aWH,2),'valid');

	% New detrended series
	xt = Y./h13;


	%%%%%%%%%%%%%% Step 6. Apply an S3×5 seasonal filter.
	% S3x5 seasonal filter 
	% Symmetric weights
	sW5 = [1/15;2/15;repmat(1/5,3,1);2/15;1/15];
	% Asymmetric weights for end of series
	aW5 = [.150 .250 .293;
		   .217 .250 .283;
		   .217 .250 .283;
		   .217 .183 .150;
		   .133 .067    0;
		   .067   0     0];

	% Apply filter to each month
	shat = NaN*Y;
	for i = 1:s
		Ns = length(sidx{i});
		first = 1:6;
		last = Ns-5:Ns;
		dat = xt(sidx{i});
		
		sd = conv(dat,sW5,'same');
		sd(1:3) = conv2(dat(first),1,rot90(aW5,2),'valid');
		sd(Ns-2:Ns) = conv2(dat(last),1,aW5,'valid');
		shat(sidx{i}) = sd;
	end

	% 13-term moving average of filtered series
	sW13 = [.5/time_scale;repmat(1/time_scale,time_scale-1,1);.5/time_scale];
	sb = conv(shat,sW13,'same');
	%sb(1:6) = sb(s+1:s+6); 
	sb(1:(time_scale/2)) = sb(s+1:s+time_scale/2); 
	%sb(T-5:T) = sb(T-s-5:T-s);
	sb((T-time_scale/2+1):T) = sb((T-time_scale/2+1-s):T-s);


	% Center to get final estimate
	s35 = shat./sb;

	% Deseasonalized series
	dt1 = Y./s35;
	% add the NaN
	Y = Y_bk;
	Y(noNaN) = dt1;
	dt = Y;

	plot(1:length(Y),dt,1:length(Y),Y_bk)
