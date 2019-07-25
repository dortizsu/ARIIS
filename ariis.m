function A = ariis(m,inputs,constants)
	% A = ariis(m,inputs,constants)
	%
	% version: 2.2, 7/2/2019
	%
	% ... preamble
	% Algorithm for Robustly Identifying the Inertial Subrange (ARIIS)
	% Based on Kolmogorov's 1941 hypotheses, ARIIS finds the most plausible bandwidth for
	% the inertial subrange in inputted velocity data OR computed autovariance spectra.
	% ARIIS is built around 2 major sub-routines:
	% 1) Isotropic bandwidth convergence
	% 2) Robust log-scaled linear regression
	% (1) Finds the inertial subrange bandwidth in frequency space. (2) determines the 
	% power law for each flow component over the identified subrange.
	%
	% Arguments:
	% m --> mode of operation, 0 if giving velocities, 1 if giving spectra.
	% inputs --> 
	% m = 0)
	% inputs.u = instantaneous stream-wise velocity [m/s]
	% inputs.v = instantaneous cross-stream velocity [m/s]
	% inputs.w = instantaneous vertial velocity [m/s]
	% constants.Uadv = mean, advection velocity [m/s], !!required!!
	% constants.z = height of measurement in surface layer [m], !!required!!
	% constants.dt = sampling interval [s], !!required!!
	% constants.smoothing = spectral smoothing yes (1) or no (0), optional (default = 1)
	% constants.nfa = level of smoothing, optional (default = 8)
	% constants.fcutoff = high-frequency cut-off for spectra, optional (default = Nyquist)
	% constants.Suwratio = Cut-off for isotropy convergence, optional (default = 4/3)
	% constants.Rfcutoff = low frequency limit for inertial subrange, optional (default = U/z)
	% constants.a = fraction of identified inner bandwidth to keep, optional (default 0.9)
	% constants.Dthresh = Cook's Distance threshhold numerator, option (default 4)
	%.......
	% m = 1)
	% NO spectral filtering or smoothing is done, this must all be done a priori
	% inputs.n = frequency bins [Hz]
	% inputs.Suu = stream-wise velocity autovariance spectrum [(m^2/s^2)/Hz]
	% inputs.Svv = cross-stream velocity autovariance spectrum [(m^2/s^2)/Hz]
	% inputs.Sww = vertical velocity autovariance spectrum [(m^2/s^2)/Hz]
	% constants.ust = friction velocity [m/s], !!required!!
	% constants.varw = vertical velocity variance [m^2/s^2], !!required!!
	% (the rest are the same as mode 0)
	%
	% Output: A, 1 x 24 row vector
    % 1 --> number of frequency bins in the identified subrange
	% 2 --> [Hz] low-frequency limit of inertial subrange
	% 3 --> [rad/m] low-wavenumber limit of inertial subrange
	% 4 --> [Hz] high-frequency limit of inertial subrange
	% 5 --> [rad/m] high-wavenumber limit of inertial subrange
	% 6 --> u-w Isotropy coefficient
	% 7 --> u-w Isotropy coefficient, only from f-bins used in fitting
	% 8 --> u-v Isotropy coefficient
	% 9 --> u-v Isotropy coefficient, only from f-bins used in fitting
	% 10--> power (p) of Suu ~ A*f^{p}
	% 11--> linear coefficient of determination from fit (i.e Pearson's correlation)
	% 12--> delta for p, as in p +/- delta spans 95% confidence interval
	% 13--> power (p) of Svv ~ A*f^{p}
	% 14--> linear coefficient of determination from fit (i.e Pearson's correlation)
	% 15--> delta for p, as in p +/- delta spans 95% confidence interval
	% 16--> power (p) of Sww ~ A*f^{p}
	% 17--> linear coefficient of determination from fit (i.e Pearson's correlation)
	% 18--> delta for p, as in p +/- delta spans 95% confidence interval
	%
	flags = [];
	% References:
	% Kolmogorov, A. The Local Structure of Turbulence in Incompressible Viscous Fluid for Very Large Reynolds Numbers. Dokl. Akad. Nauk SSSR 1941.
	% Kolmogorov, A. Dissipation of Energy in Locally Isotropic Turbulence. Comptes Rendus l’Acadamie des Sci. al’U.R.S.S. 1941, 32, 16.
	% Jimenez, J.; Wray, A.A.; Saffman, P.G.; Rogallo, R.S. The structure of intense vorticity in homogeneous isotropic turbulence. Cent. Turbul. Res. Proc. Summer Progr. 1992, p. 24.
	% Chatterjee, S.; Hadi, A.S. Influential Observations, High Leverage Points, and Outliers in Linear Regression, 1988. doi:10.2307/2245477.
	% http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals_print.html
	%
	% Disclaimer: 
	% author: D.G. Ortiz-Suslow
	% revised: 01/31/19
	% disclaimer: This software is provided as-is with no guarantee of functionality or fitness for the task for which it was designed to complete.
	%
	% Updates:
	% v1.1 --> changed default low frequency cut-off for R convergence to U/z, instead of 1/3 Hz; more general.
	% v2 --> changed to fitting non-scaled spectra (expected p = -5/3) and added 95% confidence bounds for p
    % v2.1 --> include number of frequency bins in subrange as output
	% v2.2 --> fixed bug in 95% confidence bound calculation, requires modifying local copy of rlogfit; adjust updated output, to only give delta (not full error bounds) and to NOT give amplitude
	%..................................................................................
	%...Prepare ARIIS..................................................................
	%..................................................................................
	%..pre-define defaults
	defaulted_fields = {...
		'smoothing',...		
		'nfa',...		
		'fcutoff',...
		'Suwratio',...
		'Rfcutoff',...
		'a',...
		'Dthresh'};
	defaulted_vals = {...
		'1',...	
		'8',...	
		'(1/dt)/2',...
		'4/3',...
		'Uadv/z',... % largest wall-bounded eddy that can be resolved has frequency n = Uadv/z
		'0.9',...
		'4'};
	Ndefaults = length(defaulted_fields);
	fittype ='loglog';
	spectrafits = {'Fuu' 'Fvv' 'Fww'};
	%..ARIIS flags set to 0
	no_isotropic_convergence = 0;
	subrange_too_short = 0;
	subrange_too_noisy = [0 0 0];
	flags = [no_isotropic_convergence subrange_too_short subrange_too_noisy];
	if m == 0
		%..ARIIS will calculate autovariance spectra for you.
		%..Implications: you trust ARIIS spectral smoothing
		%
		%...Grab input data
		u = inputs.u;
		v = inputs.v;
		w = inputs.w;
		%..get constants needed for ARIIS, most have default values
		%..non-defaulted values
		Uadv = constants.Uadv;
		z = constants.z;
		dt = constants.dt;
		%..values with ARIIS defaults
		for ii = 1:Ndefaults
			if isfield('constants',['constants.' defaulted_fields{ii}])
				evalc([defaulted_fields{ii} ' = constants.' defaulted_fields{ii}]);
            else
				evalc([defaulted_fields{ii} ' = ' defaulted_vals{ii}]);
			end
		end
		%..................................................................................
		%...Step 1: Initialization.........................................................
		%..................................................................................
		%..define some scaling values
		ust = ustar(u,v,w);
		ust2 = ust^2;
		varw = var(w);
		nscale = z/Uadv;
		%..calculate autovariance spectra
		Cuw = spectf(u,w,dt,1,1);
		Cuv = spectf(u,v,dt,1,1);
		n = Cuw(:,1);
		Suu = Cuw(:,2);
		Sww = Cuw(:,3);
		Svv = Cuv(:,3);
		%..if needed, can define a high frequency cut-off
		Suu(n>fcutoff) = [];
		Sww(n>fcutoff) = [];
		Svv(n>fcutoff) = [];
		n(n>fcutoff) = [];
		%..if needed, can smooth spectra using a log-uniform bin mean averager
		if true(smoothing)
			% smooth
			[smoothed,~] = logSmooth([n Suu Svv Sww],nfa);
			n = smoothed(:,1);
			Suu = smoothed(:,2);
			Svv = smoothed(:,3);
			Sww = smoothed(:,4);
        end
		%..ARIIS guts below................................................................
		%..................................................................................
		%...Step 2: Isotropic Convergence Sub-Routine......................................
		%..................................................................................
		R = Suu./Sww;
		iR = 1;
		ix = 0;
		%..find low frequency limit of isotropic subrange
		while ix < 1 && iR <= length(R)
			ix = R(iR) < Suwratio;
			if ix == 1 && n(iR) <= Rfcutoff
				ix = 0;
            end
			iR = iR + 1;
        end
		%..find high-frequency limit of isotropic subrange
		if iR < length(R)
			ix = find(R(iR:end) < Suwratio,1,'last');
			if isempty(ix)
				iR(2) = iR;
			else
				iR(2) = ix + iR-1;
			end
			iR(1) = iR(1)-1;
			IX = iR(1):iR(2);
			%..if needed, define what fraction of final inertial subrange is kept
			trim = round((1 - a)*length(IX));
			IX = IX(trim+1:end-trim);
		elseif iR >= length(R)
			% warning('No isotropic convergence subrange.')
			flags(:,1) = 1;
		end
		if length(IX) < 4
			% warning('Subrange too narrow.')
			flags(:,2) = 1;
        end
		%..................................................................................
		%...Step 3: Robust fitting sub-routine.............................................
		%..................................................................................
		if any(flags)
			warning('Flag tripped in isotropic convergence routine. No output.')
			A = [];
			return
        else
			%..Isotropic coefficient, Jimenez et al. 1993
			k = 2*pi*n(IX)/Uadv;
			dEuudk = rate(Suu(IX))./rate(k);
			Iuw = (Suu(IX) - k.*dEuudk)./(2*Sww(IX));
			Iuv = (Suu(IX) - k.*dEuudk)./(2*Svv(IX));
			%..Iterative fitting method for finding power law over Sxx(IX)
			%..scale spectra in surface layer coordinates
			f = n(IX); %n(IX)*nscale;
			Fuu = Suu(IX); %Fuu = n(IX).*Suu(IX)/ust2;
			Fvv = Svv(IX); %n(IX).*Svv(IX)/varw;
			Fww = Sww(IX); %n(IX).*Sww(IX)/varw;
			%..robust iterative fitting
			for ii = 1:length(spectrafits)
				evalc(['[robfit(' num2str(ii) '),rawfit(' num2str(ii) ')]' '= rlogfit(f,' spectrafits{ii} ',fittype,Dthresh)']);
				if sum(robfit(ii).indx) < 3
					flags(:,ii+2) = 1;
				end
			end
			I = [median(Iuw) median(Iuw(robfit(1).indx)) median(Iuv) median(Iuv(robfit(1).indx))];
			%.................................................................................
			%...Step 4: Output................................................................
			%.................................................................................
			if any(flags)
				warning('Too much spectral noise for confident slope determination. No output.')
				A = [];
			else
				A = [...
                numel(IX),...%....................1
				n(IX(1)),...%.....................2
				k(1),...%.........................3
				n(IX(end)),...%...................4
				k(end),...%.......................5
				I,...%............................6/7/8/9
				robfit(1).coeffs(2),...%..........10
				robfit(1).r2,...%.................11
				range(robfit(1).ci(:,2))/2,...%...12
				robfit(2).coeffs(2),...%..........13
				robfit(2).r2,...%.................14
				range(robfit(2).ci(:,2))/2,...%...15
				robfit(3).coeffs(2),...%..........16
				robfit(3).r2,...%.................17
				range(robfit(3).ci(:,2))/2];%.....18
			end
		end
	elseif m == 1
		%..ARIIS will use the spectra you give it.
		%..Implications: you don't want ARIIS spectral smoothing
		%
		%...Grab input data
		n = inputs.n;
		Suu = inputs.Suu;
		Svv = inputs.Svv;
		Sww = inputs.Sww;
		%..get constants needed to ARIIS, most have default values
		%..non-defaulted values
		ust = constants.ust;
		varw = constants.varw;
		Uadv = constants.Uadv;
		z = constants.z;
		dt = constants.dt;
		%..values with ARIIS defaults
		for ii = 1:Ndefaults
			if isfield('constants',['constants.' defaulted_fields{ii}])
				evalc([defaulted_fields{ii} ' = constants.' defaulted_fields{ii}]);
            else
				evalc([defaulted_fields{ii} ' = ' defaulted_vals{ii}]);
			end
		end
		%..................................................................................
		%...Step 1: Initialization.........................................................
		%..................................................................................
		%..define some scaling values
		ust2 = ust^2;
		nscale = z/Uadv;
		%..ARIIS guts below................................................................
		%..................................................................................
		%...Step 2: Isotropic Convergence Sub-Routine......................................
		%..................................................................................
		R = Suu./Sww;
		iR = 1;
		ix = 0;
		%..find low frequency limit of isotropic subrange
		while ix < 1 && iR <= length(R)
			ix = R(iR) < Suwratio;
			if ix == 1 && n(iR) <= Rfcutoff
				ix = 0;
            end
			iR = iR + 1;
        end
		%..find high-frequency limit of isotropic subrange
		if iR < length(R)
			ix = find(R(iR:end) < Suwratio,1,'last');
			if isempty(ix)
				iR(2) = iR;
			else
				iR(2) = ix + iR-1;
			end
			iR(1) = iR(1)-1;
			IX = iR(1):iR(2);
			%..if needed, define what fraction of final inertial subrange is kept
			trim = round((1 - a)*length(IX));
			IX = IX(trim+1:end-trim);
		elseif iR >= length(R)
			% warning('No isotropic convergence subrange.')
			flags(:,1) = 1;
		end
		if length(IX) < 4
			% warning('Subrange too narrow.')
			flags(:,2) = 1;
        end
		%..................................................................................
		%...Step 3: Robust fitting sub-routine.............................................
		%..................................................................................
		if any(flags)
			warning('Flag tripped in isotropic convergence routine. No output.')
			A = [];
			return
        else
			%..Isotropic coefficient, Jimenez et al. 1993
			k = 2*pi*n(IX)/Uadv;
			dEuudk = rate(Suu(IX))./rate(k);
			Iuw = (Suu(IX) - k.*dEuudk)./(2*Sww(IX));
			Iuv = (Suu(IX) - k.*dEuudk)./(2*Svv(IX));
			%..Iterative fitting method for finding power law over Sxx(IX)
			%..(old)scale spectra in surface layer coordinates
			f = n(IX); %n(IX)*nscale;
			Fuu = Suu(IX); %Fuu = n(IX).*Suu(IX)/ust2;
			Fvv = Svv(IX); %n(IX).*Svv(IX)/varw;
			Fww = Sww(IX); %n(IX).*Sww(IX)/varw;
			%..robust iterative fitting
			for ii = 1:length(spectrafits)
				evalc(['[robfit(' num2str(ii) '),rawfit(' num2str(ii) ')]' '= rlogfit(f,' spectrafits{ii} ',fittype,Dthresh)']);
				if sum(robfit(ii).indx) < 3
					flags(:,ii+2) = 1;
				end
			end
			I = [median(Iuw) median(Iuw(robfit(1).indx)) median(Iuv) median(Iuv(robfit(1).indx))];
			%.................................................................................
			%...Step 4: Output................................................................
			%.................................................................................
			if any(flags)
				warning('Too much spectral noise for confident slope determination. No output.')
				A = [];
			else
				A = [...
                numel(IX),...%....................1
				n(IX(1)),...%.....................2
				k(1),...%.........................3
				n(IX(end)),...%...................4
				k(end),...%.......................5
				I,...%............................6/7/8/9
				robfit(1).coeffs(2),...%..........10
				robfit(1).r2,...%.................11
				range(robfit(1).ci(:,2))/2,...%...12
				robfit(2).coeffs(2),...%..........13
				robfit(2).r2,...%.................14
				range(robfit(2).ci(:,2))/2,...%...15
				robfit(3).coeffs(2),...%..........16
				robfit(3).r2,...%.................17
				range(robfit(3).ci(:,2))/2];%.....18
			end
		end
	else
		error('unknown (m)ode...')
	end
end
%
%.................................................................................
%...APPENDIX......................................................................
%.................................................................................
function ust = ustar(u,v,w)
    % get friction velocity from 3 velocity components
	%
	% author: D.G. Ortiz-Suslow
	% revised: 01/31/19
	% disclaimer: This software is provided as-is with no guarantee of functionality or fitness for the task for which it was designed to complete.
    u = detrend(demean(u));
    v = detrend(demean(v));
    w = detrend(demean(w));
    ii = mean(-u.*w);
    jj = mean(-v.*w);
	ust = sqrt(sqrt(ii^2 + jj^2));
end
%
function rated = rate(x,flag)
	% rated = rate(x,flag)
	% flag = 1 (output rate vector), flag = 2 (output mean rate), flag = 3 (output [mean std])
	% calculate forward difference of a column vector. Based on flag output is a column vector of same dimension as input, the mean difference, or the mean and stan. dev.
	%
	% author: D.G. Ortiz-Suslow
	% revised: 01/31/19
	% disclaimer: This software is provided as-is with no guarantee of functionality or fitness for the task for which it was designed to complete.
	%
	if ~iscolumn(x)
	    x = x';
	end
	if nargin < 2
	    rated = [diff(x(1:2));diff(x(:))];
	elseif flag == 1
	    rated = mean(diff(x));
	elseif flag == 2
	    rated = [mean(diff(x)) std(diff(x))];
	else
	    disp('you never see this')
	end
end
%
function [rfitted,raw] = rlogfit(x,y,varargin)
	% [rfitted,raw] = rlogfit(x,y)
	% [rfitted,raw] = rlogfit(x,y,Dthresh)
	% [rfitted,raw] = rlogfit(x,y,fittype,Dthresh)
	%
	% Iterative, robust linear least-squares regression.
	% Method uses a Cook's Distance algorithm to flag/remove suspicious and influential points from final regression model.
	% function is based on logfit.m (Matlab fileexchange)
	%
	% author: D.G. Ortiz-Suslow
	% revised: 01/31/19
	% disclaimer: This software is provided as-is with no guarantee of functionality or fitness for the task for which it was designed to complete.
	%
	% Initial regression model fitting
	if nargin < 4
		[fittype,M,B,MSE,R2,S] = logfit(x,y,'nograph');
		Dthresh = 4; % default value
	elseif nargin == 4
		fittype = varargin{1};
		Dthresh = varargin{2};
		[M,B,MSE,R2,S] = logfit(x,y,fittype,'nograph');
	else
		error('Check number of inputs....')
	end
	if nargin == 3
		Dthresh = varargin{1};
	end
	% initial values
	N = length(x);
	p = 1; % set for completeness, but rlogfit cannot handle more than 1 predictor, y
	Dthresh = Dthresh/(N - p - 1); % Chatterjee and Hadi 1988
	alpha = 0.05; % 95% percentile for two-sided t statistic used in defining confidence interval
	ci = tinv(1-alpha,S.df-1);
	%
	% calculate fitted response variable & complete residual
	if strcmp('linear',fittype)
		Yfit = B + M*x;
		SE = stderr(x,y,Yfit);
		coeffs = [B M];
		bnds = [...
			B-SE(2)*ci M-SE(1)*ci;...
			B+SE(2)*ci M+SE(1)*ci];
	elseif strcmp('logx',fittype)
		Yfit = B + M*log10(x);
		SE = stderr(log10(x),y,Yfit);
		coeffs = [B M];
		bnds = [...
			B-SE(2)*ci M-SE(1)*ci;...
			B+SE(2)*ci M+SE(1)*ci];
	elseif strcmp('logy',fittype)
		Yfit = (10^B)*(10^M).^x;
		SE = stderr(x,log10(y),log10(Yfit));
		coeffs = [10^B 10^M];
		bnds = [...
			10^(B-SE(2)*ci) 10^(M-SE(1)*ci);...
			10^(B+SE(2)*ci) 10^(M+SE(1)*ci)];
	elseif strcmp('loglog',fittype)
		Yfit = (10^B)*x.^M;
		SE = stderr(log10(x),log10(y),log10(Yfit));
		coeffs = [10^B M];
		bnds = [...
			10^(B-SE(2)*ci) 10^(M-SE(1)*ci);...
			10^(B+SE(2)*ci) 10^(M+SE(1)*ci)];
	end
	% output to structure
	raw.model = fittype;
	raw.coeffs = coeffs; % [B M], re-scaled
	raw.ci = bnds; % col 1 = lower (upper) bound for B, col 2 = lower (upper) bound for M | re-scaled
	raw.yfit = Yfit;	
	raw.mse = MSE;
	raw.r2 = R2;
	%
	%
	% begin Cook's distance algorithm
	Yr = Yfit - y; % raw residual
	D = zeros(size(x));
	for j = 1:N
		nfit = 1:N;
		nfit(j) = [];
		% regression of x and y excluding jth sample
		[m,b,~,~,~] = logfit(x(nfit),y(nfit),fittype,'nograph');
		if strcmp('linear',fittype)
			yfit = b + m*x;
		elseif strcmp('logx',fittype)
			yfit = b + m*log10(x);
		elseif strcmp('logy',fittype)
			yfit = (10^b)*(10^m).^x;
		elseif strcmp('loglog',fittype)
			yfit = (10^b)*x.^m;
		end
		yr = yfit - Yfit;
		D(j) = sum(yr.^2)/(p*MSE); % Cook's Distance, sum over N excluding sample j
	end
	% Determine influential points and re-run regression analysis on kept samples
	x(D > Dthresh) = [];
	y(D > Dthresh) = [];
	if length(x) < 3 % this condition could be removed, but the assumption is a curve cannot be fit through less than 3 points
		rfitted.indx = D<Dthresh;
		rfitted.D = mean(D(D<Dthresh))/Dthresh;
		rfitted.coeffs = coeffs;
		rfitted.yfit = Yfit;
		rfitted.ci = [inf inf;inf inf];
		rfitted.mse = MSE;
		rfitted.r2 = 0;
	else
		[m,b,mse,r2,s] = logfit(x,y,fittype,'nograph');
		ci = tinv(1-alpha,s.df-1); % gives 95% confidence interval multiplier for 2-sided distribution
		if strcmp('linear',fittype)
			yfit = b + m*x;
			SE = stderr(x,y,yfit);
			coeffs = [b m];
			bnds = [...
				b-SE(2)*ci m-SE(1)*ci;...
				b+SE(2)*ci m+SE(1)*ci];
		elseif strcmp('logx',fittype)
			yfit = b + m*log10(x);
			SE = stderr(log10(x),y,yfit);
			coeffs = [b m];
			bnds = [...
				b-SE(2)*ci m-SE(1)*ci;...
				b+SE(2)*ci m+SE(1)*ci];
		elseif strcmp('logy',fittype)
			yfit = (10^b)*(10^m).^x;
			SE = stderr(x,log10(y),log10(yfit));
			coeffs = [10^b 10^m];
			bnds = [...
				10^(b-SE(2)*ci) 10^(m-SE(1)*ci);...
				10^(b+SE(2)*ci) 10^(m+SE(1)*ci)];
		elseif strcmp('loglog',fittype)
			yfit = (10^b)*x.^m;
			SE = stderr(log10(x),log10(y),log10(yfit));
			coeffs = [10^b m];
			bnds = [...
				10^(b-SE(2)*ci) (m-SE(1)*ci);...
				10^(b+SE(2)*ci) (m+SE(1)*ci)];
		end
		% output to structure
		rfitted.indx = D<Dthresh;
		rfitted.D = mean(D(D<Dthresh))/Dthresh;
		rfitted.coeffs = coeffs;
		rfitted.ci = bnds;
		rfitted.yfit = yfit;
		rfitted.mse = mse;
		rfitted.r2 = r2;
	end
end
%
function sep = stderr(x,y,yf)
	% sep = stderr(x,y,yf)
	% Calculate slope/intercept standard errors for linear regression
	N = length(x);
	if N ~= length(y) || N ~= length(yf)
		error('Inconsistent dimensions of inputs.')
	end
	R = y-yf;
	SXX = sum(var(x));
	MSE = sum(R.^2)/(N-2);
	sep = [...
	sqrt(MSE*sum(x.^2)/(N*SXX)),...
	sqrt(MSE/SXX)];
end
%
function [sf,stdsf]=logSmooth(f,varargin)
	%
	%  [sf,stdsf]=logSmooth(f:array,a)
	%
	% Log-uniform smoothing of f using bin-averaging technique.
	% Inputs:
	% f --> column vector needing smoothing. If f is array, operates on columns.
	% !! logSmooth assumes column vector or that inputs are to be operated over columns !!
	% a (optional) --> control width of smoothing windows, increasing a decreases smooth (default = 4)
	%
	% Outputs:
	% sf --> smoothed f (or array).
	% stdsf --> standard error of the means over each bin, x1.96 gives 95% confidence interval.
	%
	% Original code: Ioannis Kalogiros (5/3/2000) in Matlab v5.3
	% Including "a": D. Ortiz-Suslow (2018) Matlab v2018a/
	%
	if size(f,1) == 1
		f = f(:);
		warning('logSmooth only takes column vectors or operates over columns')
	end
	% Initialization
	a = 4; % default
	if nargin > 1
		a = varargin{1};
	end
	nf = size(f,1); % number of frequency bins
	m = a*(log(nf)/log(2)-1); 
	m = round(m) - a; % number of uniformly spaced averaging windows
	if nf<=a || m<=0 % cannot work on very narrow spectra
		warning('Asking to smooth a spectrum that is shorter than your smoothing window...')
		sf=f;
		return
	end
	% bin-averaging routine
	dl = log(nf-a)/m;
	sf(1:a,:) = 0.5*(f(1:a,:) + f(2:a+1,:)); % leading bin central difference
	stdsf(1:a,:) = sqrt(0.5*(f(1:a,:) - f(2:a+1,:)).^2); % leading bin central deviation
	l1p=0; l2p=0; np=[];
	for n = 1:m
	   l1 = (n-1)*dl;
	   l2 = l1+dl;
	   l1 = round(exp(l1))+a;
	   l2 = round(exp(l2))+a;
	   l1 = max(l1,1);
	   l1 = min(l1,nf);
	   l2 = max(l2,1);
	   l2 = min(l2,nf);
	   if l2==l1 && l2<nf
		   l2=l2+1;
	   end
	   k = l1:l2;
	   sf(n+a,:) = mean(f(k,:),1);
	   stdsf(n+a,:) =  std(f(k,:),1);
	   stdsf(n+a,:) = sem(stdsf(n+a,:),'flags',length(k));
	   if l1==l1p && l2==l2p 
		   np=[np;n+a];
	   end
	   l1p=l1;
	   l2p=l2;
	   bins(n,:) = [l1 l2];
	end
	sf(np,:)=[];
	stdsf(np,:) = [];
end
%
%.................................................................................
% below is an unaltered copy/paste of logfit.m
%.................................................................................
%% function [slope, intercept,MSE, R2, S] = logfit(x,y,varargin)
% This function plots the data with a power law, logarithmic, exponential
% or linear fit.
%
%   logfit(X,Y,graphType),  where X is a vector and Y is a vector or a
%               matrix will plot the data with the axis scaling determined
%               by graphType as follows: graphType-> xscale, yscale
%                  loglog-> log, log
%                    logx -> log, linear
%                    logy -> linear, log
%                  linear -> linear, linear
%               A line is then fit to the scaled data in a least squares
%               sense.
%               See the 'notes' section below for help choosing a method.
% 
%   logfit(X,Y), will search through all the possible axis scalings and
%               finish with the one that incurs the least error (with error
%               measured as least squares on the linear-linear data.)
% 
%   [slope, intercept, MSE, R2, S] = logfit(X,Y,graphType), returns the following:
%                slope: The slope of the line in the log-scale units.
%            intercept: The intercept of the line in the log-scale units.
%                  MSE: The mean square error between the 'y' data and the
%                       approximation in linear units.
%                   R2: The coefficient of determination
%                    S: This is returned by 'polyfit' and it allows you to
%                       be much fancier with your error estimates in the
%                       following way: (see polyfit for more information)
%                    >> S contains fields R, df, and normr, for the
%                    >> triangular factor from a QR decomposition of the
%                    >> Vandermonde matrix of x, the degrees of freedom,
%                    >> and the norm of the residuals, respectively. If the
%                    >> data y are random, an estimate of the covariance
%                    >> matrix of p is (Rinv*Rinv')*normr^2/df, where Rinv
%                    >> is the inverse of R. If the errors in the data y
%                    >> are independent normal with constant variance,
%                    >> polyval produces error bounds that contain at least
%                    >> 50% of the predictions.
% 
%   [graphType, slope, intercept, MSE, R2, S] = logfit(X,Y), if you choose
%                       not to pass a 'graphType' variable, then it will go
%                       ahead and select the one with the least square
%                       error. The firt parameter returned will be the
%                       graphType, with the following parameters in the
%                       usual order.
%               
%   logfit(X,Y,'PropertyName',PropertyValue), or
%   logfit(X,Y,graphType,'PropertyName',PropertyValue)
% 
%               see parameter options below
%__________________________________________________________________________ 
% USER PARAMETERS:
% 
% For skipping part of the data set:
%       'skip': skip 'n' rows from the beginning of the data set when
%               calculating the linear fit. Must be integer. Pass a negative
%               number to skip rows from the end instead of from the
%               beginning. All points will be plotted. 'num2skip'
%  'skipBegin': skip 'n' rows from the beginning when calculating the
%               linear fit similar to skip n. 'beginSkip'
%    'skipEnd': skip 'n' rows from the end, similar to skip -n 'endSkip'
% 
%__________________________________________________________________________ 
% For plotting in different styles
%   'fontsize': The fontsize of the axis, for axis tick labels and legend.
%               'font','fsize'
% 'markersize': The size of the marker for the points, 
% 'markertype': The type of marker for the points, such as 'o--' or '.r'
%               'markerstyle','markertype','marker'
% 
%  'linewidth': The width of the dashed line for the approximation
% 
%       'ftir': The approximation is plotted for a range around the
%               endpoints of the data set. By default it is 1/20 of the
%               range of the points. You may change this default by using
%               this parameter.
%               'fraction_to_increase_range','fractiontoincreaserange'
%__________________________________________________________________________ 
% Note the following sytax may also be used to specify 'graphtype'
%         'loglog','log','powerlaw'
%         'logx','logarithmic'
%         'logy','exponential','exp'
%         'linear','lin'
%__________________________________________________________________________ 
% Notes:
% The notes here will explain what the output means in terms of fitting
% functions depending on which method you use,
% 
% A power law relationship
% [slope, intercept] = logfit(x,y,'loglog');
%            yApprox = (10^intercept)*x.^(slope);
% 
% An exponential relationship
% [slope, intercept] = logfit(x,y,'logy');
%            yApprox = (10^intercept)*(10^slope).^x;
% 
% A logarithmic relationship
% [slope, intercept] = logfit(x,y,'logx');
%            yApprox = (intercept)+(slope)*log10(x);
% 
% A linear relationship
% [slope, intercept] = logfit(x,y,'linear');
%            yApprox = (intercept)+(slope)*x;
% 
%__________________________________________________________________________ 
% Examples:
% A power law, power 'a'
% a=2;
% x=(1:20)+rand(1,20); y=x.^a;
% power = logfit(x,y);
% % 
% A exponential relationship 
% a=3; x=(1:30)+10*rand(1,30); y=a.^x+100*rand(1,30);
% [graphType a] = logfit(x,y)
% base = 10^(a)
%
%............................................................
% NON-CRITICAL FUNCTIONING CHANGE MADE BY D. ORTIZ-SUSLOW
% include 'nograph' as on of the inputs and logfit.m will NOT plot the results of the fit
%...........................................................
%
% Thanks to Aptima inc. for  for giving me a reason to write this function.
% Thanks to Avi and Eli for help with designing and testing logfit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jonathan Lansey updated 2013,                                           %
%                   questions/comments welcome to Lansey at gmail.com     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [slope, intercept,MSE, R2, S, extra] = logfit(x,y,varargin)
	% The 'extra' is here in case 'graphtype' is not passed and needs to be
	% returned.
	extra=[];

	%% Check user inputed graphType, and standardize its value
	k=1;
	if isempty(varargin)
	    [slope, intercept,MSE, R2, S, extra] = findBestFit(x,y);
	    return;
    
	else % interpret all these possible user parameters here, so we can be more specific later.
	    switch lower(varargin{1}); % make all lowercase in case someone put in something different.
	        case {'logy','exponential','exp'}
	            graphType = 'logy';
	        case {'logx','logarithmic'}
	            graphType = 'logx';
	        case {'loglog','log','powerlaw'}
	            graphType = 'loglog';
	        case {'linear','lin'}
	            graphType = 'linear';
	        otherwise            
	            [slope, intercept, MSE, R2, S, extra] = findBestFit(x,y,varargin{:});
	            return;
	    end
	    k=k+1; % will usually look at varargin{2} later because of this
	end

	%% Set dynamic marker type defaults
	% for example, 'o' or '.' as well as size

	yIsMatrixFlag = size(y,1)>1 && size(y,2)>1; % There is more than one data point per x value
	markerSize=3;
	markerType = '.';
	markerColor = [0 0 0];

	if ~yIsMatrixFlag % check how many points there are
	    if length(y)<80 % relatively few points
	        markerType = 'o';
	        markerSize=5;
	    %   the following will overwrite markersize
	        if length(y)<30 % this number '30' is completely arbitrary
	            markerSize=7; % this '12' is also rather arbitrary
	        end
	%   else % there are many points, keep above defaults
	%         lineWidth=1;
	%         markerSize=5;
	    end
	end

	% markerLineWidth is always 2.

	%% Set static some defaults
	% before interpreting user parameters.
	fSize=15;
	num2skip=0; skipBegin = 0; skipEnd=0;
	ftir=20; %  = fraction_To_Increase_Range, for increasing where the green line is plotted
	lineColor = [.3 .7 .3]; % color of the line
	lineWidth=2;  % width of the approximate line

	%% Interpret extra user parameters

	while k <= length(varargin) && ischar(varargin{k})
	    switch (lower(varargin{k}))
	%       skipping points from beginning or end        
	        case {'skip','num2skip'}
	            num2skip = varargin{k+1};
	            k = k + 1;
	        case {'skipbegin','beginskip'}
	            skipBegin = varargin{k+1};
	            k = k + 1;
	        case {'skipend','endskip'}
	            skipEnd = varargin{k+1};
	            k = k + 1;

	%       Adjust font size
	        case {'fontsize','font','fsize'}
	            fSize = varargin{k+1};
	            k = k+1;

	%       Approx, line plotting
	        case {'ftir','fraction_to_increase_range','fractiontoincreaserange'}
	            ftir = varargin{k+1};
	            k = k+1;
                       
	%       Plotting style parameters        
	        case {'markersize'}
	            markerSize = varargin{k+1}; %forceMarkerSizeFlag=1;
	            k = k + 1;
	        case {'markercolor'}
	            markerColor = varargin{k+1};
	            k = k + 1;
	        case {'markertype','markerstyle','marker'}
	            markerType = varargin{k+1}; %forceMarkerTypeFlag=1;
	            k = k+1;
	        case {'linecolor','color'}
	            lineColor = varargin{k+1};
	            k = k+1;
	        case 'linewidth'
	            lineWidth = varargin{k+1};
	            k = k+1;
	        case 'nograph'
	            nographflag = 1;
	        otherwise
	            warning(['user entered parameter ''' varargin{k} ''' not recognized']);
	    end
	    k = k + 1;
	end

	%% Checks for user mistakes in input

	% data size and skip related errors/warnings
	    % Check they skipped an integer number of rows.
	    if round(skipBegin)~=skipBegin || round(skipEnd)~=skipEnd || round(num2skip)~=num2skip 
	        error('you can only skip an integer number of data rows');
	    end
	    if (skipEnd~=0 || skipBegin~=0) && num2skip~=0
	        warning('you have entered ambigious parameter settings, ''skipBegin'' and ''skipEnd'' will take priority');
	        num2skip=0;
	    end

	    if num2skip>0
	        skipBegin=num2skip;
	    elseif num2skip<0
	        skipEnd=-num2skip;
	    % else
	    %     num2skip==0; % so do nothing
	    end

	    % Check that the user has not skipped all of his/her data
	    if length(x)<1+skipEnd+skipBegin
	        error('you don''t have enough points to compute a linear fit');
	    end
	    if length(x)<3+skipEnd+skipBegin
	        warning('your data are meaningless, please go collect more points');
	    end
    
	% Data formatting errors and warnings    
	    % Check that 'x' is a vector
	    if size(x,1)>1 && size(x,2)>1 % something is wrong
	        error('Your x values must be a vector, it cannot be a matrix');
	    end

	    if yIsMatrixFlag % There is more than one data point per x value
	        if size(y,1)~=length(x)
	            error('the length of ''x'' must equal the number of rows in y');
	        end
	    else % y and x must be vectors by now
	        if length(x)~=length(y)
	            error('the length of ''x'' must equal the length of y');
	        end
	    end
    
	    if ~isnumeric(markerSize)
	        error('marker size must be numeric');
	    end

	% Helpful warning
	    if markerSize<=1
	        warning(['Your grandma will probably not be able to read your plot, '...
	                 'the markersize is just too small!']);
	    end



	%% Prepare y data by making it a properly oriented vector
	% skip rows as requested and create standard vectors (sometimes from matrices)

	x=x(:);
	x2fit=x(skipBegin+1:end-skipEnd);

	if yIsMatrixFlag % There is more than one data point per x value
	% note the '+1' so it can be used as an index value
	% This is the ones that will be used for fitting, rather than for plotting.
	    y2fit = y(skipBegin+1:end-skipEnd,:);
    
	    [x2fit,y2fit]= linearify(x2fit,y2fit);
	    [x,y]        = linearify(x,y);

	else % no need to linearify further
	    y=y(:);
	    y2fit=y(skipBegin+1:end-skipEnd);
	%     Note that 'x' is already forced to be a standard vector above
	end

	%% Check here for data that is zero or negative on a log scaled axis.
	% This is a problem because log(z<=0) is not a real number
	% This cell will remove it with a warning and helpful suggestion.
	% 
	% This warning can suggest you choose a different plot, or perhaps add 1 if
	% your data are large enough.
	% 
	% Note that this is done in order, so if by removing the 'y==0' values, you
	% also delete the 'x==0' values, then the 'x' warning won't show up. I
	% don't think this is of any concern though.
	% 
	switch graphType
	    case {'logy','loglog'}
	        yMask=(y<=0);
	        if sum(yMask)>0
	            yNegMask=(y<0);
	            if sum(yNegMask)>0 % there are proper negative values
	                warning(['values with y<=0 were removed.'...
	                         'Are you sure that ''logy'' is smart to take? '...
	                         'some ''y'' values were negative in your data.']);
            
	            else % just some zero values
	                if sum(y<10)/length(y) < (1/2) % if less than half your data is below than 10.
	                    warning(['values with y==0 were removed. '...
	                             'you may wish to add +1 to your data to make these points visible.']);
	                else % The numbers are pretty small, you don't want to add one.
	                    warning(['values with y==0 were removed. '...
	                             'Nothing you can really do about it sorry.']);
	                end
                
	            end
            
	            y=y(~yMask); y2Mask=(y2fit<=0); y2fit=y2fit(~y2Mask);
	            x=x(~yMask);                    x2fit=x2fit(~y2Mask);
	%             warning('values with y<=0 were removed. It may make suggest you add 1 to your data.')
	        end
	end

	switch graphType
	    case {'logx','loglog'}
	        xMask=(x<=0);
	        if sum(xMask)>0
            
	            xNegMask=(x<0);
	            if sum(xNegMask)>0 % there are proper negative values
	                warning(['values with x<=0 were removed.'...
	                         'Are you sure that ''logx'' is smart to take? '...
	                         'some ''x'' values were negative in your data.']);
            
	            else % just some zero values
	                if sum(x<10)/length(x) < (1/2) % if less than half your data is below than 10.
	                    warning(['values with x==0 were removed. '...
	                             'you may wish to add +1 to your data to make these points visible.']);
	                else % The numbers are pretty small, you don't want to add one.
	                    warning(['values with x==0 were removed. '...
	                             'Nothing you can really do about it sorry.']);
	                end
                
	            end
            
	            x=x(~xMask); x2Mask=(x2fit<=0); x2fit=x2fit(~x2Mask);
	            y=y(~xMask);                    y2fit=y2fit(~x2Mask);
	        end
	end

	%% FUNCTION GUTS BELOW
	%% set and scale the data values for linear fitting
	switch graphType
	    case 'logy'
	        logY=log10(y2fit);
	        logX=x2fit;
	    case 'logx'
	        logX=log10(x2fit);
	        logY=y2fit;
	    case 'loglog'
	        logX=log10(x2fit); logY=log10(y2fit);
	    case 'linear'
	        logX=x2fit; logY=y2fit;
	end

	%% Set the range that the approximate line will be displayed for

	if isempty(x2fit) || isempty(y2fit)
	    warning(['cannot fit any of your points on this ' graphType ' scale']);
	    slope=NaN; intercept=NaN; MSE = NaN; R2= NaN;
	    S=inf; % so that this is not used.
	    return;
	end


	range=[min(x2fit) max(x2fit)];
	% make this compatible with skipping some points.... don't know how yet....
	switch graphType
	    case {'logx','loglog'}
	        logRange=log10(range);
	        totRange=diff(logRange)+10*eps; % in case its all zeros...
	        logRange = [logRange(1)-totRange/ftir, logRange(2)+totRange/ftir];
	        ex = linspace(logRange(1),logRange(2),100); % note this is in log10 space

	    otherwise % logy, linear
	        totRange=diff(range);
	        range= [range(1)-totRange/ftir, range(2)+totRange/ftir];        
	        ex=linspace(range(1),range(2),100);
	end

	%% Do the linear fitting and evaluating

	[p, S] = polyfit(logX,logY,1);
	yy = polyval(p,ex);
	estY=polyval(p,logX); % the estimate of the 'y' value for each point.

	%% rescale the approximation results for plotting
	switch lower(graphType)
	    case 'logy'
	        yy=10.^yy;
	        estY=10.^estY; logY=10.^logY;% need to do this for error estimation
	    case 'logx'
	        ex=10.^ex;
	    case 'loglog'
	        yy=10.^yy;
	        ex=10.^ex;
	        estY=10.^estY; logY=10.^logY;% need to do this for error estimation
	    case 'linear'
	%         'do nothing';
	    otherwise
	%         'There is no otherwise at this point';
	end

	%% Calculate MSE and R2
	% Note that this is done after the data re-scaling is finished.
	MSE = mean((logY-estY).^2); % mean squared error.

	%     COVyhaty    = cov(estY,y); % = cov(estimated values, experimental values)
	%     R2        = (COVyhaty(2).^2) ./(var(estY).*var(y));
	%     
	tmp = corrcoef(estY,y).^2;
	R2 = tmp(2);

	if exist('nographflag','var')
	    %% set output data
	    % before returning
	    slope=p(1);
	    intercept = p(2);
	else
	    %% Ready the axis for plotting
	    % create or grab an axis before setting the scales
	    a=gca;
	    set(a,'fontsize',fSize);
	    holdState=ishold;

	    %% Plot the data
	    % This one is just to get the legend right
	    plot(x,y,markerType,'markersize',markerSize,'linewidth',2,'color',markerColor);

	    %% Plot the approximate line
	    hold('on'); % in case hold off was on before
	    plot(ex,yy,'--','linewidth',lineWidth,'color',lineColor);

	    %% Plot the points
	    % This time again just so it appears on top of the other line.
	    h = plot(x,y,markerType,'markersize',markerSize,'linewidth',2,'color',markerColor);
	    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 

	    %% Set the axis and to scale correctly
	    switch graphType
	        case 'logy'
	            set(a,'yscale','log');
	         case 'logx'
	            set(a,'xscale','log');
	        case 'loglog'
	            set(a,'xscale','log','yscale','log');
	        case 'linear'
	            set(a,'xscale','linear','yscale','linear');
	    end

	    %% Finish up some graph niceties
	    % fix the graph limits.
	    % no idea why this is always needed
	    axis('tight');

	    legend('data',[graphType ' fit'],'location','best'); legend('boxoff');
    
	    % reset hold state
	    if ~holdState
	        hold('off');
	    end
    
	    %% set output data
	    % before returning
	    slope=p(1);
	    intercept = p(2);
	end
end % function logfit over

%% linearify
% This function will take a vector x, and matrix y and arrange them so that
% y is a vector where each number in the i'th row of y has the value of the
% i'th number in 'x'
% This only works when the number of rows in y equals the number of
% elements in x. The new 'x' vector will be have length(y(:)) elements
function [x,y] = linearify(x,y)
	x=x(:); % just in case its not already a vector pointing this way.
	x=repmat(x,size(y,2),1);
	y=y(:);
	% if length(y)~=length(x)
	%     warning(['Look what you doin son, the length of ''x'' must equal the '...
	%            'number of rows in y to make this function useful'           ]);
	% end    
end

%% this checks to see which type of plot has the smallest error
% Then it will return and plot the results from the one with the least
% error. Note that 'graphType' is returned first, making all the following
% outputs shifted.
function [graphType, slope, intercept, MSE, R2, S] = findBestFit(x,y,varargin)
	% List of graph types to check
	testList={'loglog','logx','logy','linear'};
	MSE=zeros(4,1);

	warning('off'); hold('off'); % just so you don't have it repeating the warnings a million times
	for ii=1:4
	    [a,b,MSE(ii),c]=logfit(x,y,testList{ii},varargin{:});
	end
	warning('on')

	%% check for winning graphtype
	% the one with the minimum error wins.

	graphType=testList(MSE==min(MSE));
	switch length(graphType)
	    case 1
	%         no warning, nothing
	    case 2
	        warning([graphType{1} ' and ' graphType{2} ' had equal error, so ' graphType{1} ' was chosen)']);
	    case 3
	        warning([graphType{1} ', ' graphType{2} ' and ' graphType{3} ' had equal errors, so ' graphType{1} ' was chosen)']);
	    otherwise
	%         wow this will probably never happen
	        warning(['all graph types had equal error, ' graphType{1} ' was chosen']);
	end
	graphType=graphType{1};

	%% run it a last time to get results
	[slope, intercept, MSE, R2, S]=logfit(x,y,graphType,varargin{:});

end
%.................................................................................
%.................................................................................
