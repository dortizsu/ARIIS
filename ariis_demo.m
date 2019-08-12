%...ariis_demo..................................................
%...............................................................
% This is short program that demonstrates the way to call ARIIS.
% the test data set comes from FLIP during CASPER-West..........
% author: D.G. Ortiz-Suslow
% revised: 02/14/19
% disclaimer: This software is provided as-is with no guarantee of functionality or fitness for the task for which it was designed to complete.
clear all; clf; close all
% load test data
load('./ariis_test.mat')
% get sampling interval in seconds
dt = mean(diff(t))*(3600*24);
answ = input('Mode for ARIIS, time series (0) or spectral (1)? ','s');
if str2num(answ) == 0
	% build ARIIS input structures
	inputs.u = u;
	inputs.v = v;
	inputs.w = w;
	constants.Uadv = Uadv;
	constants.z = zlevel;
	constants.dt = dt;
	constants.fcutoff = 16;
	m = str2num(answ);
	% call ARIIS
	A = ariis(m,inputs,constants);
	%...............................................................
	%..Below is a brief demo of the ARIIS output relative to input time series
	%..initialize power spectrum call
	Nfa = 1 % don't do any frequency averaging
	Cuw = spectf(u,w,dt,Nfa); % see below for spectf.m
	Cvw = spectf(v,w,dt,Nfa);
	%..ARIIS processes using the nondimensional spectrum (see Kaimal et al. 1972 or Miyake et al. 1970)
	% friction velocit\y
	ust2 = (mean(detrend(u).*detrend(w))^2 + mean(detrend(v).*detrend(w))^2)^(1/2);
	% surface layer frequency
	f = Cuw(:,1)*zlevel/Uadv;
	% nondimensional autovariance wind velocity spectra
	Fuu = Cuw(:,1).*Cuw(:,2)/ust2;
	Fww = Cuw(:,1).*Cuw(:,3)/var(w);
	%..convert ARIIS frequency bandwidth to surface layer frequency
	lb = A(2)*zlevel/Uadv; % lower frequency bound
	ub = A(4)*zlevel/Uadv; % upper frequency bound
	%...............................................................
	% plotting routine
	figure(1)
	ax1 = subplot(2,1,1); hold(ax1,'on')
	hp = plot(t,[u v w],'-','LineWidth',2);
	set(ax1,'Box','on','FontSize',20,'FontName','Times New Roman','FontAngle','italic')
	datetick('x',15,'keepticks','keeplimits')
	xlabel('\sl{t}','FontSize',25,'FontName','Times New Roman')
	ylabel('\sl{U}','FontSize',25,'FontName','Times New Roman')
	title('ARIIS results for 30 mins of data','FontSize',25,'FontName','Times New Roman')
	legend(hp,' u',' v',' w')
	axis tight

	ax2 = subplot(2,1,2); hold(ax2,'on')
	hl = loglog(f,[Fuu Fww],'LineWidth',2);
	set(hl(2),'Color',get(hp(3),'Color'))
	set(ax2,'Box','on','FontSize',20,'FontName','Times New Roman','FontAngle','italic','MinorGridAlpha',0.15,...
		'YScale','log','XScale','log','YLim',[1e-5 1e1],'Xlim',[1e-5 1e2])
	xlabel('$$f = nU_{adv}/z$$','FontSize',25,'FontName','Times New Roman','Interpreter','Latex')
	ylabel('$$nS_{uu}/u_{*}^{2} \, | \, nS_{ww}/\sigma^{2}_{w}$$','FontSize',25,'FontName','Times New Roman','Interpreter','Latex')
	% add a Kolmogorov-like slope to axes 2
	kx = logspace(log10(lb),log10(ub));
	ky = 0.1*kx.^(-2/3);
	% add ARIIS output to axis 2
	vl = vline([lb ub],'r--'); set(vl,'LineWidth',4)
	px = [lb ub ub lb];
	py = [1e-6 1e-6 1e1 1e1];
	hpa = patch(px,py,'r','FaceAlpha',0.5); uistack(hpa,'bottom')
	hkol = loglog(kx,ky,'k-','LineWidth',4);
	% some additionals
	text(1.6363,2.6745,{'ARIIS-derived';'\sl{inertial subrange}'},...
		'FontSize',14,'FontName','Times New Roman',...
		'FontWeight','Bold','Interpreter','latex','Color',[0.5 0 0],...
		'BackGroundColor','w','EdgeColor',[0.5 0 0])
	leglabel = {...
		[' Kolmogorov = -5/3'],...
		[' u: ' num2str(A(10),'%.3f') ' \pm ' num2str(A(12),'%.3f')],...
		[' w: ' num2str(A(16),'%.3f') ' \pm ' num2str(A(18),'%.3f')]};
	legend([hkol;hl],leglabel,'Location','NorthWest')
	set(gcf,'Position',[723 9 957 946])
	disp(['ARIIS smooths spectra, Figure shows un-smoothed spectra.'])
else
	varw = var(w);
	up = detrend(u);
	vp = detrend(v);
	wp = detrend(w);
    ii = mean(-up.*wp);
    jj = mean(-vp.*wp);
	ust = sqrt(sqrt(ii^2 + jj^2));
	%..initialize power spectrum call
	Nfa = 1; % do some frequency averaging (smoothing between modes 0 and 1 are DIFFERENT)
	fcutoff = 16; % only accept data below Nyquist, here we'll be extra conservative
	Cuw = spectf(u,w,dt,Nfa); % see below for spectf.m
	Cvw = spectf(v,w,dt,Nfa);
	Cuw(Cuw(:,1)>fcutoff,:) = [];
	Cvw(Cvw(:,1)>fcutoff,:) = [];
	% smooth
	[smoothed,~] = logSmooth([Cuw(:,1) Cuw(:,2) Cvw(:,2) Cuw(:,3)],8);
	% build ARIIS input structures
	inputs.n = smoothed(:,1);
	inputs.Suu = smoothed(:,2);
	inputs.Svv = smoothed(:,3);
	inputs.Sww = smoothed(:,4);
	constants.ust = ust;
	constants.varw = varw;
	constants.Uadv = Uadv;
	constants.z = zlevel;
	constants.dt = dt;
	m = str2num(answ);
	% call ARIIS
	A = ariis(m,inputs,constants);
	% surface layer frequency
	f = smoothed(:,1)*zlevel/Uadv;
	% nondimensional autovariance wind velocity spectra
	Fuu = smoothed(:,1).*smoothed(:,2)/(ust^2);
	Fww = smoothed(:,1).*smoothed(:,4)/varw;
	%..convert ARIIS frequency bandwidth to surface layer frequency
	lb = A(2)*zlevel/Uadv; % lower frequency bound
	ub = A(4)*zlevel/Uadv; % upper frequency bound
	%...............................................................
	hl = loglog(f,[Fuu Fww],'LineWidth',2);
	hold(gca,'on')
	set(gca,'Box','on','FontSize',20,'FontName','Times New Roman','FontAngle','italic','MinorGridAlpha',0.15,...
		'YScale','log','XScale','log','YLim',[1e-5 1e1],'Xlim',[1e-5 1e2])
	grid on
	xlabel('$$f = nU_{adv}/z$$','FontSize',25,'FontName','Times New Roman','Interpreter','Latex')
	ylabel('$$nS_{uu}/u_{*}^{2} \, | \, nS_{ww}/\sigma^{2}_{w}$$','FontSize',25,'FontName','Times New Roman','Interpreter','Latex')
	% add a Kolmogorov-like slope to axes 2
	kx = logspace(log10(lb),log10(ub));
	ky = 0.1*kx.^(-2/3);
	% add ARIIS output to axis 2
	vl = vline([lb ub],'r--'); set(vl,'LineWidth',4)
	px = [lb ub ub lb];
	py = [1e-6 1e-6 1e1 1e1];
	hpa = patch(px,py,'r','FaceAlpha',0.25); uistack(hpa,'bottom')
	hkol = loglog(kx,ky,'k-','LineWidth',4);
	% some additionals
	ht = text(6.5256,0.356,{'ARIIS-derived';'\sl{inertial subrange}'},...
		'FontSize',14,'FontName','Times New Roman',...
		'FontWeight','Bold','Color',[0.5 0 0],...
		'BackGroundColor','w','EdgeColor',[0.5 0 0]);
	leglabel = {...
		[' Kolmogorov = -5/3'],...
		[' u: ' num2str(A(10),'%.3f') ' \pm ' num2str(A(12),'%.3f')],...
		[' w: ' num2str(A(16),'%.3f') ' \pm ' num2str(A(18),'%.3f')]};
	lh = legend([hkol;hl],leglabel,'Location','SouthWest');
	set(gcf,'Position',[1570 510 786 474])
	xlim([Cuw(1,1)/5 Cuw(end,1)*5])
	set(ht,'Position',[2.1958 0.59835 0])
	disp(['Smoothing used by ARIIS matches smoothing in figure.'])
end
disp('ARIIS outputs....')
disp(['Subrange lo frequency bound = ' num2str(A(2),'%.4f') ' Hz'])
disp(['Subrange hi frequency bound = ' num2str(A(4),'%.1f') ' Hz'])
disp(['Isotropy (uw) = ' num2str(A(7),'%.2f')])
disp(['Isotropy (uv) = ' num2str(A(9),'%.2f')])
disp('........................')
%...............................................................
% demo complete
%...............................................................
%...............................................................
%...............................................................
%...below is a function for computing the power spectral density from a time series
%****D. Ortiz-Suslow did not write this.****
function S = spectf(varargin)
	% S = SPECTf(x,dt)
	% S = SPECTf(x,dt,Nfa)
	% S = SPECTf(x,y,dt)
	% S = SPECTf(x,y,dt,Nfa)
	% S = SPECTF(x,y,dt,Nfa,a0)
	% calculate power spectral density and cross-sepctrum if needed
	% contributors: M. Donelan, W. Drennan, with modifications/cleaning from DOS
	Nfa = 31;
	a0 = 0;
	Nargsin = numel(varargin);
	x = varargin{1}; x = x(:)';
	N  = numel(x);
	window = blackhar(N)'; % window function
	% window=hamming(N)';
	%
	% figure out what you gave us...
	if Nargsin == 2
		dt = varargin{2};
	elseif Nargsin == 3
		if numel(varargin{2}) == 1
			dt = varargin{2};
			Nfa = varargin{3};
		else
			y = varargin{2};
			dt = varargin{3};
		end
	elseif Nargsin > 3
		y = varargin{2}; y = y(:)';
		dt = varargin{3};
		Nfa = varargin{4};
		if Nargsin > 4
			a0 = varargin{5};
		end
	end
	%
	% proceed to spectral calculations
	%
    Xx = fft(window.*detrend(x));
	% number of points in FFT
    Nfft = length(Xx);
    maxb = floor(Nfft/2)+1;
    Xx(maxb+1:Nfft) = [];
    Xx(maxb) = Xx(maxb)/2;
    C = dt/(Nfa*pi*norm(window)^2);    % Scaling coefficient
    df = 2*pi/(dt*Nfft);
	%
	if Nfa == 1
		f = [a0:maxb-1]*df;
		Pxx = (abs(Xx(a0+1:maxb)).^2)*C;
	else
		if Nfa > 20 
			% When Nfa is large enough this is as fast as vectorized
			% averaging and it requires far less memory    
	        m = 0;
			a = a0+1;
			b = a0+Nfa;
	        while b <= maxb
				m = m+1;
				Pxx(m) = sum(abs(Xx(a:b)).^2)*C;
				f(m) = df*((a+b-2)/2);
				a = a + Nfa;
				b = b + Nfa;
			end
		else
			m = fix((maxb - a0)/Nfa);
	        f = ((1:m)*Nfa + (a0 - 0.5 - Nfa/2))*df;
	        b = a0 + m*Nfa;
	        sx = zeros(Nfa,m);
	        sx(:) = abs(Xx(a0+1:b)).^2;  
	        Pxx=sum(sx)*C;
		end
        a = a0+1+m*Nfa;
        if a <= maxb
           m = m+1;
           c = maxb+1-a; 
           Pxx(m) = sum(abs(Xx(a:maxb)).^2)*C*Nfa/c;
           f(m) = df*(a+maxb-2)/2;
        end
	end
	S = [f/2/pi;2*pi*Pxx]';
	%
	% cross spectrum if necessary
	if exist('y','var')
 	   Yy = fft(window.*detrend(y));
 	   Yy(maxb) = Yy(maxb)/2;
 	   Yy(maxb+1:Nfft) = [];
	   if Nfa==1
		   Pyy = (abs(Yy(a0+1:maxb)).^2)*C;
		   Pxy = (conj(Xx(a0+1:maxb)).*Yy(a0+1:maxb))*C;
	   else
		   if Nfa > 20
			   m=0;
			   a=a0+1;
			   b=a0+Nfa;
			   while b <= maxb
				   m = m+1;
				   Pyy(m) = sum(abs(Yy(a:b)).^2)*C;
				   Pxy(m) = sum(conj(Xx(a:b)).*Yy(a:b))*C;
				   a = a + Nfa;
				   b = b + Nfa;
			   end
		   else
			   m = fix((maxb - a0)/Nfa);
			   b = a0 + m*Nfa;
			   sx = zeros(Nfa,m);
			   sx(:) = abs(Yy(a0+1:b)).^2;  
			   Pyy=(sum(sx)*C);
			   sx(:) = conj(Xx(a0+1:b)).*Yy(a0+1:b);  
			   Pxy=(sum(sx)*C);
			   a = a0 + 1 + m*Nfa;
		  end
		  if a <= maxb
		     m = m+1; 
		     c = maxb+1-a; 
		     Pyy(m) = sum(abs(Yy(a:maxb)).^2)*C*Nfa/c;
		     Pxy(m) = sum(conj(Xx(a:maxb)).*Yy(a:maxb))*C*Nfa/c;
		  end
		end
 	    phase = atan2d(-imag(Pxy),real(Pxy));
 	    coh = abs(Pxy./sqrt(Pxx.*Pyy));
 	    S = [S (2*pi*Pyy)' (2*pi*Pxy)' phase' coh'];
	end
end
% needed for spectf.m to run
function w = blackhar(n)
    %
    %  BLACKHAR(N) returns the N-point Blackman-Harris window as a column vector
    %  K.Kahma 1989-07-20
    m = (0:n-1)' * ( 2*pi/(n-1) ) ;
    w = (.35875 - .48829*cos(m) + .14128*cos(2*m) - 0.01168*cos(3*m) ) ;
end
function [sf,stdsf]=logSmooth(f,varargin)
	%
	% [sf,stdsf]=logSmooth(f,a)
	% [sf,stdsf]=logSmooth(f)
	%
	% Log-uniform smoothing of f using a bin-averaging technique.
	%
	% Inputs:
	% f --> column vector needing smoothing. If f is array, operates on columns.
	% a (optional) --> control width of smoothing windows, increasing "a" decreases smoothness (default = 4)
	%
	% Outputs:
	% sf --> smoothed array.
	% stdsf --> standard error of the means over each bin, x1.96 gives 95% confidence interval.
	% Contributor: John Kalogiros (2000) & DOS
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
	   stdsf(n+a,:) = stdsf(n+a,:)/sqrt(length(k));
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