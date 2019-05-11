%...ariis_demo..................................................
%...............................................................
% This is short program that demonstrates the way to call ARIIS.
% the test data set comes from FLIP during CASPER-West..........
% author: D.G. Ortiz-Suslow
% revised: 02/14/19
% disclaimer: This software is provided as-is with no guarantee of functionality or fitness for the task for which it was designed to complete.
clear all; clf; close all
% load test data
load('./test.mat')
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
	Nfa = 1; % don't do any frequency averaging
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
	set(ax1,'Box','on','FontSize',25,'FontName','Times New Roman')
	datetick('x',15,'keepticks','keeplimits')
	xlabel('t','FontSize',25,'FontName','Times New Roman')
	ylabel('U','FontSize',25,'FontName','Times New Roman')
	title('ARIIS results for 30 mins of data','FontSize',25,'FontName','Times New Roman')
	legend(hp,' u',' v',' w')
	axis tight

	ax2 = subplot(2,1,2); hold(ax2,'on')
	hl = loglog(f,[Fuu Fww],'LineWidth',2);
	set(hl(2),'Color',get(hp(3),'Color'))
	set(ax2,'Box','on','FontSize',25,'FontName','Times New Roman',...
		'YScale','log','XScale','log','YLim',[1e-5 1e1],'Xlim',[1e-5 1e2])
	xlabel('f = nUadv/z','FontSize',25,'FontName','Times New Roman')
	ylabel('nS_{xx}/u_{*}^{2}(\sigma^{2}_{w})','FontSize',25,'FontName','Times New Roman')
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
		[' u: ' num2str(A(11),'%.4f') ' \pm ' num2str(range(A(13:14))/2,'%.6f')],...
		[' w: ' num2str(A(21),'%.4f') ' \pm ' num2str(range(A(23:24))/2,'%.6f')]};
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
	Nfa = 32; % do some frequency averaging (smoothing between modes 0 and 1 are DIFFERENT)
	fcutoff = 16; % only accept data below Nyquist, here we'll be extra conservative
	Cuw = spectf(u,w,dt,Nfa); % see below for spectf.m
	Cvw = spectf(v,w,dt,Nfa);
	Cuw(Cuw(:,1)>fcutoff,:) = [];
	Cvw(Cvw(:,1)>fcutoff,:) = [];
	% build ARIIS input structures
	inputs.n = Cuw(:,1);
	inputs.Suu = Cuw(:,2);
	inputs.Svv = Cvw(:,2);
	inputs.Sww = Cuw(:,3);
	constants.ust = ust;
	constants.varw = varw;
	constants.Uadv = Uadv;
	constants.z = zlevel;
	constants.dt = dt;
	m = str2num(answ);
	% call ARIIS
	A = ariis(m,inputs,constants);
	% surface layer frequency
	f = Cuw(:,1)*zlevel/Uadv;
	% nondimensional autovariance wind velocity spectra
	Fuu = Cuw(:,1).*Cuw(:,2)/(ust^2);
	Fww = Cuw(:,1).*Cuw(:,3)/varw;
	%..convert ARIIS frequency bandwidth to surface layer frequency
	lb = A(2)*zlevel/Uadv; % lower frequency bound
	ub = A(4)*zlevel/Uadv; % upper frequency bound
	%...............................................................
	hl = loglog(f,[Fuu Fww],'LineWidth',2);
	hold(gca,'on')
	set(gca,'Box','on','FontSize',25,'FontName','Times New Roman',...
		'YScale','log','XScale','log','YLim',[1e-5 1e1],'Xlim',[1e-5 1e2])
	xlabel('f = nUadv/z','FontSize',25,'FontName','Times New Roman')
	ylabel('nS_{xx}/u_{*}^{2}(\sigma^{2}_{w})','FontSize',25,'FontName','Times New Roman')
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
	text(6.5256,0.356,{'ARIIS-derived';'\sl{inertial subrange}'},...
		'FontSize',14,'FontName','Times New Roman',...
		'FontWeight','Bold','Interpreter','latex','Color',[0.5 0 0],...
		'BackGroundColor','w','EdgeColor',[0.5 0 0])
	leglabel = {...
		[' Kolmogorov = -5/3'],...
		[' u: ' num2str(A(11),'%.4f') ' \pm ' num2str(range(A(13:14))/2,'%.6f')],...
		[' w: ' num2str(A(21),'%.4f') ' \pm ' num2str(range(A(23:24))/2,'%.6f')]};
	legend([hkol;hl],leglabel,'Location','NorthWest')
	set(gcf,'Position',[723 9 957 946])
	xlim([Cuw(1,1)/5 Cuw(end,1)*5])
	disp(['Smoothing used by ARIIS matches smoothing in figure.'])
end

disp('ARIIS outputs....')
disp(['Subrange lo frequency bound = ' num2str(A(2),4) ' Hz'])
disp(['Subrange hi frequency bound = ' num2str(A(4),4) ' Hz'])
disp(['Isotropy (uw) = ' num2str(A(7),3)])
disp(['Isotropy (uv) = ' num2str(A(9),3)])
disp('........................')


%...............................................................
% demo complete
%...............................................................
%...............................................................
%...............................................................
%...below is a function for computing the power spectral density from a time series
%****D. Ortiz-Suslow did not write this.****
function S = spectf(x,y,dt,Nfa,a0)                              
	%        S = SPECTf(x,dt,Nfa)
	%        S = SPECTf(x,y,dt,Nfa)   cross spectrum
	% 
	% Frequency averaged power spectrum estimate,  GEOPHYSICAL NORMALIZATION
	% Trend is removed, Blackman-Harris window is used. K.K.Kahma 1990-05-19
	%
	%     x , y  = data vectors
	%     dt = sampling interval in seconds
	%     Nfa = number of elementary frequency bands which are averaged
	%
	%     S(:,1) = f      (1/second == Hz)
	%     S(:,2) = Sxx    (unit*unit*second)
	%
	% If cross spectrum is calculated 
	%     S(:,3) = Syy
	%     S(:,4) = Sxy
	%     S(:,5) = phase angle = 180/pi*atan2(-imag(Sxy),real(Sxy))
	%     S(:,6) = coherence   = abs(Sxy./sqrt(Sxx.*Syy))
	%
	%     positive phase means x leads y

	%        S = SPECTF(x,y,dt,Nfa,a0)
	% Elementary frequency bands 0:a0-1 (matlab index 1:a0) are ignored. 
	% Default a0 is 0, i.e. all bands including zero (mean value) are included.
	% 

	x = x(:).';      	% Make sure x is a row vector
	N  = max(size(x));      % Number of data points
	window=blackhar(N).';

	if max(size(y)) ~= N,
	   if (max(size(y)) == 1) || (nargin < 5)

	% ***************
	   % Spectrum 
	% ***************

	   if (nargin < 4), Nfa = 0; end    % default a0
	   if (nargin < 3), dt = 31; end    % default Nfa
	   a0 = Nfa; Nfa = dt; dt = y; 

	   Nfft=0;
	   maxb=0;
	   C=0; 
	   df=0;  % To define these variables before Xx'
	   Xx = fft(window.*detrend(x));
	   Nfft = length(Xx);                 % Number of points in FFT
	   maxb = floor(Nfft/2)+1;
	   Xx(maxb+1:Nfft)=[];
	   Xx(maxb) = Xx(maxb)/2;

	   C = dt/(Nfa*pi*norm(window)^2);    % Scaling coefficient
	   df = 2*pi/(dt*Nfft);
  
	   if Nfa==1
	      f = [a0:maxb-1]*df;
	      Pxx = (abs(Xx(a0+1:maxb)).^2)*C;
	   else
	      if Nfa > 20 
	        m=0; a=a0+1; b=a0+Nfa;
	        while b <= maxb
	           m=m+1;
	           Pxx(m) = sum(abs(Xx(a:b)).^2)*C;
	           f(m) = df*((a+b-2)/2);
	           a=a+Nfa; b=b+Nfa;
	        end
	      else
	        m=fix((maxb-a0) / Nfa);
	        f=([1:m]*Nfa+(a0-0.5-Nfa/2))*df;
	        b=a0+m*Nfa;
	        sx=zeros(Nfa,m);
	        sx(:) = abs(Xx(a0+1:b)).^2;  
	        Pxx=sum(sx)*C;
	      end
	      a=a0+1+m*Nfa;
	      if a <= maxb
	         m=m+1;
	         c = maxb+1-a; 
	         Pxx(m) = sum(abs(Xx(a:maxb)).^2)*C*Nfa/c;
	         f(m) = df*(a+maxb-2)/2;
	      end
	    end 
	    clear Xx window
	    S = [f/2/pi;2*pi*Pxx].';
 
	   else

	  error('x and y are not of same size'); end

	else

	% **********************
	   % Cross spectrum
	% **********************

	   if (nargin < 5), a0 = 0; end    % default a0
	   if (nargin < 4), Nfa = 31; end  % default Nfa

	   y = y(:).';
	   Nfft=0; maxb=0; C=0; df=0;
	   Xx = fft(window.*detrend(x));
	   Nfft = length(Xx);                 % Number of points in FFT
	   maxb = floor(Nfft/2)+1;
	   Xx(maxb+1:Nfft)=[];
	   Xx(maxb) = Xx(maxb)/2;

	   C = dt/(Nfa*pi*norm(window)^2);    % Scaling coefficient
	   df = 2*pi/(dt*Nfft);

	   Yy = fft(window.*detrend(y));
	   Yy(maxb) = Yy(maxb)/2;
	   Yy(maxb+1:Nfft)=[];

	   if Nfa==1
	      f = [a0:maxb-1]*df;
	      Pxx = (abs(Xx(a0+1:maxb)).^2)*C;
	      Pyy = (abs(Yy(a0+1:maxb)).^2)*C;
	      Pxy = (conj(Xx(a0+1:maxb)).*Yy(a0+1:maxb))*C;
	   else
	      if Nfa > 20
	         m=0; a=a0+1; b=a0+Nfa;
	         while b <= maxb
	            m=m+1;
	            Pxx(m) = sum(abs(Xx(a:b)).^2)*C;
	            Pyy(m) = sum(abs(Yy(a:b)).^2)*C;
	            Pxy(m) = sum(conj(Xx(a:b)).*Yy(a:b))*C;
	            f(m) = df*((a+b-2)/2);
	            a=a+Nfa; b=b+Nfa;
	         end
	      else
	         m=fix((maxb-a0) / Nfa);
	         f=([1:m]*Nfa+(a0-0.5-Nfa/2))*df;
	         b=a0+m*Nfa;
	         sx=zeros(Nfa,m);
	         sx(:) = abs(Xx(a0+1:b)).^2;  
	         Pxx=(sum(sx)*C);
	         sx(:) = abs(Yy(a0+1:b)).^2;  
	         Pyy=(sum(sx)*C);
	         sx(:) = conj(Xx(a0+1:b)).*Yy(a0+1:b);  
	         Pxy=(sum(sx)*C);
	         a=a0+1+m*Nfa;
	      end
   
	      if a <= maxb
	         m=m+1; 
	         c = maxb+1-a; 
	         Pxx(m) = sum(abs(Xx(a:maxb)).^2)*C*Nfa/c;
	         Pyy(m) = sum(abs(Yy(a:maxb)).^2)*C*Nfa/c;
	         Pxy(m) = sum(conj(Xx(a:maxb)).*Yy(a:maxb))*C*Nfa/c;
	         f(m) = df*(a+maxb-2)/2;
	      end
	   end
	   phase = atan2d(-imag(Pxy),real(Pxy)); % original code had (-) sign for imag
	   coh   = abs(Pxy./sqrt(Pxx.*Pyy));
	   clear Xx Yy window sx sy sxy
	   S = [f/2/pi;2*pi*Pxx;2*pi*Pyy;2*pi*Pxy;phase;coh].';
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