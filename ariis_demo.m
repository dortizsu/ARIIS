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
dt = mean(diff(tn))*(3600*24);
disp('....')
disp('ARIIS has two primary modes of operation:')
disp('Mode 0: input time series of 3D velocities (and scalars) or...')
disp('Mode 1: input 1 autovariance spectra from velocities (and scalars).')
disp('In mode 0, you can choose to apply some spectral smoothing (recommended).')
disp('In mode 1, spectra SHOULD have been smoothed before passing to ARIIS.')
disp('The demo allows you see ARIIS results after operating with mode 0 or 1...')
disp('....')
answ = input('Mode for ARIIS, time series (0) or spectral (1)? ','s');
if str2num(answ) == 0
	% build ARIIS input structures
	inputs.u = u;
	inputs.v = v;
	inputs.w = w;
	inputs.T = T + 273.15;
	inputs.C = T; % no concentration data in the example.
	constants.Uadv = Uadv;
	constants.WD = 0; % this is just for example purposes.
	constants.z = zlevel;
	constants.dt = dt;
	constants.htype = 1; % example uses IRGASON w/ CSAT-type sonic anemometer
	constants.dpl = 0.1537; % optical path-length of IRGASON gas analyzer
	constants.fcutoff = 5;
	constants.smoothing = 0; % log-uniform smoothing
	constants.nfa = 12; % degree of smoothing
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
	hp = plot(tn,[u v w],'-','LineWidth',2);
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
elseif str2num(answ) == 1
	varw = var(w);
	up = detrend(u);
	vp = detrend(v);
	wp = detrend(w);
    ii = mean(-up.*wp);
    jj = mean(-vp.*wp);
	ust = sqrt(sqrt(ii^2 + jj^2));
	%..initialize power spectrum call
	Nfa = 1; % do some frequency averaging (smoothing between modes 0 and 1 are DIFFERENT)
	Cuw = spectf(u,w,dt,Nfa); % see below for spectf.m
	Cvw = spectf(v,w,dt,Nfa);
	Ctt = spectf(T + 273.15,dt,Nfa);
	% smooth
	[smoothed,~] = logSmooth([Cuw(:,1) Cuw(:,2) Cvw(:,2) Cuw(:,3) Ctt(:,2)],12);
	% build ARIIS input structures
	inputs.n = smoothed(:,1);
	inputs.Suu = smoothed(:,2);
	inputs.Svv = smoothed(:,3);
	inputs.Sww = smoothed(:,4);
	inputs.Stt = smoothed(:,5);
	inputs.Scc = smoothed(:,5); % no concentration data in the example.
	constants.ust = ust;
	constants.varw = varw;
	constants.Uadv = Uadv;
	constants.WD = 0; % this is just for example purposes.
	constants.z = zlevel;
	constants.dt = dt;
	constants.htype = 1; % example uses IRGASON w/ CSAT-type sonic anemometer
	constants.dpl = 0.1537; % optical path-length of IRGASON gas analyzer
	constants.fcutoff = 5; % don't consider any frequencies above this limit
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
else
	error(['Unknown ARIIS mode selected...answer = ' answ])
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
% 
function [sf,sef]=logSmooth(f,varargin)
	%
	% [sf,sef] = logSmooth([f P])
	% [sf,sef] = logSmooth([f P],a)
	% [sf,sef] = logSmooth([f P],a,avgmode)
	%
	% Logarithmic smoothing sub-routine that can apply two techniques:
	% (0) the default method applies log-uniform mean smoothing. Smoothing is controlled by "a" (default value 4).
	% (1) bin-averaging is applied through uniform decade segmentation. Segmentation is controlled by "a"
	%
	% Inputs:
	% [f P] --> M x N array where "f" is the frequency amplitudes of the unsmoothed P, where P can be M x N-1
	% a --> smoothing factor. 
	% * for method (0): "a" is the slope of the frequency-dependent PDF of the number of samples per averaging bin used to calculate the smooth spectrum. a = 4 has twice as steep a slope as a = 8, therefore a = 4 smooths twice as fast as 8. Given the inverse proportionality of a to slope, a saturates in effectiveness at ~20-24. Therefore, values [2,24] is recommended.
	% * for method(1): "a"-1 is the number of segementations per decade of "f". This number is uniformly applied across all decades.
	%
	% Outputs:
	% sf --> smoothed [f P], all smoothing is done using a mean of the unsmoothed amplitudes within each smoothing bin/segment.
	% sef --> standard error within each bin. sef = 0 means that there was only 1 amplitude in a given bin.
	%
	% Original code (method 0) from Ioannis Kalogiros (2000) in Matlab v5.3
	% Additional code (method 1 & minor changes) D. Ortiz-Suslow (2020)
	%
	%
	if size(f,1) == 1
		f = f(:);
		warning('logSmooth only takes column vectors or operates over columns')
	end
	% Initialization
	a = 4; % default
	avgmode = 0; % default
	if nargin > 1
		a = varargin{1};
		if length(varargin) == 2
			avgmode = varargin{2};
		elseif length(varargin) > 2
			warning('Ignoring extra input...')
		end
	end
	nf = size(f,1); % number of frequency bins
	if avgmode == 0
		% original routine from JK
		% bin-averaging routine
		m = a*(log(nf)/log(2)-1); 
		m = round(m) - a; % number of uniformly spaced averaging windows
		if nf<=a || m<=0 % cannot work on very narrow spectra
			warning('Asking to smooth a spectrum that is shorter than your smoothing window...')
			sf=f;
			return
		end
		dl = log(nf-a)/m;
		sf(1:a,:) = 0.5*(f(1:a,:) + f(2:a+1,:)); % leading bin central difference
		sef(1:a,:) = sqrt(0.5*(f(1:a,:) - f(2:a+1,:)).^2); % leading bin central deviation
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
		   sef(n+a,:) =  std(f(k,:),1)/sqrt(length(k));
		   if l1==l1p && l2==l2p 
			   np=[np;n+a];
		   end
		   l1p=l1;
		   l2p=l2;
		   bins(n,:) = [l1 l2];
		end
		sf(np,:)=[];
		sef(np,:) = [];
	elseif avgmode == 1
		% uniform frequeny amplitudes per decade routine by DOS
		a = a + 1;
		n = log10(f(:,1));
		lb = floor(min(n));
		ub = ceil(max(n));
		n = lb:ub;
		nf = [];
		sf = [];
		sef = [];
		for ii = 1:length(n)-1
			lb = n(ii);
			ub = n(ii+1);
			dumn = logspace(lb,ub,a)';
			for jj = 1:a-1
				dumsf(jj,:) = mean(f(dumn(jj) <= f(:,1) & f(:,1) <= dumn(jj+1),2:end));
				dumstf(jj,:) = sem(f(dumn(jj) <= f(:,1) & f(:,1) <= dumn(jj+1),2:end));
			end
			nf = [nf;dumn(1:end-1) + diff(dumn)/2];
			sf = [sf;dumsf];
			sef = [sef;dumstf];
		end
		sef(isnan(sf)) = [];
		nf(isnan(sf)) = [];
		sf(isnan(sf)) = [];
		sf = [nf sf];
	else
		error('Averaging mode request not recognized: 0 or 1')
	end
end
% 
function hhh=vline(x,in1,in2)
	% function h=vline(x, linetype, label)
	% 
	% Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
	% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
	% label appears in the same color as the line.
	%
	% The line is held on the current axes, and after plotting the line, the function returns the axes to
	% its prior hold state.
	%
	% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
	% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
	% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
	% overridden by setting the root's ShowHiddenHandles property to on.
	%
	% Example:
	% h = vline(42,'g','The Answer')
	%
	% % returns a handle to a green vertical line on the current axes at x=42, and 
	% % creates a text object on the current axes, close to the line, which reads 
	% % "The Answer".
	%
	% % vline also supports vector inputs to draw multiple lines at once:
	%
	% vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
	%
	% % draws three lines with the appropriate labels and colors.
	% close all;

	% 
	% By Brandon Kuczenski for Kensington Labs.
	% brandon_kuczenski@kensingtonlabs.com
	% 8 November 2001

	Nx = length(x);
	if Nx==0,
	  h = [];
	elseif Nx>1  % vector input
	    for I=1:length(x)
	        switch nargin
	        case 1
	            linetype='r:';
	            label='';
	        case 2
	            if ~iscell(in1)
	                in1={in1};
	            end
	            if I>length(in1)
	                linetype=in1{end};
	            else
	                linetype=in1{I};
	            end
	            label='';
	        case 3
	            if ~iscell(in1)
	                in1={in1};
	            end
	            if ~iscell(in2)
	                in2={in2};
	            end
	            if I>length(in1)
	                linetype=in1{end};
	            else
	                linetype=in1{I};
	            end
	            if I>length(in2)
	                label=in2{end};
	            else
	                label=in2{I};
	            end
	        end
	        h(I)=vline(x(I),linetype,label);
	    end
	else
	    switch nargin
	    case 1
	        linetype='r:';
	        label='';
	    case 2
	        linetype=in1;
	        label='';
	    case 3
	        linetype=in1;
	        label=in2;
	    end

    
    
    
	    g=ishold(gca);
	    hold on

	    y=get(gca,'ylim');
	    h=plot([x x],y,linetype);
	    if length(label)
	        xx=get(gca,'xlim');
	        xrange=xx(2)-xx(1);
	        xunit=(x-xx(1))/xrange;
	        if xunit<0.8
	            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
	        else
	            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
	        end
	    end     

	    if g==0
	    hold off
	    end
	    set(h,'tag','vline','handlevisibility','off')
	end % else

	if nargout
	    hhh=h;
	end
end