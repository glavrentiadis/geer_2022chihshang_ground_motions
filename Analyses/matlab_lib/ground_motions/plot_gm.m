function [fig1, fig2] = plot_gm(titletext,gm_raw,gm_prc,gm_corr,acc_units,leg_labels,leg_colors,figid)
% plot raw vs processed time history

%suppress extra legend entities
warning('off','MATLAB:legend:IgnoringExtraEntries')

%default inputs
if nargin < 3; gm_prc    = [];  end
if nargin < 4; gm_corr   = [];  end
if nargin < 5; acc_units = 'g'; end
if nargin < 6; leg_labels = {'Seed Motion','Processed Motion','Correction'}; end
if nargin < 7; leg_colors = {'#D95319','#0072BD','#77AC30'}; end
if nargin < 8; figid = 1; end

%colors
col_raw  = leg_colors{1}; %"#D95319";
if ~isempty(gm_prc);  col_prc  = leg_colors{2}; end %"#0072BD";
if ~isempty(gm_corr); col_corr = leg_colors{3}; end %"#77AC30";

% ------------------------------------------------------------------------- %
% -------------  Plot Time Histories and Spectra  ------------------------- %
% ------------------------------------------------------------------------- %
fig1 = figure(figid); clf;
orient tall;

%subfig locations
s_ind = [1 18; 21 36; 39 55; 63 98];


% accel time his
subplot(100,1,[s_ind(1,1):s_ind(1,2)])
plot(gm_raw.time,gm_raw.acc,'-','color',col_raw,'linewidth',1); 
hold on;
if ~isempty(gm_prc);  plot(gm_prc.time,gm_prc.acc,'-','color',col_prc,'linewidth',0.7); end
if ~isempty(gm_corr); plot(gm_corr.time,gm_corr.acc,'--','color',col_corr,'linewidth',0.7); end
legend(leg_labels,'Location','NorthEast');
ylabel(['Acc. (',acc_units,')']);
dum=ylim; 
mx = max(max(dum),-min(dum));
ylim([-mx mx]);
% set(gca,'YTick',[-mx -mx/2 0 mx/2 mx]);
set(gca,'XTickLabel','')
graygrid()
title(titletext,'FontWeight','bold')

% velocity time his
subplot(100,1,[s_ind(2,1):s_ind(2,2)])
plot(gm_raw.time,gm_raw.vel,'-','color',col_raw,'linewidth',1);
hold on;
if ~isempty(gm_prc);  plot(gm_prc.time,gm_prc.vel,'-','color',col_prc,'linewidth',0.7); end
if ~isempty(gm_corr); plot(gm_corr.time,gm_corr.vel,'--','color',col_corr,'linewidth',0.7); end
ylabel('Vel (cm/s)');
dum=ylim; 
mx = max(max(dum),-min(dum));
ylim([-mx mx]);
% set(gca,'YTick',[-mx -mx/2 0 mx/2 mx]);
set(gca,'XTickLabel','')
graygrid();

% displacement time his
subplot(100,1,[s_ind(3,1):s_ind(3,2)])
plot(gm_raw.time,gm_raw.dis,'-','color',col_raw,'linewidth',1); 
hold on; 
if ~isempty(gm_prc);  plot(gm_prc.time,gm_prc.dis,'-','color',col_prc,'linewidth',0.7); end
if ~isempty(gm_corr); plot(gm_corr.time,gm_corr.dis,'--','color',col_corr,'linewidth',0.7); end
xlabel('Time (s)');
ylabel('Disp (cm)');
dum=ylim; mx = max(max(dum),-min(dum));
ylim([-mx mx]);
% set(gca,'YTick',[-mx -mx/2 0 mx/2 mx]);
graygrid();

% acc resp spectra
subplot(100,1,[s_ind(4,1):s_ind(4,2)])
loglog(gm_raw.per,gm_raw.ars,'-','color',col_raw,'linewidth',1); 
hold on;
if ~isempty(gm_prc); loglog(gm_prc.per,gm_prc.ars,'b','linewidth',0.7); end
xlabel('Period (s)')
Tmax = max (10, max(gm_raw.per)) ;
xlim ([0.01 Tmax]);
ylabel(['Spectral Acceleration (',acc_units,')'])
% ylabel(['Spectral Acceleration - 5% Damping (',acc_units,')'])
% text(0.011,0.08,'5% Damping','BackgroundColor','white','EdgeColor','black');
% text(0.014,0.26,'5% Damping','BackgroundColor','white','EdgeColor','black','Units','normalized');
% text(0.86,0.82,'5% Damping','BackgroundColor','white','EdgeColor','black','Units','normalized');
text(0.022,0.21,'5% Damping','BackgroundColor','white','EdgeColor','black','Units','normalized');
legend('Seed Motion','Processed Motion','Location','SouthWest');
xlabel('Spectral Periods (s)');
graygrid();


% ------------------------------------------------------------------------- %
% ---------------  Plot Husid Plots  -------------------------------------- %
% ------------------------------------------------------------------------- %
fig2 = figure(figid+1);clf
orient tall

%husid plot
subplot(100,1,[1 30])
plot(gm_raw.time,gm_raw.IaNorm_ofT,'-','color',col_raw,'linewidth',1); 
hold on;
if ~isempty(gm_prc); plot(gm_prc.time,gm_prc.IaNorm_ofT,'-','color',col_prc,'linewidth',0.7); end
legend(leg_labels,'Location','NorthEast');
ylabel('Normalized Arias Acc Intensity');
title(titletext,'FontWeight','bold')
ylim([0 1])
graygrid()
tmp=xlim;
space = 0.07;
%seed motion iffo
text(tmp(2)*1.01,1.0,leg_labels{1},...
	'Color', col_raw, 'FontWeight','bold',...
	'FontSize',8)
text(tmp(2)*1.01,1.0-1*space,...
	['Ia = ',num2str(gm_raw.Ia,'%1.1f'),' m/s'], ...
	'Color', col_raw, 'FontWeight','bold',...
	'FontSize',8)
text(tmp(2)*1.01,1.0-2*space,...
	['D_5_-_9_5 = ',num2str(gm_raw.D5_95,'%1.1f'),' s'], ...
	'Color', col_raw, 'FontWeight','bold',...
	'FontSize',8)
text(tmp(2)*1.01,1.0-3*space,...
	['D_5_-_7_5 = ',num2str(gm_raw.D5_75,'%1.1f'),' s'], ...
	'Color', col_raw, 'FontWeight','bold',...
	'FontSize',8)
text(tmp(2)*1.01,1.0-4*space,...
	['pga = ',num2str(gm_raw.pga,'%1.1f'),acc_units], ...
	'Color', col_raw, 'FontWeight','bold',...
	'FontSize',8)
text(tmp(2)*1.01,1.0-5*space,...
	['pgv = ',num2str(gm_raw.pgv,'%1.1f'),' cm/s'], ...
	'Color', col_raw, 'FontWeight','bold',...
	'FontSize',8)
text(tmp(2)*1.01,1.0-6*space,...
	['pgd = ',num2str(gm_raw.pgd,'%1.1f'),' cm'], ...
	'Color', col_raw, 'FontWeight','bold',...
	'FontSize',8)
%processed motion info
if ~isempty(gm_prc)
    text(tmp(2)*1.01,0.5,leg_labels{2},...
	    'Color', col_prc, 'FontWeight','bold',...
	    'FontSize',8)
    text(tmp(2)*1.01,0.5-1*space,...
	    ['Ia = ',num2str(gm_prc.Ia,'%1.1f'),' m/s'], ...
	    'Color', col_prc, 'FontWeight','bold',...
	    'FontSize',8)
    text(tmp(2)*1.01,0.5-2*space,...
	    ['D_5_-_9_5 = ',num2str(gm_prc.D5_95,'%1.1f'),' s'], ...
	    'Color', col_prc, 'FontWeight','bold',...
	    'FontSize',8)
    text(tmp(2)*1.01,0.5-3*space,...
	    ['D_5_-_7_5 = ',num2str(gm_prc.D5_75,'%1.1f'),' s'], ...
	    'Color', col_prc, 'FontWeight','bold',...
	    'FontSize',8)
    text(tmp(2)*1.01,0.5-4*space,...
	    ['pga = ',num2str(gm_prc.pga,'%1.1f'),acc_units], ...
	    'Color', col_prc, 'FontWeight','bold',...
	    'FontSize',8)
    text(tmp(2)*1.01,0.5-5*space,...
	    ['pgv = ',num2str(gm_prc.pgv,'%1.1f'),' cm/s'], ...
	    'Color', col_prc, 'FontWeight','bold',...
	    'FontSize',8)
    text(tmp(2)*1.01,0.5-6*space,...
	    ['pgd = ',num2str(gm_prc.pgd,'%1.1f'),' cm'], ...
	    'Color', col_prc, 'FontWeight','bold',...
	    'FontSize',8)
end

%cumulative vel plot
subplot(100,1,[35 64])
plot(gm_raw.time,gm_raw.IvNorm_ofT,'-','color',col_raw,'linewidth',1); 
hold on;
if ~isempty(gm_prc); plot(gm_prc.time,gm_prc.IvNorm_ofT,'-','color',col_prc,'linewidth',0.7); end
ylabel('Normalized Arias Vel Intensity');
ylim([0 1])
graygrid()

%cumulative disp plot
subplot(100,1,[69 98])
plot(gm_raw.time,gm_raw.IdNorm_ofT,'-','color',col_raw,'linewidth',1); 
hold on;
if ~isempty(gm_prc); plot(gm_prc.time,gm_prc.IdNorm_ofT,'-','color',col_prc,'linewidth',0.7); end
ylabel('Normalized Arias Disp Intensity');
xlabel('Time (s)');
ylim([0 1])
graygrid()

end