function [fig1] = plot_gm2(titletext,gm_raw,gm_prc,cmp,acc_units,leg_labels,leg_colors,figid)
% plot raw vs processed time history

%suppress extra legend entities
warning('off','MATLAB:legend:IgnoringExtraEntries')

%default inputs
if nargin < 3; gm_prc     = {};                                 end
if nargin < 4; cmp        = {'Z','N','E'};                      end
if nargin < 5; acc_units  = 'g';                                end
if nargin < 6; leg_labels = {'Seed Motion','Processed Motion'}; end
if nargin < 7; leg_colors = {'#D95319','#0072BD'};              end
if nargin < 8; figid      = 1;                                  end

%colors
col_raw  = leg_colors{1}; %"#D95319"; 
if ~isempty(gm_prc); col_prc  = leg_colors{2}; end %"#0072BD";

% ------------------------------------------------------------------------- %
% -------------  Plot Time Histories and Spectra  ------------------------- %
% ------------------------------------------------------------------------- %
fig1 = figure(figid); clf;
orient tall;

%subfig locations
s_ind = {[01 10; 12 21; 23 32], ...
         [34 43; 45 54; 56 65] + 1, ...
         [67 76; 78 87; 89 98] + 2};

%iterate ground-motion components
for j = 1:length(gm_raw)
    % accel time his
    subplot(100,1,[s_ind{j}(1,1):s_ind{j}(1,2)])
    plot(gm_raw{j}.time,gm_raw{j}.acc,'-','color',col_raw,'linewidth',1); %rgb values for mathcad's brown is 128 64 0
    hold on;
    if ~isempty(gm_prc); plot(gm_prc{j}.time,gm_prc{j}.acc,'color',col_prc,'linewidth',0.7); end
    if j == 1; legend(leg_labels,'Location','NorthEast'); end
    % ylabel(['Acc. ',cmp{j},' (',acc_units,')']);
    ylabel({['Acc. ',cmp{j}],['(',acc_units,')']});
    dum=ylim; mx = max(max(dum),-min(dum));
    ylim([-mx mx]);
    % set(gca,'YTick',[-mx -mx/2 0 mx/2 mx]);
    set(gca,'XTickLabel','')
    graygrid()
    if j == 1; title(titletext,'FontWeight','bold'); end

    % velocity time his
    subplot(100,1,[s_ind{j}(2,1):s_ind{j}(2,2)])
    plot(gm_raw{j}.time,gm_raw{j}.vel,'-','color',col_raw,'linewidth',1);
    hold on;
    if ~isempty(gm_prc); plot(gm_prc{j}.time,gm_prc{j}.vel,'color',col_prc,'linewidth',0.7); end
    % ylabel(['Vel. ',cmp{j},' (cm/s)']);
    ylabel({['Vel. ',cmp{j}],'(cm/s)'});
    dum=ylim; mx = max(max(dum),-min(dum));
    ylim([-mx mx]);
    % set(gca,'YTick',[-mx -mx/2 0 mx/2 mx]);
    set(gca,'XTickLabel','')
    graygrid();

    % displacement time his
    subplot(100,1,[s_ind{j}(3,1):s_ind{j}(3,2)])
    plot(gm_raw{j}.time,gm_raw{j}.dis,'-','color',col_raw,'linewidth',1); 
    hold on; 
    if ~isempty(gm_prc); plot(gm_prc{j}.time,gm_prc{j}.dis,'color',col_prc,'linewidth',0.7); end
    if j == length(gm_raw); xlabel('Time (s)'); end
    % ylabel(['Disp. ',cmp{j},' (cm)']);
    ylabel({['Disp. ',cmp{j}],'(cm)'});
    dum=ylim; mx = max(max(dum),-min(dum));
    ylim([-mx mx]);
    % set(gca,'YTick',[-mx -mx/2 0 mx/2 mx]);
    if j ~= length(gm_raw); set(gca,'XTickLabel',''); end
    graygrid();
end


end