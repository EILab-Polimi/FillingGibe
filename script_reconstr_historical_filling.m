clc
clear

%% load data and define parameters
addpath('data')
addpath('utils')
load 'reconstruct_hist_filling_database.mat'
gibe_masl = 660;    % [m] Gibe III level at the bottom of the dam 
lag = 16;           % [days] water travel time between gibe and turkana
delta = 3600*24;    % [sec/day]
MEF = 65;           % [m3/s] minimum environmental flow
Turb_disch = 102;   % [m3/s] turbine discharge capacity

%% initialize vectors for simulation
v_G(1:lag +1) = interp1qr(lsv_gibe(1,:), lsv_gibe(3,:), obs_level_gibe(1));  %gibe III volume
v_T(1:lag +1) = interp1qr(lsv_Turkana(1,:), lsv_Turkana(3,:), obs_level_turkana(1)); %turkana volume
v_T_natural(1:lag+1) = v_T(lag +1); % Turkana volume in a hypothetical simulation not comprising Gibe

%% 3 SEASON STRATEGY
% Bega (October to January), Belg (February to May) and Kiremt (June to September)
bega = zeros(size(inflow_gibe));
belg = zeros(size(inflow_gibe));
kiremt = zeros(size(inflow_gibe));

r_G(1:31) = MEF; % we consider dam filling starting with Belg season, in Feb 2015, while we only release MEF during January 

idx_begin = find(date_day == '1-Feb-2015');
r_G(32:idx_begin) = 0; 
bega(1:idx_begin) = 1;

idx_end_belg = find(date_day == '31-May-2015');
r_G(idx_begin+1 :idx_end_belg) = MEF; 
belg(idx_begin+1 :idx_end_belg) = 1;

idx_end_kiremt = find(date_day == '30-Sep-2015');
r_G(idx_end_belg + 1 : idx_end_kiremt) = Turb_disch*3; 
kiremt(idx_end_belg + 1 : idx_end_kiremt) = 1;

idx_end_bega = find(date_day == '31-Jan-2016');
r_G(idx_end_kiremt+1 :idx_end_bega) = Turb_disch*2;
bega (idx_end_kiremt+1 :idx_end_bega) = 1;

idx_end_belg = find(date_day == '31-May-2016');
r_G(idx_end_bega+1 :idx_end_belg) = Turb_disch*2; 
belg(idx_end_bega +1 :idx_end_belg) = 1;

idx_end_kiremt = find(date_day == '30-Sep-2016');
r_G(idx_end_belg + 1 : idx_end_kiremt) = Turb_disch*6;
kiremt (idx_end_belg + 1 : idx_end_kiremt) = 1;

idx_end_bega = find(date_day == '31-Jan-2017');
r_G(idx_end_kiremt+1 :idx_end_bega) = Turb_disch*3; 
bega (idx_end_kiremt+1 :idx_end_bega) = 1;

idx_end_belg = find(date_day == '31-May-2017');
r_G(idx_end_bega+1 :idx_end_belg) = Turb_disch*3; 
belg(idx_end_bega +1 :idx_end_belg) = 1;

idx_end_kiremt = find(date_day == '30-Sep-2017');
r_G(idx_end_belg + 1 : idx_end_kiremt) = Turb_disch*5; 
kiremt (idx_end_belg + 1 : idx_end_kiremt) = 1;

idx_end_bega = find(date_day == '31-Jan-2018');
r_G(idx_end_kiremt+1 :idx_end_bega) = Turb_disch*4; 
bega (idx_end_kiremt+1 :idx_end_bega) = 1;

idx_end_belg = find(date_day == '31-May-2018');
r_G(idx_end_bega+1 :idx_end_belg) = Turb_disch*4; 
belg(idx_end_bega +1 :idx_end_belg) = 1;

idx_end_kiremt = find(date_day == '30-Sep-2018');
r_G(idx_end_belg + 1 : idx_end_kiremt) = Turb_disch*7; 
kiremt (idx_end_belg + 1 : idx_end_kiremt) = 1;

r_G(idx_end_kiremt + 1 :length(inflow_gibe)) = Turb_disch*2;
bega(idx_end_kiremt  + 1 :length(inflow_gibe)) = 1;

%% simulation of the release strategy
H = length(inflow_gibe);
for i = lag+1:H - 1 
    v_G(i+1) = v_G(i) + (inflow_gibe(i+1) - r_G(i+1) - evap_gibe(i+1))*delta;
    v_T(i+1) = v_T(i) + (r_G(i+1 - lag) +  inflow_turkana(i+1) - evap_turkana(i+1))*delta;
    v_T_natural(i+1) = v_T_natural(i) + (inflow_gibe(i+1 - lag) - evap_gibe(i+1-lag) + inflow_turkana(i+1) - evap_turkana(i+1))*delta ; 
end

l_G = interp1(lsv_gibe(3,:), lsv_gibe(1,:), v_G) + gibe_masl; 
l_T = interp1(lsv_Turkana(3,:), lsv_Turkana(1,:), v_T); 
l_T_natural = interp1(lsv_Turkana(3,:), lsv_Turkana(1,:), v_T_natural); 

%% plot results from Jan 1st 2015 
% (nb. simulations begin in Dec 2014 to account for lag time in initialization)
id1Jan = find(date_day == '1-Jan-2015');
l_T = l_T(id1Jan:end);
l_G = l_G(id1Jan:end);
r_G = r_G(id1Jan:end);
l_T_natural = l_T_natural(id1Jan:end);
date_day = date_day(id1Jan:end); 
bega = bega(id1Jan:end);
belg = belg(id1Jan:end);
kiremt = kiremt(id1Jan:end);
inflow_gibe = inflow_gibe(id1Jan:end);

%% plot simulation results

figure; 
subplot(3,1,1)
date_day.Format = 'MM/yy';

area(date_day, bega*2.8e+03, ...
    'FaceAlpha', 0.3, 'EdgeColor','none', ...
    'FaceColor',[0 0.447058826684952 0.74117648601532]); %[0.929411768913269 0.694117665290833 0.125490203499794]);

hold on

area(date_day, belg*2.8e+03, ...
    'FaceAlpha', 0.45, 'EdgeColor','none', ...
    'FaceColor',[0 0.447058826684952 0.74117648601532]); %[0.929411768913269 0.694117665290833 0.125490203499794]);

area(date_day, kiremt*2.8e+03, ...
    'FaceAlpha', 0.8, 'EdgeColor','none', ...
    'FaceColor',[0 0.447058826684952 0.74117648601532]); %[0.929411768913269 0.694117665290833 0.125490203499794]);

area(date_day, inflow_gibe,'DisplayName','inflow',...
    'FaceColor', [1 1 1], ...% [0 0.376470595598221 0.745098054409027],...
    'EdgeColor','k',...
    'LineWidth',0.5);
ylabel('inflow and release [m^3/s]')
hold on
plot(date_day, r_G,'k', 'LineWidth',2)
legend('Bega', 'Belg', 'Kiremt', 'inflow', 'release', 'Location', 'northwest')
title('Gibe III inflow and release [m^3/s]')
set(gca,'FontSize', 14)
ylim([0 2.5e+03]) 
xl = [datetime(2015,1,1), datetime(2018,10,30)];
xlim(xl)
xtick = datetime(2015, 1, 1)+calmonths(0:4:12*4-1);
xtick.Format = 'MM/yy';
set(gca, 'xtick', xtick)

subplot(3,1,2)
plot(date_month,obs_level_gibe + gibe_masl,'Color',  rgb('ForestGreen'), 'LineWidth',2)
hold on 
plot(date_month(29:5:34),obs_level_gibe(29:5:34) + gibe_masl,'Color', rgb('ForestGreen'),'LineStyle', '--', 'LineWidth',2,'HandleVisibility','off')  % observed data unavailable for eccessive cloud cover
ylabel('Level [masl]')
plot(date_day, l_G, 'k', 'LineWidth',2 );
plot(date_month, repmat(861, length(date_month), 1), 'k--', 'LineWidth',0.4) 
plot(date_month, repmat(800, length(date_month), 1), 'k-.', 'LineWidth',0.4) 
plot(date_month, repmat(892, length(date_month), 1), 'k-.', 'LineWidth',0.4 ,'HandleVisibility','off') 

legend('observed','simulated', 'normal operating level', 'min-max operating level',  'Location','southeast')
title('Gibe III level')
set(gca,'FontSize', 14)
xl = [datetime(2015,1,1), datetime(2018,10,30)];
xlim(xl)
grid on
set(gca, 'xtick', xtick)

dt.Format = 'MM/yy';
subplot(3,1,3)
plot(date_day, obs_level_turkana,'Color', rgb('ForestGreen'), 'LineWidth',2)
hold on
plot(date_day, l_T, 'k', 'LineWidth',2)
hold on

plot(date_day, l_T_natural, 'k:', 'LineWidth',2)
legend('observed', 'simulated', 'pre Gibe III',  'Location','northwest')
ylabel('Level [masl]')
title('Lake Turkana level')
set(gca,'FontSize', 14)
grid on
ylim([362.5 366])
xl = [datetime(2015,1,1), datetime(2018,10,30)];
xlim(xl)
set(gca, 'xtick', xtick)



