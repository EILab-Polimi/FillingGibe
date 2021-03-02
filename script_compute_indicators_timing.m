% script for computing indicators of alternative filling timings, and generating figures 3b and S6

clc
clear

load data/alternative_timings_dataset.mat

delta = 3600*24;
lag_g3K = 5;
lag_KT = 7;
qmaxGIII = 10*95;
effic = 0.88;  
Nyst = 2;


for year =  2007:2016 
%% select filling year
start = ['1-Jan-', num2str(year)];
endy = ['31-Dec-', num2str(year+Nyst-1) ];
ini = find(d_ == start);
fin = find(d_ == endy);

inflow_gibe_ = inflow_gibe(ini:fin);
inflow_lateral_ = inflow_lateral(ini:fin);
inflow_turKer_ = inflow_turKer(ini:fin);

%% initialize simulation variables
v(lag_g3K + lag_KT +1) = 434137000; %  
h(lag_g3K + lag_KT +1) = interp1(lsv_gibe(3,:), lsv_gibe(1,:), v(lag_g3K + lag_KT +1));
s(lag_g3K + lag_KT +1) = interp1(lsv_gibe(3,:), lsv_gibe(2,:), v(lag_g3K + lag_KT +1));
volTurk(lag_g3K + lag_KT +1) = interp1(lsv_Turkana(1,:), lsv_Turkana(3,:),turkana_obs_level(ini+365*7));
s_T(lag_g3K + lag_KT +1) = interp1(lsv_Turkana(3,:), lsv_Turkana(2,:), volTurk(lag_g3K + lag_KT +1));
delta_inflow = nan(size(volTurk));

%% Gibe III and Turkana simulations
for i = lag_g3K + lag_KT +1:length(inflow_gibe_) - 1 

    v(i+1) = (v(i) + inflow_gibe_(i+1)*delta - release_gibe(i+1)*delta) - evap_gibe(i+1)*s(i);
    e(i+1) = evap_gibe(i+1)*s(i);
    s(i+1) = interp1(lsv_gibe(3,:), lsv_gibe(2,:), v(i+1));
    h(i+1) = interp1(lsv_gibe(3,:), lsv_gibe(1,:), v(i+1));
    delta_inflow(i+1) = release_gibe(i+1 - lag_g3K - lag_KT) + inflow_turKer_(i+1) + inflow_lateral_(i+1);
    
    volTurk(i+1) = volTurk(i) + (release_gibe(i+1 - lag_g3K - lag_KT) +  inflow_turKer_(i+1) + inflow_lateral_(i+1))*delta - evap_turk(i+1)*s_T(i);
    s_T(i+1) = interp1(lsv_Turkana(3,:), lsv_Turkana(2,:), volTurk(i+1));
 
end


%% compute objectives

% hydropower
r_G = release_gibe(lag_g3K + lag_KT +1:365*Nyst); 
h = h(lag_g3K + lag_KT +1:365*Nyst); 
q_turb=min(release_gibe(lag_g3K + lag_KT +1:365*Nyst), qmaxGIII);
g_hyd_GIII = 24*effic.*9.81.*q_turb.*h./10^3;   % MWh 
Ghyd(year-2006) = nansum(g_hyd_GIII )/Nyst/1000;  %GWh/year 

% turkana level drop
l_T = interp1(lsv_Turkana(3,:), lsv_Turkana(1,:), volTurk); 
DeltaTurk(year-2006) = l_T(end) - l_T(lag_g3K + lag_KT +1);

% flood pulse
[~, idmaxrec] = max(omorateRegime);
Grec(year-2006, :) = [ max( delta_inflow(idmaxrec-15:idmaxrec+15) ),max( delta_inflow(idmaxrec-30 + 365 :idmaxrec+30 + 365) ) ];
Grec_m = mean(Grec');

% final Gibe III level
hend(year-2006)=h(end);

end


lab = {'Hydropower production [GWh/y]', 'Final Gibe III level [masl]', 'Turkana level drop [m]', 'Flood pulse [m^3/s]'}; 
ylab = {'GWh/y','masl','m','m^3/s'}; 
obj = [Ghyd', hend'+660 , DeltaTurk', Grec_m']; % all to be maximized


%% plots - all years

lab_y = {'2007','2008','2009','2010','2011','2012','2013','2014','2015', '2016' };

figure;
subplot(211)
plot( datemExt(97:216), h_sum(97:216), 'k', 'LINEWIDTH', 2 )
ylabel('[mm]')
title('Harmonical Precipitation trend')
set(gca, 'FontSize', 14); grid on; box on;

    
subplot(212)
scatter( obj(:,1),  obj(:,3), 100, 'filled')
hold on; 
scatter( obj([2 3 4 9 10],1),  obj([2 3 4 9 10],3), 100, 'r', 'filled')
xlabel('HP production [GWh/year]')
ylabel('Turkana Level Drop [m]')
set(gca, 'FontSize', 15); grid on; box on;
obj = [obj, [2007:2016]' ];

%% plots - all indicators

lab_y = {'2007','2009','2013','2015'};
years = [ 2007 2009 2013 2015 ]-2006; 
obj = obj(years, :);

figure; 
for i = 1:length(years)
    subplot(2,2,i)
    for j = 1:4
    bar( j, obj(j,i), 'FaceColor', color_sch(j,:) ); 
    set(gca,'XTick', 1:j, 'XTickLabel',lab_y(1:j)); hold on;
    
    end
    ylabel(ylab(i))
    if i == 2
        %ylim([0 900])
    else if i == 4
            ylim([0 1200])
        end
    end
    set(gca, 'FontSize', 12); grid on;
    title(lab(i))
end



