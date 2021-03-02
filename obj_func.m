function J = obj_func( k )

% Constrain k so that fraction r_g3d during dry seasons is lower than
% fraction of normal season, in turn lower than wet season. 

const = 1;
if k(2)>= k(1)
    if k(3)>= k(2) 
            const = 1;
    end
end

J = [Inf Inf Inf];

if const == 1

global  reservoir_specifics sim policy SPEI sys_param spei_date
norm_operating_level = reservoir_specifics.norm_operating_level;
min_rel              = reservoir_specifics.min_rel;
max_rel              = reservoir_specifics.max_rel;
lag                  = reservoir_specifics.lag; 
dead_storage_g3      = reservoir_specifics.dead_storage_g3; 
lsv_gibeIII          = reservoir_specifics.lsv_gibe; 
lsv_Turkana          = reservoir_specifics.lsv_turkana; 
effic                = reservoir_specifics.effic;
min_operating_level  = reservoir_specifics.min_operating_level;
g3_tail              = reservoir_specifics.g3_tail;

inflow_gibe          = sim.inflow_gibe;
inflow_turkana       = sim.inflow_turkana;
evap_gibe            = sim.evap_gibe; 
evap_turkana         = sim.evap_turkana; 
cyclo_inflow_delta   = sim.cyclo_inflow_delta;
cyclo_inflow_gibe    = sim.cyclo_inflow_gibe;
date_day             = sim.date_day;


%% set r_g3 as fraction k of the cyclostationary inflow in Gibe III
j = 1;

for i = 1:length(SPEI)
    idx = find(spei_date(i) == date_day);
    k_fract(j:idx) =  k(SPEI(i));
    j = idx + 1;
end


r_g3 = cyclo_inflow_gibe(1:length(k_fract)).*k_fract;


%% G3 simulation

delta=60*60*24; % [seconds/day]
r_g3 = [zeros(1, lag+1), r_g3];

% pre-allocate vectors
v_g3 = nan(size(r_g3)); 
v_T  = nan(size(r_g3));
h_g3 = nan(size(r_g3));
inflowDelta = nan(size(r_g3));

% initial conditions
v_g3(lag) = dead_storage_g3;
h_g3(lag) = interp1qr(lsv_gibeIII(3,:), lsv_gibeIII(1,:), v_g3(lag ));
v_T(lag)  = interp1qr(lsv_Turkana(1,:), lsv_Turkana(3,:), 365);

i = lag;
while ( (i < length(r_g3)) + (h_g3(i)<norm_operating_level) )==2

    r = interp1qr(min_rel(1,:), min_rel(2,:), h_g3(i));
    R = interp1qr(max_rel(1,:), max_rel(2,:), h_g3(i));
    r_g3(i+1)=max( min(r_g3(i+1), R),r);
    
    if v_g3(i)<=dead_storage_g3
        r_g3(i+1) = 0; 
    end
    
    v_g3(i+1) = v_g3(i) + (inflow_gibe(i+1) - r_g3(i+1) - evap_gibe(i+1))*delta;
    h_g3(i+1) = interp1qr(lsv_gibeIII(3,:), lsv_gibeIII(1,:), v_g3(i+1));
    
    inflowDelta(i+1) = r_g3(i+1 - lag) + inflow_turkana(i + 1) ;
    v_T(i+1) = v_T(i) + (inflowDelta(i+1) - evap_turkana(i+1))*delta ;
    i = i+1;
    
end

% once operating level is reached, a regime operating policy is
% implemented, designed via SDP.

if h_g3(i)>=norm_operating_level
    h_init = h_g3(i);
    q_sim = inflow_gibe(i+1:end);
    s_init = interp1qr(lsv_gibeIII(1,:), lsv_gibeIII(3,:), h_init) ;
    [~,h,r] = simulateSDPGibe3( q_sim, s_init, i );
    h_g3(i+1:end) = h(2:end);
    r_g3(i+1:end) = r(2:end);
    inflowDelta(i+1:end) = r_g3(i+1 - lag:end-lag) + inflow_turkana(i + 1:end)  ;
end


%% remove lag times in variables needed to compute objectives
h_g3 = h_g3(lag +1 : end);
inflowDelta = inflowDelta(lag +1 : end);
r_g3 = r_g3(lag +1 : end);
Ny = 4;

%% hydropower
qmaxGIII = 10*95; 
q_turb=min(r_g3, qmaxGIII);
g_hyd_GIII = effic.*9.81.*q_turb*24.*(h_g3-g3_tail)./10^3;  
Ghyd = - sum(g_hyd_GIII )/Ny/1000;

%% environment
g_env     = (cyclo_inflow_delta - inflowDelta).^2; 
Genv      = sum(g_env)/length(inflowDelta)/10000;

%% additional objective to reach operating level h = 201
Glevel = h_g3(end); 
if Glevel < min_operating_level
    Glevel = -Inf; 
end

J = [Ghyd, Genv, -Glevel];
end

