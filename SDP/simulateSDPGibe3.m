function [ s_GIII, h_GIII, r_GIII] = simulateSDPGibe3( q, s_in, i )
    
global policy sys_param;

sys_param.integration_substep=24;
% day of the year
q_sim = [ nan; q' ];
H = length(q_sim) - 1;
T = 365;
% vectors initialization 
s_GIII = nan(length(q_sim),1);
h_GIII = nan(length(q_sim),1);
r_GIII = nan(length(q_sim),1);
% GibeIII/Turkana characteristics
lsv_g3 = sys_param.lsvG3;
% initial conditions
s_GIII(1) = s_in ;
h_GIII(1) = interp1qr(lsv_g3(3,:), lsv_g3(1,:), s_GIII(1));


% simulation loop
for t = 1:H
    
    doy = mod(i+t-1,T)+1 ;
    
    % release decision
      discr_s = sys_param.algorithm.discr_s;
      discr_q = sys_param.algorithm.discr_q;
      discr_u = sys_param.algorithm.discr_u;
      
      min_rel = sys_param.algorithm.min_rel;
      max_rel = sys_param.algorithm.max_rel;
      
      [ ~ , idx_q ] = min( abs( discr_q - q_sim(t+1) ) );
      
      % Minimum and maximum release for current storage and inflow:
    
    v = interp_lin_scalar( discr_s , min_rel( : , idx_q, doy ) , s_GIII(t) );
    sys_param.simulation.vv = repmat( v, 1, length(discr_q) );
    
    V =  interp_lin_scalar( discr_s , max_rel( : , idx_q, doy ) , s_GIII(t) );
    sys_param.simulation.VV = repmat( V, 1, length(discr_q) );
    [ ~, idx_u ] = Bellman_sdp( policy.H(:,doy) , s_GIII(t), doy );
    % Choose one decision value (idx_u can return multiple equivalent decisions)
    uu = extractor_ref( idx_u , discr_u );
    
    [s_GIII(t+1), r_GIII(t+1)] = massBalance( s_GIII(t), uu, q_sim(t+1), doy );
    h_GIII(t+1) = interp_lin_scalar(lsv_g3(3,:), lsv_g3(1,:), s_GIII(t+1));
end



end


