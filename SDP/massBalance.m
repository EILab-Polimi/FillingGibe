function [s1,r1] = massBalance( s, u, q, doy )

% Output:
%      s1 - final storage. 
%      r1 - release over the 24 hours. 
%
% Input:
%       s - initial storage. 
%       u - release decision 
%       q - inflow 
%       dt - day of the year
%
% See also MIN_RELEASE, MAX_RELEASE

global sys_param lsv_gibeIII;

HH = sys_param.integration_substep;
delta = sys_param.simulation.delta/HH;
s_ = nan(HH+1,1);
h_ = nan(HH+1,1);
r_ = nan(HH+1,1);
A_ = nan(HH+1,1);

eG3 = sys_param.simulation.ev  ;
lsv_gibeIII = sys_param.lsv;
s_(1) = s;
h_(1) = interp_lin_scalar(lsv_gibeIII(3,:), lsv_gibeIII(1,:), s_(1));
A_(1) = interp_lin_scalar(lsv_gibeIII(3,:), lsv_gibeIII(2,:), s_(1));
for i=1:HH
    [r_min, r_max] = min_max_relG3( h_(i), doy, sys_param ) ;
    r_(i+1) = min( r_max , max( r_min , u ) );
    s_(i+1) = s_(i) + delta*( q - r_(i+1) ) - eG3(doy)*A_(i)/1000/HH ;
    h_(i+1) = interp_lin_scalar(lsv_gibeIII(3,:), lsv_gibeIII(1,:), s_(i+1));
    A_(i+1) = interp_lin_scalar(lsv_gibeIII(3,:), lsv_gibeIII(2,:), s_(i+1));
end

s1 = s_(HH);
r1 = mean(r_(2:end));
