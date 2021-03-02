function [mr, MR] = min_max_relG3(h, doy, sys_param)

% min release
if h < 140
    mr = 0.0;
elseif h >= 232
    mr = interp_lin_scalar( sys_param.maxRelG3(1,:), sys_param.maxRelG3(2,:), h);
else
    mr = sys_param.MEF(doy);
end

% max release
MR = interp_lin_scalar( sys_param.maxRelG3(1,:), sys_param.maxRelG3(2,:), h);

end