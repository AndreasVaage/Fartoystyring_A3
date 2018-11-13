function [U_d,X_d] = constant_bearing_guidance(p,pos_t,vel_t)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
delta = 1000;
U_amax = 10;
p_err = p - pos_t;
dist = 400; % [m]
p_err_dist = max(0,norm(p_err) - dist);

k = U_amax* p_err_dist/sqrt(p_err_dist^2 + delta^2);

v_a = -k*p_err/norm(p_err);
v_d = vel_t + v_a;

U_d = norm(v_d);
X_d = atan2(v_d(2),v_d(1));
end

