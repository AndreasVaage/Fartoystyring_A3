function [pos_x, pos_y, vel_x, vel_y] = target_motion(WP,U_t,t)

    psiTemp=atan2(WP(2,2)-WP(2,1),WP(1,2)-WP(1,1));
    vel_x = U_t*cos(psiTemp);
    vel_y = U_t*sin(psiTemp);
    pos_x = WP(1,2) + vel_x*t;
    pos_y = WP(2,2) + vel_y*t;
    
end
