function quat = depnp_rotm2quat(R)
    quat = zeros(4,1);
    quat(1) = 0.5*sqrt(trace(R)+1);
    quat_quart_rev = 1/(4*quat(1));
    quat(2) = (R(3,2) - R(2,3))*quat_quart_rev;
    quat(3) = (R(1,3) - R(3,1))*quat_quart_rev;
    quat(4) = (R(2,1) - R(1,2))*quat_quart_rev;
end