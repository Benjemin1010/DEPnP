function R = depnp_quat2Rotm(q)

                       b2 = q(2)*q(2);c2 = q(3)*q(3);d2 = q(4)*q(4);
        ab = q(1)*q(2)depnp;ac = q(1)*q(3);ad = q(1)*q(4);
        bc = q(2)*q(3);bd = q(2)*q(4);
        cd = q(3)*q(4);

        R = [
                 [1-2*(c2+d2)  2*(bc-ad)  2*(ac+bd)];
                 [2*(bc+ad)  1-2*(b2+d2)  2*(cd-ab)];
                 [2*(bd-ac)  2*(ab+cd)  1-2*(b2+c2)];
            ];
end