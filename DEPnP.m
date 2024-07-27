function [q_DEPnP,T_DEPNP,depth_err] = DEPnP(Xw,impts,Xc_true)
%% new equation
n = size(Xw,2);
pts_mean = mean(Xw,2);
pts_uncent = Xw - pts_mean;
if n < 6
    Cw=define_control_points();
    Alpha=compute_alphas(Xw',Cw);
    M = kron(Alpha,[1 0 -1; 0 1 -1]);
    M(:,[3,6,9,12]) =  M(:,[3,6,9,12]) .* impts(:);
else
    Cw = DEPnP_define_control_points_pca(pts_mean,pts_uncent);
    Alpha=compute_alphas(Xw',Cw);
    M = kron(Alpha,[1 0 -1; 0 1 -1]);
    M(:,[3,6,9,12]) =  M(:,[3,6,9,12]) .* impts(:);
end

%initial
[Cc0s,~]=eig(M'*M);
[Cc0,Xc0,~] = DEPnP_compute_norm_sign_scaling_factor(Cc0s(:,1),Alpha,pts_uncent');
Alpha = Alpha';Cw = Cw';
%compute Xc0 and q0
% METHOD 1
% A0 = Xc0 - mean(Xc0,2);
% B0 = Xw - mean(Xw,2);
% CB0 = Cc0 - mean(Cc0,2);
% CA0 = Cw - mean(Cw,2);
% W0 = computeW(B0,A0);
% Wc0 = computeW(CB0,CA0);
% [q0, Sw0, Vw0] = DEPnP_kernel_noise([W0;Wc0],1);
% q0 = q0/norm(q0);
% METHOD 2
R0 = DEPnP_getRot(pts_uncent',Xc0);
% q0 = (rotm2quat(R0))';
q0 = DEPnP_rotm2quat(R0);

x0 = [Cc0(:,3);q0];
%compute CB
CB = Cw - mean(Cw,2);
%optimize
Alphr =  Alpha'/(Alpha*Alpha');
imptPart = compteJacobian_imptsPart([impts;ones(1,n)],Alpha,Alphr);
[x_opt_gn,~,~,depth_Xc] = DEPnP_gauss_newton(CB,Alpha,impts,Alphr,imptPart,x0);
%calculate q t
q_DEPnP = x_opt_gn(5:8)';
R_DEPnP = DEPnP_quat2Rotm(q_DEPnP);
T_DEPnP = mean(depth_Xc,2) - R_DEPnP*pts_mean;

Cc_true = Xc_true*Alphr;
depth_err = norm(x_opt_gn(1:4)'-Cc_true(end,:))/norm(x_opt_gn(1:4))*100;

end


