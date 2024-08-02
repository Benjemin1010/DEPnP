%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=depnp_getRot(pts_uncent,cpts)

ccent=mean(cpts);
% wcent=mean(wpts);

A = cpts - ccent;
% B = wpts - wcent;

M = A'*pts_uncent;

[U S V]=svd(M);
R=U*V';

if det(R)<0
  R=-R;
end

% T=ccent'-R*wcent';

end