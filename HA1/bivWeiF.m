function x = bivWeiF(v1, v2, k, lambda)

alpha = 0.638;
p = 3;
q = 1.5;
WF1 = WeiF(v1, k, lambda);
WF2 = WeiF(v2, k, lambda);
x = WF1.*WF2.*(1+alpha.*(1-WF1.^p).^(q).*(1-WF2.^p).^(q));

end 