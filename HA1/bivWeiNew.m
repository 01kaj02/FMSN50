function x = bivWeiNew(v1, v2, k, lambda)

alpha = 0.638;
p = 3;
q = 1.5;
f1 = density(v1,k,lambda);
f2 = density(v2,k,lambda);
WF1 = WeiF(v1, k, lambda);
WF2 = WeiF(v2, k, lambda);
x = f1 .* f2 .* (1 + alpha .* (1 - WF1.^p).^(q-1) .* (1 - WF2.^p).^(q-1).*(WF1.^p .* (1 + p .* q) - 1) .* (WF2.^p .* (1 + p .* q) - 1));

end 