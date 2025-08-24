function x = density(v, k, lambda)

x = (k/lambda).* (v/lambda).^(k-1).*exp(-(v/lambda).^k);
end 