function f = alpha(d,lambda,e,n,nnew, t, tnew)

val = zeros(d,1);

for i=1:d % e måste vara d+1 lång 
   if i < d
           interval = e(i+1) - e(i); % Compute interval width
   else
            interval = 1; 
   end 
   HL1 = -lambda(i)*interval;
   HL2 = log(t(i+1)-t(i));
   HL3 = -log(tnew(i+1)-tnew(i)); 
   HL4 = (nnew(i)-n(i))*log(lambda(i));
   HL = HL1+HL2+HL3+HL4;
   val(i) = HL;
end

f = exp(sum(val));
f = min(f,1);



