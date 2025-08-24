function n = nTau(data,min,max)

n = 0;
for i=1:length(data)
    val = data(i);
    if  val >= min && val < max 
        n = n+1;
    end 

end

n = round(n);