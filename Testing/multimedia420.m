n = 6;
a =reshape(1:n^2,n,[])';

x = [];

for j = 1:size(a,1)
    for k = 1:size(a,2)
        indHor = mod(ceil((size(a,2)*(j-1)+k)/2)+1,2)+1+2*floor((j-1)/2);
        indVer = mod(k-1,2)+1 + mod(2*floor((size(a,2)*(j-1)+k-1)/4),size(a,2));
        x(indHor,indVer) = a(j,k);
        
        fprintf('(%d,%d) --> (%d,%d) \n',j,k,indHor,indVer)
    end
end
a
x
