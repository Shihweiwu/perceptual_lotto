function y = rand_vm(mu,k,n)
%
% 
% mu: mean
% k: variance
% n: sample size
%
% Shih-Wei Wu.


x = 1:1:360;
mu
mPdfs = vmPdfs(x,mu,k,'norm');
cum_p = cumsum(mPdfs);

rand_unif = rand(n,1);

y = zeros(n,1);
for i=1:n
    index=find(cum_p<=rand_unif(i));
    if sum(index)==0
        rand_index = 1;
    else
        rand_index = max(index)+1;
    end
    
    y(i) = x(rand_index);
    
end

figure(1);clf
hist(y,50);
title(['von Mises, mu=' num2str(mu) ' k=' num2str(k) ' sample size=' num2str(n)]);


