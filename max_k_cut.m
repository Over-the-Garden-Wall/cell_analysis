function [bestcut, cutval] = max_k_cut(B, K, iterations)


if ~exist('iterations','var') || isempty(iterations)
    iterations = 100;
end

try 
    mytime(1);
catch    
    run ~/cvx/sdpt3/install_sdpt3.m
end

% [blk,Avec,C,b,objval,X] = max_kcut(W,k,1)

% Copyright (c) 1997 by
% K.C. Toh, M.J. Todd, R.H. Tutuncu
% Last modified: 18 May 07
if ~isreal(B); error('only real B allowed'); end;
n = length(B); e = ones(n,1);
n2 = n*(n-1)/2;
C{1} = -(1-1/K)/2*(spdiags(B*e,0,n,n)-B);
b = e;
blk{1,1} = 's'; blk{1,2} = n;
blk{2,1} = 'l'; blk{2,2} = n2;
A = cell(1,n);
for j = 1:n; A{j} = sparse(j,j,1,n,n); end;
Avec = svec(blk(1,:),A,1);
tmp = speye(n*(n+1)/2);
idx = cumsum([1:n]);
Atmp = tmp(:,setdiff([1:n*(n+1)/2],idx));
Avec{1,1} = [Avec{1,1}, Atmp/sqrt(2)];
Avec{2,1} = [sparse(n2,n), -speye(n2,n2)];
b = [b; -1/(K-1)*ones(n2,1)];
C{2,1} = zeros(n2,1);
%
[obj,X,y,Z] = sqlp(blk,Avec,C,b);
X = X{1};

maxval = 0;
for t = 1:iterations
    z = randn(n,K);
    [dummy, mcut] = max(X * z, [], 2);

    cutval = 0;
    for i = 1:K
        cutval = cutval + sum(sum(B.* (double(mcut==i)*double(mcut'~=i))));
    end
    
    if cutval > maxval
        maxval = cutval;
        bestcut = mcut;
    end
end
    
end