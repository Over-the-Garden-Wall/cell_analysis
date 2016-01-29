function [bestcut maxval] = max_cut(W, iterations)

if ~exist('iterations','var') || isempty(iterations)
    iterations = 100;
end

try
    cvx_check_dimension(1);
catch    
    run ~/cvx/cvx_setup.m
end


n=size(W,1);
% L = 1/4 * (diag(W*ones(n,1))-W);

  
disp('starting CVX')
% Construct and solve the model
cvx_begin sdp
    variable X(n,n) symmetric
    minimize( trace(  W*X  ) )
    subject to 
        diag(X) == ones(n,1);
    	X >= 0;
cvx_end
% 
% [V D] = eig(Y);
% % Look at the rank
% d = diag(D);
% plot(diag(D))
% % This matrix has rank that's effectively 5
% r = 5;
% X = V(:,end-(r-1):end)*diag(sqrt(d(end-(r-1):end)));

[V, D] = eig(X);

maxval = 0;
bestcut = [];
for t = 1:iterations
    mcut = sign((V * diag(sqrt(diag(D)))) * randn(n,1));

    cutval = .5*sum(sum(W.*(1-mcut*mcut')));
    
    if cutval > maxval
        maxval = cutval;
        bestcut = mcut;
    end
end

end