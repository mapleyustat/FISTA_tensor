function X_bar = t_blkdiag( X_hat )

[ ~ , ~ , n3] = size(X_hat)            ;
X_bar         = X_hat(:,:,1)           ;
for i = 2:n3
    X_bar = blkdiag(X_bar,X_hat(:,:,i));
end

end

