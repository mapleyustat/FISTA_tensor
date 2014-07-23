function objV = TNN(X)
% calculate the tensor nuclear norm
[~, S ,~] = ntsvd(X,0);
objV    = sum(sum(sum(S)));