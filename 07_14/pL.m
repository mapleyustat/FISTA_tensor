function [X TNN] = pL( Y , lambda )

[U S V] = svd(Y);
X       = U* (S .* (S>lambda)) * V' ;
TNN     = sum(sum(S .* (S>lambda))) ;

end

