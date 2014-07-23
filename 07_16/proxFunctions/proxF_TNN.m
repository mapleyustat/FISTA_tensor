function [x, objV] = proxF_TNN(y,rho,Tsize)

Y  = reshape(y,Tsize);
n1 = Tsize(1);
n2 = Tsize(2);
n3 = Tsize(3);

% this returned this STILL FFT'd version
[U S V]=ntsvd(Y,1);

sTrueV =zeros(min(n1,n2),n3);
for i = 1:n3
    s = S(:,:,i);  % tube
    s = diag(s);
    sTrueV(:,i) = s;
end
%%
[sTrueV objV] = proxF_l1(sTrueV,rho);

%%
for i = 1:min(n1,n2) 
    for j = 1:n3
        S(i,i,j) = sTrueV(i,j);
    end
end
    
U = ifft(U,[],3);
S = ifft(S,[],3);
V = ifft(V,[],3);

% X = tprod( tprod(U,S), ttrans_HO(V));

X = tprod( tprod(U,S), tran(V));

x = X(:);
