function X = t_iblkdiag( X_bar , tensor_size )
n1 = tensor_size(1);
n2 = tensor_size(2);
n3 = tensor_size(3);
X = zeros(tensor_size)  ;
for i = 1:n3
    X(:,:,i) = X_bar(n1*(i-1)+1:n1*i,...
        n2*(i-1)+1:n2*i) ;
end


end

