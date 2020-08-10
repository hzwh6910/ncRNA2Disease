function [L_H] = construct_Hypergraphs_knn(W,k_nn)

L_H=[];
n_size = size(W,1);
n_vertex = n_size;
n_edge = n_size;


H = zeros(n_vertex,n_edge);
%build Association matrix of Hypergraphs
	for i=1:n_vertex
		ll = W(i,:);
		ll(i)=[];
		[B,index_i] = sort(ll,'descend');
		k_ii = index_i(1:k_nn);
		H(i,k_ii) = 1;
	end

d_v = sum(H');
d_e = sum(H);

D_v = diag(d_v);
D_e_I = diag(d_e.^(-1));
D_e_I(D_e_I==inf)=0;

A = eye(n_edge);

I = eye(n_size);


Thta = A*D_e_I;
Thta = H*Thta*H';
Thta = (D_v^(-0.5))*Thta*(D_v^(-0.5));

L_H = I - Thta;
end