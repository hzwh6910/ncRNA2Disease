function [LapA] = grtmf(W1,W2,y,lamda1,lamda2,k1,k2,k_nn,IsHG)
LapA=[];
KK=[];
[n,m] = size(y);
%S1=preprocess_PNN(W1,p_nearest_neighbor);
%S2=preprocess_PNN(W2,p_nearest_neighbor);
L_H1=[];
L_H2=[];
if IsHG ==1
	%S_1 =preprocess_PNN(W1,k_nn);
	%S_2 =preprocess_PNN(W2,k_nn);
	[L_H1] = construct_Hypergraphs_knn(W1,k_nn);
	[L_H2] = construct_Hypergraphs_knn(W2,k_nn);
else

	S_1 =preprocess_PNN(W1,k_nn);
	d_1 = sum(S_1);
	D_1 = diag(d_1);
	L_D_1 = D_1 - S_1;
	d_tmep_1=eye(n)/(D_1^(1/2));
	L_H1 = d_tmep_1*L_D_1*d_tmep_1;
	
	S_2 =preprocess_PNN(W2,k_nn);
	d_2 = sum(S_2);
	D_2 = diag(d_2);
	L_D_2 = D_2 - S_2;
	d_tmep_2=eye(m)/(D_2^(1/2));
	L_H2 = d_tmep_2*L_D_2*d_tmep_2;

end

[U1,S_k,V1] = svds(W1,k1);
G1 = U1*(S_k^0.5);  


[U2,S_k,V2] = svds(W2,k2);
G2 = U2*(S_k^0.5); 

A = G1;
B = G2;

inv_A = pinv(A);
t_inv_B = pinv(B');
ki = eye(n);
%A = W1;
%B = W2;
	a = inv_A*(ki + lamda1*L_H1)*A;
	b = B'*L_H2*t_inv_B;  b = lamda2*b;
	c = inv_A*y*t_inv_B;
	KK = sylvester(a,b,c);
	
	%reconstruct Y*
LapA = A*KK*B';

end


function S=preprocess_PNN(S,p)
%preprocess_PNN sparsifies the similarity matrix S by keeping, for each
%drug/target, the p nearest neighbors and discarding the rest.
%
% S = preprocess_PNN(S,p)

    NN_mat = zeros(size(S));

    % for each drug/target...
    for j=1:length(NN_mat)
        row = S(j,:);                           % get row corresponding to current drug/target
        row(j) = 0;                             % ignore self-similarity
        [~,indx] = sort(row,'descend');         % sort similarities descendingly
        indx = indx(1:p);                       % keep p NNs
        NN_mat(j,indx) = S(j,indx);             % keep similarities to p NNs
        NN_mat(j,j) = S(j,j);                   % also keep the self-similarity (typically 1)
    end

    % symmetrize the modified similarity matrix
    S = (NN_mat+NN_mat')/2;

end