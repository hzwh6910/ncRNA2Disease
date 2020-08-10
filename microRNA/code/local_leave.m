%tju cs for bioinformatics 
clear;
load('../data/dataset.mat');
y_train = miRNA_disease_Y;
[II,JJ] = find(miRNA_disease_Y == 1);
K1 = [];
K1(:,:,1)=miRNA_Function_S;
K1(:,:,2)=miRNA_Sequences_Needle_S;
K2 = [];
K2(:,:,1)=disease_Function_S;
K2(:,:,2)=disease_Sem_S;
[n_miRNA,n_Diseas]=size(y_train);
Pre_value = [];
k_nn = 10;IsHG =0;
%k1 = 250;k2 = 400;
k1 = 300;k2=1100;
lambda1 = 2^0;lambda2=2^0;
for i = 1:n_Diseas  
    y_train = miRNA_disease_Y;
    y_train(:,i) = 0;
	
	K1(:,:,3)=kernel_corr(y_train,1,0,1);
	K2(:,:,3)=kernel_corr(y_train,2,0,1);

	% 2. multiple kernel 
	weight_v1 = randn(size(K1,3),1);
	weight_v1(1:size(K1,3)) = 1/size(K1,3);
	weight_v2 = weight_v1;
	
	
	[weight_v1] = cka_kernels_weights(K1,y_train,1);
	K_COM1 = combine_kernels(weight_v1, K1);
	
	
	[weight_v2] = cka_kernels_weights(K2,y_train,2);
	K_COM2 = combine_kernels(weight_v2, K2);
	
	
	%3.GRTMF
	[F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);
    Pre_value = [Pre_value,F_1(:,i)];
     if mod(i,50) == 0
         i             
     end
end

[X_1,Y_1,tpr,aupr_1] = perfcurve(miRNA_disease_Y(:), Pre_value(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(miRNA_disease_Y(:), Pre_value(:),1);


