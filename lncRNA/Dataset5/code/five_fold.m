%tju cs for bioinformatics 
clear;
seed = 12345678;
rand('seed', seed);
load('../data/disease_sim_2017.mat');
load('../data/lncR_disease_2017.mat');
load('../data/lncR_sim_2017.mat');
lncSim = lncR_sim_matrix;
disSim_Jaccard = disease_sim_matrix;
interMatrix = lncR_disease_matrix;
y = interMatrix;

K1 = [];
mat = process_kernel(lncSim);
K1(:,:,1)=Knormalized(mat);
K2 = [];
mat = process_kernel(disSim_Jaccard);
K2(:,:,1)=Knormalized(mat);


nfolds =5; 
k_nn = 10;IsHG =1;
%k1 = 100;k2=98; %avg,grtmf
k1 = 100;k2 = 400;
lambda1 = 2^0;lambda2=2^0;
crossval_idx = crossvalind('Kfold',y(:),nfolds);

for fold = 1:nfolds
    y_train = interMatrix;
	test_idx  = find(crossval_idx==fold);
	y_train(test_idx) = 0;
	%K1(:,:,2)=kernel_corr(y_train,1,0,1);
	%K2(:,:,2)=kernel_corr(y_train,2,0,1);
	[KD, KL] = GaussianKernel(y_train', 1, 1)
	K1(:,:,2)=KL;
	K2(:,:,2)=KD;
	[KD, KL] = consine(y_train');
	K1(:,:,3)=KL;
	K2(:,:,3)=KD;
	
	%K1 = K1(:,:,3);
	%K2 = K2(:,:,3);
		
		
	% 2. multiple kernel 
	weight_v1 = randn(size(K1,3),1);
	weight_v1(1:size(K1,3)) = 1/size(K1,3);
	weight_v2 = weight_v1;
	
	
	[weight_v1] = cka_kernels_weights(K1,y_train,1);
	%weight_v1 = ones(size(K1,3),1);weight_v1 = ones(size(K1,3),1)/size(K1,3);
	K_COM1 = combine_kernels(weight_v1, K1);
	
	
	[weight_v2] = cka_kernels_weights(K2,y_train,2);
	%weight_v2 = ones(size(K2,3),1);weight_v2 = ones(size(K2,3),1)/size(K2,3);
	K_COM2 = combine_kernels(weight_v2, K2);
	
	
	%3.GRTMF
	[F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);
	%[F_1] = cmf(K_COM1,K_COM2,y_train,k1,k2,lambda1,lambda2,0.1);
	%[F_1] = grmf(K_COM1,K_COM2,y_train,k1,k2,0.1,0.001,0.001,10);
	%[F_1] = tmf(K_COM1,K_COM2,y_train,1,50,50);
	%[F_1] = nrlmf(K_COM1,K_COM2,y_train,k1,k2,0.125,0.125,3,0.0001,0.0001,0.001,10);
    y(test_idx)= F_1(test_idx);
end

[X_1,Y_1,tpr,aupr_F_1] = perfcurve(interMatrix(:),y(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_F_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(interMatrix(:),y(:),1);

