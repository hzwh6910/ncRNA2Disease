%tju cs for bioinformatics 
clear;
seed = 12345678;
rand('seed', seed);
lncRNAsimilarity = importdata('../data/lncRNAsimilarity.txt');
diseasesimilarity=importdata('../data/diseasesimilarity.txt');
geneDis = load('../data/known_gene_disease_interaction.txt');
interMatrix = load('../data/known_lncRNA_disease_interaction.txt');
lncSim = lncRNAsimilarity;
disSim_Jaccard = diseasesimilarity;
y = interMatrix;

K1 = [];
mat = process_kernel(lncSim); 
K1(:,:,1)=Knormalized(mat);
K1(:,:,2) = Knormalized(process_kernel(ncRNASS(interMatrix, disSim_Jaccard)));
K2 = [];
mat = process_kernel(disSim_Jaccard);
K2(:,:,1)=Knormalized(mat);
K2(:,:,2)= Knormalized(process_kernel(kernel_corr(geneDis,2,0,1)));


nfolds =5; 
k_nn = 10;IsHG =1;
k1 = 200;k2 = 300; 
lambda1 = 2^0;lambda2=2^0;
crossval_idx = crossvalind('Kfold',y(:),nfolds);

for fold = 1:nfolds
    y_train = interMatrix;
	test_idx  = find(crossval_idx==fold);
	y_train(test_idx) = 0;
	K1(:,:,3) = kernel_corr(y_train,1,0,1)
	K2(:,:,3) = kernel_corr(y_train,2,0,1);
	
	%K1=K1(:,:,2);
	%K2=K2(:,:,2);
	
	
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
	%[F_1] = cmf(K_COM1,K_COM2,y_train,115,k2,lambda1,lambda2,0.1);
	%[F_1] = grmf(K_COM1,K_COM2,y_train,115,k2,0.1,0.001,0.001,10);
	%[F_1] = tmf(K_COM1,K_COM2,y_train,1,100,100);
	%[F_1] = nrlmf(K_COM1,K_COM2,y_train,k1,k2,0.125,0.125,3,0.0001,0.0001,0.001,10);
    y(test_idx)= F_1(test_idx);
end

[X_1,Y_1,tpr,aupr_F_1] = perfcurve(interMatrix(:),y(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_F_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(interMatrix(:),y(:),1);

