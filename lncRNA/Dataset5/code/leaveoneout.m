%tju cs for bioinformatics 
clear
load('../data/disease_sim_2017.mat');
load('../data/lncR_disease_2017.mat');
load('../data/lncR_sim_2017.mat');
lncSim = lncR_sim_matrix;
disSim_Jaccard = disease_sim_matrix;
interMatrix = lncR_disease_matrix;
k_nn = 10;IsHG =1;
%k1 = 250;k2 = 400;
k1 = 100;k2=400;
lambda1 = 2^0;lambda2=2^0;
y_train = interMatrix;
[II,JJ] = find(interMatrix == 1);
K1 = [];
K1(:,:,1)=lncSim;
K2 = [];
K2(:,:,1)=disSim_Jaccard;
Pre_value =interMatrix;
weight=[];
for i = 1:length(II)    
    y_train(II(i),JJ(i)) = 0;
    [KD, KL] = GaussianKernel(y_train', 1, 1)
	K1(:,:,2)=KL;
	K2(:,:,2)=KD;
	[KD, KL] = consine(y_train');
	K1(:,:,3)=KL;
	K2(:,:,3)=KD;
	

    [weight_v1] = cka_kernels_weights(K1,y_train,1);
	%weight_v1 = ones(size(K1,3),1);weight_v1 = ones(size(K1,3),1)/size(K1,3);
    K_COM1 = combine_kernels(weight_v1, K1);		
 
    [weight_v2] = cka_kernels_weights(K2,y_train,2);
	%weight_v2 = ones(size(K2,3),1);weight_v2 = ones(size(K2,3),1)/size(K2,3);
    K_COM2 = combine_kernels(weight_v2, K2);
     
    [F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);
   
   
    Pre_value(II(i),JJ(i)) = F_1(II(i),JJ(i));
    y_train = interMatrix;
    if mod(i,10) == 0
        i      
    end
end


[KD, KL] = GaussianKernel(y_train', 1, 1)
K1(:,:,2)=KL;
K2(:,:,2)=KD;
[KD, KL] = consine(y_train');
K1(:,:,3)=KL;
K2(:,:,3)=KD;


[weight_v1] = cka_kernels_weights(K1,y_train,1);
%weight_v1 = ones(size(K1,3),1);weight_v1 = ones(size(K1,3),1)/size(K1,3);
K_COM1 = combine_kernels(weight_v1, K1);		
 
[weight_v2] = cka_kernels_weights(K2,y_train,2);
%weight_v2 = ones(size(K2,3),1);weight_v2 = ones(size(K2,3),1)/size(K2,3);
K_COM2 = combine_kernels(weight_v2, K2);
     
[F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);

for i =1:length(F_1(:))
    if interMatrix(i)==0
        Pre_value(i)=F_1(i);
    end
end

[X_1,Y_1,tpr,aupr_1] = perfcurve(interMatrix(:), Pre_value(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(interMatrix(:), Pre_value(:),1);

