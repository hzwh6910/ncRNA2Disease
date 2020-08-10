%tju cs for bioinformatics 
clear;
load('../data/disSim_Jaccard.mat');
load('../data/interMatrix.mat');
load('../data/lncSim.mat');
k_nn = 10;IsHG =1;
%k1 = 250;k2 = 400;
k1 = 200;k2=600;
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
    K1(:,:,2)=kernel_corr(y_train,1,0,1);
	K2(:,:,2)=kernel_corr(y_train,2,0,1);
	

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


K1(:,:,2)=kernel_corr(y_train,1,0,1);
K2(:,:,2)=kernel_corr(y_train,2,0,1);


[weight_v1] = cka_kernels_weights(K1,y_train,1);
K_COM1 = combine_kernels(weight_v1, K1);		
 
[weight_v2] = cka_kernels_weights(K2,y_train,2);
K_COM2 = combine_kernels(weight_v2, K2);
     
[F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);

for i =1:length(F_1(:))
    if interMatrix(i)==0
        Pre_value(i)=F_1(i);
    end
end

[X_1,Y_1,tpr,aupr_1] = perfcurve(interMatrix(:), Pre_value(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(interMatrix(:), Pre_value(:),1);

