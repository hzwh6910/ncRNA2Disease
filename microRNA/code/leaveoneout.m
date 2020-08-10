%tju cs for bioinformatics 
clear
load('../data/dataset.mat');
k_nn = 10;IsHG =1;
%k1 = 250;k2 = 400;
k1 = 300;k2=1100;
lambda1 = 2^0;lambda2=2^0;
y_train = miRNA_disease_Y;
[II,JJ] = find(miRNA_disease_Y == 1);
K1 = [];
K1(:,:,1)=miRNA_Function_S;
K1(:,:,2)=miRNA_Sequences_Needle_S;
K2 = [];
K2(:,:,1)=disease_Function_S;
K2(:,:,2)=disease_Sem_S;
Pre_value =miRNA_disease_Y;
weight=[];
for i = 1:length(II)    
    y_train(II(i),JJ(i)) = 0;
    K1(:,:,3)=kernel_corr(y_train,1,0,1);
	K2(:,:,3)=kernel_corr(y_train,2,0,1);

    [weight_v1] = cka_kernels_weights(K1,y_train,1);
	%weight_v1 = ones(size(K1,3),1);weight_v1 = ones(size(K1,3),1)/size(K1,3);
    K_COM1 = combine_kernels(weight_v1, K1);		
 
    [weight_v2] = cka_kernels_weights(K2,y_train,2);
	%weight_v2 = ones(size(K2,3),1);weight_v2 = ones(size(K2,3),1)/size(K2,3);
    K_COM2 = combine_kernels(weight_v2, K2);
     
    [F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);
   
    Pre_value(II(i),JJ(i)) = F_1(II(i),JJ(i));
    y_train = miRNA_disease_Y;
    if mod(i,500) == 0
        i      
    end
end


K1(:,:,3)=kernel_corr(y_train,1,0,1);
K2(:,:,3)=kernel_corr(y_train,2,0,1);


[weight_v1] = cka_kernels_weights(K1,y_train,1);
K_COM1 = combine_kernels(weight_v1, K1);		
 
[weight_v2] = cka_kernels_weights(K2,y_train,2);
K_COM2 = combine_kernels(weight_v2, K2);
     
[F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);

for i =1:length(F_1(:))
    if miRNA_disease_Y(i)==0
        Pre_value(i)=F_1(i);
    end
end

[X_1,Y_1,tpr,aupr_1] = perfcurve(miRNA_disease_Y(:), Pre_value(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(miRNA_disease_Y(:), Pre_value(:),1);

