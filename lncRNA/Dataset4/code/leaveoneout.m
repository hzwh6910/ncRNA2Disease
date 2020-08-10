%tju cs for bioinformatics 
clear
lncRNAsimilarity = importdata('../data/lncRNAsimilarity.txt');
diseasesimilarity=importdata('../data/diseasesimilarity.txt');
geneDis = load('../data/known_gene_disease_interaction.txt');
interMatrix = load('../data/known_lncRNA_disease_interaction.txt');
lncSim = lncRNAsimilarity;
disSim_Jaccard = diseasesimilarity;
k_nn = 10;IsHG =1;
%k1 = 250;k2 = 400;
k1 = 600;k2 = 800;
lambda1 = 2^0;lambda2=2^0;
y_train = interMatrix;
[II,JJ] = find(interMatrix == 1);
index = find(interMatrix == 1);
K1 = [];
mat = process_kernel(lncSim); 
K1(:,:,1)=Knormalized(mat);
K1(:,:,2) = ncRNASS(interMatrix, disSim_Jaccard);
K2 = [];
mat = process_kernel(disSim_Jaccard);
K2(:,:,1)=Knormalized(mat);
K2(:,:,2)= kernel_corr(geneDis,2,0,1);
Pre_value =interMatrix;
weight=[];
for i = 1:length(II)    
    y_train(II(i),JJ(i)) = 0;
    K1(:,:,3) = kernel_corr(y_train,1,0,1)
	K2(:,:,3) = kernel_corr(y_train,2,0,1);
    %K2(:,:,2) = kernel_corr(y_train,2,0,1);

    [weight_v1] = cka_kernels_weights(K1,y_train,1);
	%weight_v1 = ones(size(K1,3),1);weight_v1 = ones(size(K1,3),1)/size(K1,3);
    K_COM1 = combine_kernels(weight_v1, K1);		
 
    [weight_v2] = cka_kernels_weights(K2,y_train,2);
	%weight_v2 = ones(size(K2,3),1);weight_v2 = ones(size(K2,3),1)/size(K2,3);
    K_COM2 = combine_kernels(weight_v2, K2);
	
	[lncRNA,disease]=ind2sub(size(y_train),index(i));
    if sum(y_train(lncRNA,:))==0
        y_train(lncRNA,:)=Sim_lnc(y_train,K_COM1,lncRNA);
    end
    if sum(y_train(:,disease))==0
        y_train(:,disease)=Sim_dis(y_train,K_COM2,disease);
    end
     
    [F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);
   
    Pre_value(II(i),JJ(i)) = F_1(II(i),JJ(i));
    y_train = interMatrix;
    if mod(i,10) == 0
        i      
    end
end

K1(:,:,3) = kernel_corr(y_train,1,0,1);
K2(:,:,3) = kernel_corr(y_train,2,0,1);



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

