%tju cs for bioinformatics 
clear;
dataname = '';
dataname
k_nn = 10;IsHG =0;
%k1 = 250;k2 = 400;
k1 = 300;k2 = 900;
lambda1 = 2^0;lambda2=2^0;
[y,L1,L2] = loadtabfile(['../data/interactions/' dataname 'cd_adjmat.txt']);
y_train = y;
[II,JJ] = find(y == 1);
k1_paths = {['../data/interactions/' dataname 'circRNA_Gene.txt'],...
					['../data/interactions/' dataname 'circRNA_miRNA.txt'],...
                    };
K1 = [];
for i=1:length(k1_paths)
	[y_,l1,l2] = loadtabfile(k1_paths{i});
	K1(:,:,i) = kernel_corr(y_,1,0,1);
end
k2_paths = {['../data/kernels/' dataname 'circSim.txt']
		   };
for j=1:length(k2_paths)
	i = i+1;
	[mat, labels] = loadtabfile(k2_paths{j});
	mat = process_kernel(mat);
	K1(:,:,i) = Knormalized(mat);
end
k3_paths = {['../data/interactions/' dataname 'disease_Gene.txt'],...
                   
                    ['../data/interactions/' dataname 'disease_miRNA.txt'],...
                    };
K2 = [];
for i=1:length(k3_paths)
	[y_,l1,l2] = loadtabfile(k3_paths{i});
	K2(:,:,i) = kernel_corr(y_,1,0,1);
end
k4_paths = {['../data/kernels/' dataname 'disease_sim.txt']
		   };
for k=1:length(k4_paths)
	i=i+1;
	[mat, labels] = loadtabfile(k4_paths{k});
	mat = process_kernel(mat);
	K2(:,:,i) = Knormalized(mat);
end
Pre_value =y;
weight=[];
for i = 1:length(II)    
    y_train(II(i),JJ(i)) = 0;
    K1(:,:,4)=kernel_corr(y_train,1,0,1);
	K2(:,:,4)=kernel_corr(y_train,2,0,1);

    [weight_v1] = cka_kernels_weights(K1,y_train,1);
	%weight_v1 = ones(size(K1,3),1);weight_v1 = ones(size(K1,3),1)/size(K1,3);
    K_COM1 = combine_kernels(weight_v1, K1);		
 
    [weight_v2] = cka_kernels_weights(K2,y_train,2);
	%weight_v2 = ones(size(K2,3),1);weight_v2 = ones(size(K2,3),1)/size(K2,3);
    K_COM2 = combine_kernels(weight_v2, K2);
     
    [F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);
   
    Pre_value(II(i),JJ(i)) = F_1(II(i),JJ(i));
    y_train = y;
    if mod(i,10) == 0
        i      
    end
end


K1(:,:,4)=kernel_corr(y_train,1,0,1);
K2(:,:,4)=kernel_corr(y_train,2,0,1);


[weight_v1] = cka_kernels_weights(K1,y_train,1);
K_COM1 = combine_kernels(weight_v1, K1);		
 
[weight_v2] = cka_kernels_weights(K2,y_train,2);
K_COM2 = combine_kernels(weight_v2, K2);
     
[F_1] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);

for i =1:length(F_1(:))
    if y(i)==0
        Pre_value(i)=F_1(i);
    end
end

[X_1,Y_1,tpr,aupr_1] = perfcurve(y(:), Pre_value(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(y(:), Pre_value(:),1);

