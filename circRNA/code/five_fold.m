clear;
seed = 12345678;
rand('seed', seed);
nfolds = 5; nruns=1;
dataname = '';
dataname
% load adjacency matrix
[y,L1,L2] = loadtabfile(['../data/interactions/' dataname 'cd_adjmat.txt']);

fold_aupr=[];fold_auc=[];

k_nn = 10;IsHG =1;CVS=1;%%IsHG是不是超图，1代表是超图
k1 = 300;k2 = 900;

lambda1 = 2^0;lambda2=2^0;
globa_true_y_lp=[];
globa_predict_y_lp=[];
for run=1:nruns
    % split folds
	[num_D,num_G] = size(y);
	if CVS == 1
		crossval_idx = crossvalind('Kfold',y(:),nfolds);
	elseif CVS == 2
		KP = 1:1:num_D;
		crossval_idx = crossvalind('Kfold',KP,nfolds);
	elseif CVS == 3
		KP = 1:1:num_G;
		crossval_idx = crossvalind('Kfold',KP,nfolds);
	end
   

    for fold=1:nfolds
        train_idx = find(crossval_idx~=fold);
        test_idx  = find(crossval_idx==fold);

        y_train = y;
        
		if CVS == 1
			y_train(test_idx) = 0;
		elseif CVS == 2
			y_train(test_idx,:) = 0;
		elseif CVS == 3
			y_train(:,test_idx) = 0;
		end
		

        %%  1.kernels
		%% cal kernels
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
        K1(:,:,i+1) = kernel_corr(y_train,1,0,1);
		
		%K1 compelted
		
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
        K2(:,:,i+1) = kernel_corr(y_train,2,0,1);
		
		%K1 compelted
		
		
        %% perform predictions
        %lambda = 1;
		%K1 = K1(:,:,3);
		%K2 = K2(:,:,1);


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
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
		
		%y_train = preprocess_WKNKN(y_train,K_COM1,K_COM2,3,0.5);A_cos_com=y_train;
		%3.LGC
		[A_cos_com] = grtmf(K_COM1,K_COM2,y_train,lambda1,lambda2,k1,k2,k_nn,IsHG);
        %[A_cos_com] = cmf(K_COM1,K_COM2,y_train,88,500,lambda1,lambda2,0.1);
        %[A_cos_com] = grmf(K_COM1,K_COM2,y_train,88,k2,0.1,0.001,0.001,10);
        %[A_cos_com] = nrlmf(K_COM1,K_COM2,y_train,88,500,0.125,0.125,3,0.0001,0.0001,0.001,10);
        %[A_cos_com] = tmf(K_COM1,K_COM2,y_train,1,k1,k2);
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
        %% 4. evaluate predictions
        yy=y;
    
		if CVS == 1
			test_labels = yy(test_idx);
			predict_scores = A_cos_com(test_idx);
		elseif CVS == 2
			test_labels1 = yy(test_idx,:);
			test_labels = test_labels1(:);
			predict_scores = A_cos_com(test_idx,:);
			predict_scores = predict_scores(:);
		elseif CVS == 3
			test_labels1 = yy(:,test_idx);
			test_labels = test_labels1(:);
			predict_scores = A_cos_com(:,test_idx);
			predict_scores = predict_scores(:);
		end
		[X,Y,tpr,aupr_] = perfcurve(test_labels,predict_scores,1, 'xCrit', 'reca', 'yCrit', 'prec');
		
		[X,Y,THRE,AUC_,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(test_labels,predict_scores,1);

		
		fprintf('---------------\nRUN %d - FOLD %d  \n', run, fold)

		fprintf('%d - FOLD %d - weighted_kernels_: %f \n', run, fold, aupr_)
		

		fold_aupr=[fold_aupr;aupr_];
		fold_auc=[fold_auc;AUC_];

		
		globa_true_y_lp=[globa_true_y_lp;test_labels];
		globa_predict_y_lp=[globa_predict_y_lp;predict_scores];
		
		
    end
    
    
end
RMSE = sqrt(sum((globa_predict_y_lp-globa_true_y_lp).^2)/length(globa_predict_y_lp))



mean_aupr = mean(fold_aupr)
mean_auc = mean(fold_auc)

