function [w] = cka_kernels_weights(Kernels_list,adjmat,dim)

% adjmat : binary adjacency matrix
% dim    : dimension (1 - rows, 2 - cols)

num_kernels = size(Kernels_list,3);

weight_v = zeros(num_kernels,1);

y = adjmat;
    % Graph based kernel
if dim == 1
        ga = y*y';
else
        ga = y'*y;
end
%ga = Knormalized(ga);
N_U = size(y,dim);
l=ones(N_U,1);
U = eye(N_U) - (l*l')/N_U;


M = zeros(num_kernels,num_kernels);

for i=1:num_kernels
	for j=1:num_kernels
		kk1 = U*Kernels_list(:,:,i)*U;
		kk2 = U*Kernels_list(:,:,j)*U;
		mm = trace(kk1'*kk2);
	
		m1 = trace(kk1*kk1');
		m2 = trace(kk2*kk2');
		M(i,j) = mm/(sqrt(m1)*sqrt(m2));
	end
end

a = zeros(num_kernels,1);

for i=1:num_kernels

	kk = U*Kernels_list(:,:,i)*U;
	aa = trace(kk'*ga);
	a(i) = aa*(N_U-1)^-2;
end





v = randn(num_kernels,1);






falpha = @(v)obj_function(v,M,a);
        

[x_alpha, fval_alpha] = optimize_weights(v, falpha);

w = x_alpha;

end

function [J] = obj_function(v,M,a)

    J =v'*M*v - v'*a ;
end

function [x, fval] = optimize_weights(x0, fun)
    n = length(x0);
    Aineq   = [];
    bineq   = [];
    Aeq     = ones(1,n);
    beq     = 1;
    LB      = zeros(1,n);
    UB      = ones(1,n);

    options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'notify');
    [x,fval] = fmincon(fun,x0,Aineq,bineq,Aeq,beq,LB,UB,[],options);
end