function [kd,kl] = GaussianKernel(interaction, gamadd, gamall)
    [nd,nl] = size(interaction);
    %calculate gamad for Gaussian kernel calculation
    sd = zeros(nd, 1);
    sl = zeros(nl, 1);
    %disease lncRNA
    for i = 1:nd
        sd(i) = norm(interaction(i,:))^2;
    end
    gamad = nd/sum(sd)*gamadd;

    %calculate gamal for Gaussian kernel calculation
    for i = 1:nl
        sl(i) = norm(interaction(:,i))^2;
    end
    gamal = nl/sum(sl)*gamall;

    %calculate Gaussian kernel for the similarity between disease: kd
    kd = zeros(nd);
    for i = 1:nd
        for j = 1:nd
            kd(i,j) = exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
        end
    end

    %calculate Gaussian kernel for the similarity between microbe: km
    kl = zeros(nl);
    for i = 1:nl
        for j = 1:nl
            kl(i,j) = exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
        end
    end
end