function [ kd,kl ] = consine( interaction )
    [nd,nl] = size(interaction);
    

    %calculate the consine similarity between disease: kd
    kd = zeros(nd);
    for i = 1:nd
        for j = 1:nd
            denomi = norm(interaction(i, :)) * norm(interaction(j, :));
            if denomi~=0
                kd(i,j) =interaction(i,:)*(interaction(j,:))'/denomi;
            end
        end
    end

    %calculate the consine similarity between lncRNA: km
    kl = zeros(nl);
    for i = 1:nl
        for j = 1:nl
            denomi = norm(interaction(:,i))*norm(interaction(:,j));
            if denomi~=0
                kl(i,j) =(interaction(:,i))'*interaction(:,j)/denomi;
            end
        end
    end
end

