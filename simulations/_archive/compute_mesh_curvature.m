function [H,kg] = compute_mesh_curvature(V,F)

% Get surface area of each triangle
A0 = zeros(size(F,1),1);
for i = 1:size(F,1)
    tri = V(F(i,:),:);
    A0(i) = 1/2*norm(cross(tri(1,:)-tri(2,:),tri(3,:)-tri(1,:)));
end

pairs = nchoosek(1:3,2);
kg = zeros(size(V,1),1);
H = zeros(size(V,1),1);
for i = 1:size(V,1)
    fs = find(sum(F==i,2));
    kg(i) = 0;
    H(i) = 0;
    if(~isempty(fs))
        Ai = sum(A0(fs))/3;
        p_i = V(i,:);

        % Mean curvature
        temp = F(fs,:);
        U = unique(temp);
        if(any(sum(U(:)==temp(:)',2)<2))
            H(i) = 0;
        else
            U = setdiff(U,i)';
            deltaP = 0;
            for j = U
                p_j = V(j,:);
                v_adj = temp(find(sum(temp==j,2)),:);
                k1 = setdiff(v_adj(1,:),[i,j]);
                k2 = setdiff(v_adj(2,:),[i,j]);
                p_k1 = V(k1,:);
                p_k2 = V(k2,:);
                vs1 = p_i-p_k1;
                vs2 = p_j-p_k1;
                alpha_ij = acos(dot(vs1/norm(vs1),vs2/norm(vs2)));
                vs1 = p_i-p_k2;
                vs2 = p_j-p_k2;
                beta_ij = acos(dot(vs1/norm(vs1),vs2/norm(vs2)));
                deltaP = deltaP + (cot(alpha_ij)+cos(beta_ij))*(p_j-p_i);
            end
            deltaP = deltaP/(2*Ai);
            H(i) = norm(deltaP)/2;
        end

        % Gaussian curvature
        T = 0;
        for j = 1:length(fs)
            temp = F(fs(j),:);
            temp = temp(pairs);
            idcs = find(sum(temp==i,2));
            vs = V(temp(idcs,:),:)-p_i;
            vs = vs(sum(vs==0,2)<3,:);
            vs = vs./vecnorm(vs,2,2);
            T = T+acos(dot(vs(1,:),vs(2,:)));
        end
        kg(i) = (2*pi-T)/Ai;
    end
end