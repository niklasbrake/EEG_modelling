function R = bootstrap_column(Q,M,rho)

n = size(Q,3);
idcs = randi(100,M,1);
[C,~,ic] = unique(idcs);
a_counts = accumarray(ic,1);

R = zeros(size(Q,1),size(Q,2),M);
i = 0;
for j = 1:length(a_counts)
    R(:,:,i+1:i+a_counts(j)) = phaseran(Q(:,:,C(j)),a_counts(j));
    i = i+a_counts(j);
end

if(nargin==3)
    % X0 = phaseran(Q(:,:,randi(n)),1);
    X0 = mean(R,3);
    X0 = X0./std(X0);
    V = std(R);
    R = sqrt(rho)*X0 + sqrt(1-rho)*R./V;
    R = R.*V;
end