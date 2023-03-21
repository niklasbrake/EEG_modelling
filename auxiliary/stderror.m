function se = stderror(X)
	se = nanstd(X,0,1)./sqrt(sum(~isnan(X),1));
