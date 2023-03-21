function C = getMatrixCorrelation(spikeTrain,idcs1,idcs2)
    N = length(idcs1);
    M = length(idcs2);
    % C = zeros(N,M);
    C = zeros(N,3);
    idx = 1;
    if(N*M>1e5)
        h = waitbar(0,'Computing pairwise correlations...')
        for i = 1:N
            waitbar(i/N,h);
            for j = 1:M
                d = pdist2(spikeTrain{idcs1(i)},spikeTrain{idcs2(j)})<20;
                if(sum(d)>0)
                    C(idx,:) = [i,j,sum(d(:))/numel(d)];
                    idx=idx+1;
                    if(idx>size(C,1))
                        C = [C;zeros(N,3)];
                    end
                end
                % C(i,j) = sum(d(:)<100)/numel(d);
            end
        end
        C(isnan(C(:)))=0;
        delete(h);
    else
        for i = 1:N
            for j = 1:M
                d = pdist2(spikeTrain{idcs1(i)},spikeTrain{idcs2(j)});
                % C(i,j) = sum(d(:)<100)/numel(d);
            end
        end
        C(isnan(C(:)))=0;
    end