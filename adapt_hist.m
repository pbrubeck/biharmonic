function [kk] = adapt_hist(r)
    [N,edges,bin]=histcounts(r);
    mb=max(bin);
    if(mb>1)
        kk=find(bin==mb);
    else
        [ri,kk]=max(r);
        if(ri/sum(r)<2/numel(r)), kk=[]; end
    end
end