function [kk] = adapt_hist(r)
    [N,edges,bin]=histcounts(r);
    mb=max(bin);
    if(mb>1)
        kk=find(bin==mb);
    else
        kk=[];
    end
end