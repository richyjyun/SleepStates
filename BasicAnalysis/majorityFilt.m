function filtidx = majorityFilt(idx, window)
%% Returns the majority filtered version of data, using width window

filtidx = zeros(length(idx),1);
for i = 1:length(idx)
    left = max(1,i-floor(window/2)); right = min(length(idx),i+floor(window/2));
    c = histcounts(idx(left:right),1:5);
    [~,ind] = max(c);
    if length(ind)==1
        filtidx(i) = ind;
    else
        if any(ind==idx(i))
            filtidx(i) = idx(i);
        else
            filtidx(i) = ind(1);
        end
    end
end

end