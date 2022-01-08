function c = getCoD(x,y)
%% Returns the coeffcient of determination between signals x and y

RSS = sum((x-y).^2);
TSS = sum((x-mean(x)).^2);

c = 1-(RSS/TSS);

end