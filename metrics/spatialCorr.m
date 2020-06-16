function SCC = spatialCorr(X1,X2)

SCC = NaN(size(X1,2),1);
for i = 1:size(X1,2)
    SCC(i) = corr(X1(:,i),X2(:,i));
end
    
end