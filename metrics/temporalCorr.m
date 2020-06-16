function TCC = temporalCorr(X1,X2)

TCC = NaN(size(X1,1),1);
for i = 1:size(X1,1)
    TCC(i) = corr(X1(i,:)',X2(i,:)');
end
    
end