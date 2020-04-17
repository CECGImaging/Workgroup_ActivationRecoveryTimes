function X = gaussFiltfilt(X, sigma)

if sigma <= 0
    return;
end

l = ceil(3*sigma);
b = exp(-0.5*((-l:l)/sigma).^2);
b = b./sum(b);

X = [repmat(X(:,1),1,l) X repmat(X(:,end),1,l)];
X = fftfiltfilt(b, X')';
X = X(:,1+l:end-l);

end