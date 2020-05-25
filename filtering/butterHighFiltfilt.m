function X = butterHighFiltfilt(X, fs, fc)

[b,a] = butter(2, fc/(fs/2), 'high');
pad = 10*size(X,2);
X = [repmat(X(:,1),1,pad) X repmat(X(:,end),1,pad)];
X = filtfilt(b, a, X')';
X = X(:,1+pad:end-pad);

end