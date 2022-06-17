function c3d = filteringMarkersData(c3d)

markers = btkGetMarkers(c3d);
labels = fieldnames(markers);
freq = btkGetPointFrequency(c3d);

% define butterworth filter
fc = 15;
fs = freq;
[b,a] = butter(2,fc/(fs/2));
for i = 1:length(labels)
    v(:,1+3*(i-1):3+3*(i-1)) = filter(b,a,markers.(labels{i}),[],1);
    % correct the peak at the beginning
    v(1:12,1+3*(i-1):3+3*(i-1)) = markers.(labels{i})(1:12,:);
end

btkSetMarkersValues(c3d, v);
% markersFilt = btkGetMarkers(c3d);
end