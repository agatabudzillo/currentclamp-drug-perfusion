function [ spike_intervals ] = getspikeintervals( locs )

% getspikeintervals is calculates spike intervals given an array of spike times
% in a way that ensures that there are no floating point errors (which I had 
% previously encountered when dealing with long traces with Matlab student version)
        
locs_cropping = locs;
spike_intervals = zeros(length(locs),1);
q = 1;
 
while q <= length(locs)-1
locs_cropping = locs_cropping - locs_cropping(1);
locs_cropping = locs_cropping(2:end);
spike_intervals(q) = locs_cropping(1);
q = q+1;
end

end

