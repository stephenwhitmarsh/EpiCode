function listing = dir_data(path)
%get only the data file (.ncs if neuralynx, .TRC if micromed, .eeg if
%Brainvision)

% if nargin == 0
%     name = '.';
% elseif nargin == 1
%     name = varargin{1};
% else
%     error('Too many input arguments.')
% end



listing = dir(path);

inds = [];
k    = 1;

while k <= length(listing)
    if strcmp(listing(k).name(1), '.') 
        inds(end + 1) = k;
    end
    k = k + 1;
end

listing(inds) = [];

inds = [];
k    = 1;

while k <= length(listing)
    if any(strcmp(listing(k).name(end-3:end), {'.mrk', 'vhdr', 'vmrk', '.mat', '.txt', '.bni'})) %extensions to ignore
        inds(end + 1) = k;
    end
    k = k + 1;
end

listing(inds) = [];