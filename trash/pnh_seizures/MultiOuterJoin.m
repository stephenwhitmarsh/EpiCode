function Joined = MultiOuterJoin(varargin)
Joined = varargin{1};
for k = 2:nargin
    Joined = outerjoin(Joined, varargin{k},'MergeKeys', true);
end
end