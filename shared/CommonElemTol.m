%% from: https://fr.mathworks.com/matlabcentral/answers/444501-how-to-use-intersect-command-with-a-tolerance-value

function [AI, BI] = CommonElemTol(A, B, Tol)
   A  = A(:);
   B  = B(:);
   nA = numel(A);
   M  = zeros(1, nA);
   
   % Collect the index of the first occurrence in B for every A:
   for iA = 1:nA
      dist = abs(A(iA) - B);             % EDITED: Of course abs() is needed
      Ind  = find(dist < Tol, 1);        % Absolute tolerance
      % Ind = find(dist ./ A(iA) < Tol, 1);  % Relative tolerance
      
      if ~isempty(Ind)
         M(iA) = Ind;
      end
   end
   AI = find(M);        % If any occurrence was found, this A exists
   if isempty(AI)       % prevent: Empty matrix: 0-by-1
      AI = [];
   end
   BI = M(AI);          % at this index in B
end