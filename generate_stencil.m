function distance = generate_stencil(neighborhood_size, step_size,stencil_shape)
% NOTE: We use the full 5x5 neighborhood and set unnecessary pixels with a
% negative distance.

% define negative value for invalid elements
inval = nan;-1e+05;

if strcmp(stencil_shape, 'diamond')
  
  % building the distance function for the following 5x5 stencil:
  %         6       16
  %     2   7   12  17  22
  %         8    x  18
  %     4   9   14  19  24
  %         10      20
  distance = step_size * ...
    [  inval, sqrt(5), inval, sqrt(5),   inval, ...
     sqrt(5), sqrt(2),     1, sqrt(2), sqrt(5), ...
       inval,       1, inval,       1,   inval, ...
     sqrt(5), sqrt(2),     1, sqrt(2), sqrt(5), ...
       inval, sqrt(5), inval, sqrt(5),   inval];
  
elseif strcmp(stencil_shape, 'full')
  
  % building the distance function for the following 5x5 stencil:
  %     1   6   11  16  21
  %     2   7   12  17  22
  %     3   8    x  18  23
  %     4   9   14  19  24
  %     5   10  15  20  25
  % NOTE: We use the full 5x5 neighborhood and set unnecessary pixels with a
  % negative distance.
%   distance = step_size * ...
%     [sqrt(8), sqrt(5),     2, sqrt(5), sqrt(8), ...
%      sqrt(5), sqrt(2),     1, sqrt(2), sqrt(5), ...
%            2,       1, inval,       1,       2, ...
%      sqrt(5), sqrt(2),     1, sqrt(2), sqrt(5), ...
%      sqrt(8), sqrt(5),     2, sqrt(5), sqrt(8)];
  distance = zeros(neighborhood_size);
  distance(ceil(end/2), ceil(end/2)) = 1;
  distance = double(step_size * bwdist(distance));
  distance(ceil(end/2), ceil(end/2)) = inval;
  distance = reshape(distance,[1 numel(distance)]);
  
elseif strcmp(stencil_shape, 'square')
  
  % building the distance function for the following 5x5 stencil:
  %     1   6   11  16  21
  %     2               22
  %     3        x      23
  %     4               24
  %     5   10  15  20  25
  % NOTE: We use the full 5x5 neighborhood and set unnecessary pixels with a
  % negative distance.
%   distance = step_size * ...
%     [sqrt(8), sqrt(5),     2, sqrt(5), sqrt(8), ...
%      sqrt(5),   inval, inval,   inval, sqrt(5), ...
%            2,   inval, inval,   inval,       2, ...
%      sqrt(5),   inval, inval,   inval, sqrt(5), ...
%      sqrt(8), sqrt(5),     2, sqrt(5), sqrt(8)];

  distance = zeros(neighborhood_size);
  distance(ceil(end/2), ceil(end/2)) = 1;
  distance = double(step_size * bwdist(distance));
  distance(2:end-1, 2:end-1) = inval;
  distance = reshape(distance,[1 numel(distance)]);
  
end

end

