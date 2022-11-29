% This code numerically computes the first Eigenfunction of the
% infinity Laplacian operator on the unit square [-1 1]^2 for given
% boundary conditions.
% @authors: Farid Bozorgnia, Daniel Tenbrinck
% @date: 03/10/16

% tidy up workspace
clear; close all;

% start time measurement
tic

% initialize necessary parameters
reg_accur = 1e-07;
max_iterations = 1000;
samples = 50;

% initialize grid on unit square [-1 1]^2 according to amount of sample
% points
step_size = 2 / samples;
[x,y] = meshgrid(linspace(-1-step_size,1+step_size,samples+3));

% initialize ground truth
f = x.^3-3*x.*y.^2;

% set boundary values by defining function on whole domain and erasing
% inner part then
p1 = x.^3-3*x.*y.^2;
p1(3:end-2,3:end-2) = 0;

% visualize boundary conditions
figure(1); surf(x,y,p1);
pause(0.1);

% initialize matrix of values
values = zeros(12*25, prod(size(p1) - [4 4]));

% initialize matrix for new values of p1
hh = zeros(1,size(values,2));

% building the distance function for the following 5x5 stencil:
%     1   6   11  16  21
%     2   7   12  17  22
%     3   8    x  18  23
%     4   9   14  19  24
%     5   10  15  20  25
% NOTE: We use the full 5x5 neighborhood and set unnecessary pixels with a
% negative distance.
distance = step_size * ...
    [sqrt(8), sqrt(5), 2, sqrt(5), sqrt(8), ...
    sqrt(5), sqrt(2), 1, sqrt(2), sqrt(5), ...
    2, 1, -10, 1, 2, ...
    sqrt(5), sqrt(2), 1, sqrt(2), sqrt(5), ...
    sqrt(8), sqrt(5), 2, sqrt(5), sqrt(8)];

% initialize stopping criterions
rel_change = 1;
iteration = 1;

% initialize containers
rel_changes = zeros(1,max_iterations);
L2_error = rel_changes;

% iterate until convergence or maximum number of iterations is reached
while rel_change > reg_accur && iteration <= max_iterations
    
    % generate patch rows
    patches = im2col(p1, [5 5], 'sliding');
    
    % compute values abs( (v(k) - v(l) )/( d(k)+d(l) )
    for i = 1:12
        values( 1+(i-1)*25:i*25, :) = abs(patches - circshift(patches,[i,0])) ./ ...
            repmat(distance' + circshift(distance',[i,0]), [1 size(patches,2)] );
    end
    
    % compute maximum values in each column
    max_values = max(values, [], 1);
    
    % get indices of corresponding rows (not necessarily unique!)
    indices = zeros(1,size(patches,2));
    for i = 1:size(patches,2)
        candidates = find(values(:,i) == max_values(i));
        indices(i) = candidates(1); % choose the first candidate
    end
    
    % determine corresponding index pair within neighborhood
    index_pairs = zeros(2,size(indices,2));
    index_pair1 = mod(indices,25);
    index_pair1(index_pair1 == 0) = 25;
    index_pair2 = mod(index_pair1 - floor(indices / 25) - 1, 25);
    index_pair2(index_pair2 == 0) = 25;
    index_pairs(1,:) = index_pair1;
    index_pairs(2,:) = index_pair2;
    
    % compute new values of p1 as hh = (v(k1)*d(k2)+ v(k2)*d(k1))/(d(k1)+d(k2))
    % -> this corresponds to u* in the Oberman discretization
    for i = 1:size(patches,2)
        hh(i) = ( patches(index_pairs(1,i),i) * distance(index_pairs(2,i)) + ...
            patches(index_pairs(2,i),i) * distance(index_pairs(1,i)) ) / ...
            (distance(index_pairs(1,i)) + distance(index_pairs(2,i)));
    end
    hh = reshape(hh, size(p1) - [4 4]);
    
    % compute relative change between two iterations
    rel_change = norm(p1(3:end-2,3:end-2) - hh);
    rel_changes(iteration) = rel_change;
      
    % update p1
    p1(3:end-2,3:end-2) = hh;
    
    % compute L2 error of approximation to ground truth
    L2_error(iteration) = norm(p1(:) - f(:));
    
    % increase iteration counter
    iteration = iteration + 1;
    
end

% visualize result
figure(1); surf(x,y,p1);
pause(0.1);

% show relative error
figure; plot(L2_error);
figure; plot(rel_changes);

% stop time measurement
toc