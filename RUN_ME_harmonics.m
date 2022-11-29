% This code numerically computes the first Eigenfunction of the
% infinity Laplacian operator on the unit square [-1 1]^2 for given
% boundary conditions.
% @authors: Farid Bozorgnia, Daniel Tenbrinck
% @date: 08/10/16

% tidy up workspace
clear; close all;

% start time measurement
tic

% initialize necessary parameters
reg_accur = 1e-07;
max_iterations = 1000;
samples = 50;
neighborhood_size = 11;
stencil_shape = 'full'; % valid names: diamond, full, square

% compute radius of neighborhood
nb_half = floor(neighborhood_size/2);

% initialize grid on unit square [-1 1]^2 according to amount of sample
% points
step_size = 2 / samples;
[x,y] = meshgrid(linspace(-1-step_size*nb_half,1+step_size*nb_half,samples+neighborhood_size));
% x = (x+1)/2;
% y = (y+1)/2;

% initialize ground truth
% f = x.^3-3*x.*y.^2;
f = abs(x).^(4/3) - abs(y).^(4/3);

% Z1 = sqrt(x.^2 + (y-0.5).^2);
% Z2 = sqrt((x-1).^2 + (y-0.5).^2);
% f = 0.5*(Z1-Z2+1);

% set boundary values by defining function on whole domain and erasing
% inner part then
p1 = f;
% p1(nb_half+1:end-nb_half,nb_half+1:end-nb_half) = 0;

% visualize boundary conditions
figure(1); surf(x,y,p1);
pause(0.01);

% initialize matrix of values
values = zeros(neighborhood_size^2 * (floor(neighborhood_size^2/2) - 1), prod(size(p1) - 2*[nb_half nb_half]));

% initialize matrix for new values of p1
hh = zeros(1,size(values,2));

% building the distance function for a specified 5x5 stencil
distance = generate_stencil(neighborhood_size, step_size, stencil_shape);

% initialize stopping criterions
rel_change = 1;
iteration = 1;

% initialize containers
rel_changes = zeros(1,max_iterations);
L2_error = rel_changes;

% iterate until convergence or maximum number of iterations is reached
while rel_change > reg_accur && iteration <= max_iterations
    
    % generate patches for each pixel as columns
    patches = im2col(p1, [neighborhood_size neighborhood_size], 'sliding');
    
    % compute values abs( (v(k) - v(l) )/( d(k)+d(l) )
    for i = 1:floor(neighborhood_size^2/2) - 1
        values( 1+(i-1)*neighborhood_size^2:i*neighborhood_size^2, :) = abs(patches - circshift(patches,[i,0])) ./ ...
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
    index_pair1 = mod(indices,neighborhood_size^2);
    index_pair1(index_pair1 == 0) = neighborhood_size^2;
    index_pair2 = mod(index_pair1 - floor(indices / neighborhood_size^2) - 1, neighborhood_size^2);
    index_pair2(index_pair2 == 0) = neighborhood_size^2;
    index_pairs(1,:) = index_pair1;
    index_pairs(2,:) = index_pair2;
    
    % compute new values of p1 as hh = (v(k1)*d(k2)+ v(k2)*d(k1))/(d(k1)+d(k2))
    % -> this corresponds to u* in the Oberman discretization
    for i = 1:size(patches,2)
        hh(i) = ( patches(index_pairs(1,i),i) * distance(index_pairs(2,i)) + ...
            patches(index_pairs(2,i),i) * distance(index_pairs(1,i)) ) / ...
            (distance(index_pairs(1,i)) + distance(index_pairs(2,i)));
    end
    hh = reshape(hh, size(p1) - 2*[nb_half nb_half]);
      
    % compute relative change between two iterations
    rel_change = norm(p1(nb_half+1:end-nb_half,nb_half+1:end-nb_half) - hh);
    rel_changes(iteration) = rel_change;
      
    % update p1
    p1(nb_half+1:end-nb_half,nb_half+1:end-nb_half) = hh;
    surf(p1); drawnow;
    % compute L2 error of approximation to ground truth
    L2_error(iteration) = norm(p1(:) - f(:));
    
    % increase iteration counter
    iteration = iteration + 1;
    
end

% visualize result
figure(1); surf(x,y,p1);
pause(0.1);

% restrict plots to performed iterations
L2_error = L2_error(1:iteration-1);
rel_changes = rel_changes(1:iteration-1);

% show relative errors
figure; plot(L2_error);
figure; plot(rel_changes);

max(abs(f(:) - p1(:)))
L2_error(iteration-1)

% stop time measurement
toc