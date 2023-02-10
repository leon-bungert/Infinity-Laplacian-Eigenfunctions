% This code numerically computes second eigenfunctions, using the splitting
% in positive and negative parts.
% @authors: Farid Bozorgnia, Leon Bungert, Daniel Tenbrinck
% @date: 10/02/22

% tidy up workspace
clear;
close all;

% define problem
problem = 'eigenfunction';
% problem = 'infinity-harmonic';

% define square domain
domain = @(x,y) max(abs(x),abs(y)) <= 1;

% define boundary conditions
f = @(x,y) 0.*x.^0;

% either enforce positive or negative peaks or do normalization
% WARNING! Normalization is not guaranteed to produce correct solution
high_ridge = true;
normalization = false;

% define initialization
% init = 'zero';
init = 'random';
% init = 'distance';

% start time measurement
tic

% initialize necessary parameters
TOL             = 1e-7; 
max_iterations  = 2000;
samples         = 50;       % has to be even!!!!
stencil_shape   = 'full';   % valid stencils: full, square
neighborhood_size = 7;      % at least 3, higher stencil size yields instabilities

% initialize grid on unit square [-1 1]^2 according to amount of sample
% points
step_size = 2 / samples;
[x,y] = meshgrid(linspace(-1,1,samples+1));

% set boundary values by defining function on whole domain and erasing
% inner part then
p = zeros(size(x));
p([1 end],:)=-1;
p(:,[1 end])=-1;
p(find(~domain(x,y)))=-1;


distance_function = double(bwdist(p < 0)) * step_size;
[max_dist, ~] = max(distance_function(:));
max_idx = find(distance_function == max_dist);
p(max_idx)=-1;


second_distance_function = double(bwdist(p < 0)) * step_size;
[max_dist, ~] = max(second_distance_function(:));
lambda = 1 / max_dist;

dist_inner = second_distance_function(2:end-1, 2:end-1);
max_idx = find(dist_inner == max_dist);

% set peak values in bottom left and top right corner
max_idx = max_idx([1,end]);

% initialization
switch init
    case 'random'
        phi = randn(size(distance_function));        
    case 'distance'    
        phi = distance_function;
    case 'zero'
        phi = zeros(size(distance_function));
end
phi = max_dist*phi/sqrt(max(abs(phi(:)))+0.001^2);

u = phi(2:end-1,2:end-1);

% visualize initial conditions
figure(1); surf(u);
pause(0.01);

% building the distance function for a specified stencil
distance = generate_stencil(neighborhood_size, step_size, stencil_shape);
radius = floor(neighborhood_size/2);

% initialize boundary conditions on larger grid
bcinner = f(x,y);
bc = padarray(bcinner, [radius-1, radius-1], NaN);

xinn = x(2:end-1,2:end-1); yinn = y(2:end-1,2:end-1);
mask = domain(xinn, yinn);

% initialize matrix of values
values_pos = zeros(floor(neighborhood_size^2 /2)*neighborhood_size^2, numel(u));
values_neg = zeros(floor(neighborhood_size^2 /2)*neighborhood_size^2, numel(u));

% initialize matrix for values of ustar
ustar = zeros(1,size(values_pos,2));

% initialize stopping criterions
rel_change = inf;
scheme_accuracy = inf;
iteration = 0;

% initialize containers
rel_changes = zeros(1,max_iterations);

if high_ridge
%     u(max_idx) = max_dist;
    u(max_idx) = [max_dist, -max_dist];
end
if normalization
    u_pos = max(u,0);
    u_neg = max(-u,0);    
    u_pos = max_dist * u_pos/max(u_pos(:));
    u_neg = max_dist * u_neg/max(u_neg(:));    
    u = u_pos - u_neg;
end

while max(rel_change,0.) > TOL && iteration <= max_iterations  
        
    u_pos = max(u,0);
    u_neg = max(-u,0); 

    % generate patch rows
    pad_u_pos = bc;
    pad_u_neg = bc;

    pad_u_pos(radius+1:end-radius,radius+1:end-radius) = u_pos;
    pad_u_neg(radius+1:end-radius,radius+1:end-radius) = u_neg;
    patches_pos = im2col(pad_u_pos, [neighborhood_size neighborhood_size], 'sliding');
    patches_neg = im2col(pad_u_neg, [neighborhood_size neighborhood_size], 'sliding');

    % compute values abs( (v(k) - v(l) )/( d(k)+d(l) )
    for i = 1:floor(neighborhood_size^2 /2)
        values_pos( 1+(i-1)*neighborhood_size^2:i*neighborhood_size^2, :) = abs(patches_pos - circshift(patches_pos,[i,0])) ./ ...
            repmat(distance' + circshift(distance',[i,0]), [1 size(patches_pos,2)] );
        values_neg( 1+(i-1)*neighborhood_size^2:i*neighborhood_size^2, :) = abs(patches_neg - circshift(patches_neg,[i,0])) ./ ...
            repmat(distance' + circshift(distance',[i,0]), [1 size(patches_neg,2)] );
    end


    % compute maximum values in each column
    [~, indices_pos] = max(values_pos, [], 1);
    [~, indices_neg] = max(values_neg, [], 1);

    % determine corresponding index pair within neighborhood
    index_pairs_pos = zeros(2,size(indices_pos,2));
    index_pair1_pos = mod(indices_pos,neighborhood_size^2);
    index_pair1_pos(index_pair1_pos == 0) = neighborhood_size^2;
    index_pair2_pos = mod(index_pair1_pos - ceil(indices_pos / neighborhood_size^2) , neighborhood_size^2);
    index_pair2_pos(index_pair2_pos == 0) = neighborhood_size^2;
    index_pairs_pos(1,:) = index_pair1_pos;
    index_pairs_pos(2,:) = index_pair2_pos;

    index_pairs_neg = zeros(2,size(indices_neg,2));
    index_pair1_neg = mod(indices_neg,neighborhood_size^2);
    index_pair1_neg(index_pair1_neg == 0) = neighborhood_size^2;
    index_pair2_neg = mod(index_pair1_neg - ceil(indices_neg / neighborhood_size^2) , neighborhood_size^2);
    index_pair2_neg(index_pair2_neg == 0) = neighborhood_size^2;
    index_pairs_neg(1,:) = index_pair1_neg;
    index_pairs_neg(2,:) = index_pair2_neg;

    % compute new values of p1 as hh = (v(k1)*d(k2)+ v(k2)*d(k1))/(d(k1)+d(k2))
    % -> this corresponds to u* in the Oberman discretization
    for i = 1:size(patches_pos,2)
        ustar_pos(i) = ( patches_pos(index_pairs_pos(1,i),i) * distance(index_pairs_pos(2,i)) + ...
            patches_pos(index_pairs_pos(2,i),i) * distance(index_pairs_pos(1,i)) ) / ...
            (distance(index_pairs_pos(1,i)) + distance(index_pairs_pos(2,i)));
        ustar_neg(i) = ( patches_neg(index_pairs_neg(1,i),i) * distance(index_pairs_neg(2,i)) + ...
            patches_neg(index_pairs_neg(2,i),i) * distance(index_pairs_neg(1,i)) ) / ...
            (distance(index_pairs_neg(1,i)) + distance(index_pairs_neg(2,i)));
    end

    ustar_pos = reshape(ustar_pos, size(u_pos));
    ustar_neg = reshape(ustar_neg, size(u_neg));
    n_inf_L_pos = u_pos - ustar_pos;
    n_inf_L_neg = u_neg - ustar_neg;



    % positive upwind gradients 
    slopes_pos = ( u_pos(:)'- patches_pos ) ./ distance';
    [max_slopes, local_idx] = max(slopes_pos, [], 1);
    max_slopes = reshape(max_slopes, size(u_pos));
    distances_max = reshape(distance(local_idx), size(u_pos));
    [local_subs_i, local_subs_j] = ind2sub([neighborhood_size, neighborhood_size],local_idx);
    [center_subs_i, center_subs_j] = ind2sub(size(u_pos), 1:numel(u_pos));
    global_subs_i = center_subs_i + local_subs_i-1;
    global_subs_j = local_subs_j + center_subs_j - 1;
    global_idx = sub2ind(size(u_pos) + [2*radius 2*radius] , global_subs_i, global_subs_j);
    val_pos = reshape(pad_u_pos(global_idx), size(u_pos));

    grad_term_pos = ( u_pos - val_pos - lambda.*distances_max.*u_pos );

    % negative upwind gradients 
    slopes_neg = ( u_neg(:)'- patches_neg ) ./ distance';
    [min_slopes, local_idx] = max(slopes_neg, [], 1);
    min_slopes = reshape(min_slopes, size(u_neg));
    distances_min = reshape(distance(local_idx), size(u_neg));
    [local_subs_i, local_subs_j] = ind2sub([neighborhood_size, neighborhood_size],local_idx);
    [center_subs_i, center_subs_j] = ind2sub(size(u_neg), 1:numel(u_neg));
    global_subs_i = center_subs_i + local_subs_i-1;
    global_subs_j = local_subs_j + center_subs_j - 1;
    global_idx = sub2ind(size(u_neg) + [2*radius 2*radius] , global_subs_i, global_subs_j);
    val_neg = reshape(pad_u_neg(global_idx), size(u_neg));

    grad_term_neg = ( u_neg - val_neg - lambda.*distances_min.*u_neg );

    scheme_pos = min( grad_term_pos , n_inf_L_pos) ;%+ n_inf_L_neg ;
    scheme_neg = min( grad_term_neg, n_inf_L_neg ) ;%- n_inf_L_pos ;    

    rho = .9;
    u_pos_new = u_pos - rho * scheme_pos;
    u_neg_new = u_neg - rho * scheme_neg;
    
    u_pos_new = max(u_pos_new, 0);
    u_neg_new = max(u_neg_new, 0);

    
    if normalization 
        max_pos = max(u_pos_new(:));
        max_neg = max(u_neg_new(:));

        if max_pos > 0
            u_pos_new = max_dist * u_pos_new/max_pos;
        end
        if max_neg > 0
            u_neg_new = max_dist * u_neg_new/max_neg;
        end
    end
    unew = u_pos_new - u_neg_new;
    if high_ridge
        unew(max_idx) = [max_dist, -max_dist];
    end
    % compute relative change between two iterations
    diff = abs(u-unew).*mask;
    rel_change = norm(diff, 'inf') / norm(unew.*mask, 'inf');
%     scheme_accuracy = norm(scheme.*mask, 'inf');
        
    
    % update u
    u = unew.*mask;
    
       
    % increase iteration counter
    iteration = iteration + 1;
    
    pad_u = bcinner;
    pad_u(2:end-1, 2:end-1) = unew;
    solution = pad_u;
    
    % visualize current solution
    subplot(1,2,1);
    surf(x,y,solution); title('Solution'); axis equal; view(-30,20); drawnow;
    subplot(1,2,2);
    surf(xinn,yinn,mask.*(scheme_pos+scheme_neg)); title('Scheme'); view(-30,20); drawnow;
    
    % give some output
    disp(['Current iteration: ' num2str(iteration) ' with a '... 
            'relative change of ' num2str(rel_change)]);
       
end

