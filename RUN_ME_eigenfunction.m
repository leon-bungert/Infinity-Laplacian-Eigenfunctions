% This code numerically computes the first Eigenfunction of the
% infinity Laplacian operator on the unit square [-1 1]^2 for given
% boundary conditions.
% @authors: Farid Bozorgnia, Leon Bungert, Daniel Tenbrinck
% @date: 10/01/22

% tidy up workspace
clear all;
close all;

% define problem
problem = 'eigenfunction';                       
% problem = 'infinity-harmonic';

% define domain
shape = 'square';
% shape = 'square-disk';
% shape = 'half-disk';
% shape = 'disk';
% shape = 'ellipsis';
% shape = 'L-shape';
% shape = 'trapez';
% shape = 'rectangle';
% shape = 'stadium';
% shape = 'triangle';
% shape = 'triangle2';
% shape = 'dumbell';
% shape = 'two-disks';
% shape = 'heart1';
% shape = 'heart2';

% define boundary conditions, mainly for harmonic equation
% f = @(x,y) (abs(x).^(4/3) - abs(y).^(4/3));
% f = @(x,y) x.^3-3.*x.*y.^2;
f = @(x,y) 0.*x.^0;
% f = @(x,y) (x==1).*(abs(y)<0.01);

% either enforce positive or negative peaks or do normalization
% WARNING! Normalization is not guaranteed to produce correct solution,
high_ridge      = true;   
normalization   = false;

% define initialization
init = 'zero';
% init = 'random';
% init = 'distance';

% initialize necessary parameters
alpha           = 1;      % concave approximation parameter <= 1
TOL             = 1*1e-7; 
max_iterations  = 5000;
max_broyden     = 0;
broyden         = 'bad';
samples         = 1*50;       % has to be even!!!!
stencil_shape   = 'full';   % valid stencils: full, square
neighborhood_size = 9;     % at least 3
save2disk       = true;     % save results to disk
visualize       = true;     % visualize the results
savingFreq      = 10;       % saving frequency
visFreq         = 100;       % visualization frequency

% catch normalization in case of inf-harmonic problem
if strcmpi(problem, 'infinity-harmonic')
    high_ridge = true;
    normalization = false;
end

% initialize grid on unit square [-1 1]^2 according to amount of sample
% points
step_size = 2 / samples;


[x,y] = meshgrid(linspace(-1,1,samples+1));

% define subdomain of the square domain
switch shape
    case 'square'
        domain = @(x,y) max(abs(x),abs(y)) <= 1;          
    case 'square-disk'
        domain = @(x,y) (x.^2 + y.^2 < 1)|(x > 0);
        closed_domain = @(x,y) (x.^2 + y.^2 <= 1)|(x >= 0);
    case 'half-disk'
        domain = @(x,y) (x.^2 + y.^2 < 1).*(x<0);
        closed_domain = @(x,y) (x.^2 + y.^2 <= 1).*(x<=0); 
    case 'disk'
        domain = @(x,y) x.^2 + y.^2 < 1;
        closed_domain = @(x,y) x.^2 + y.^2 <= 1;        
    case 'ellipsis'
        domain = @(x,y) x.^2 + (4*y).^2 < 1;
        closed_domain = @(x,y) x.^2 + (4*y).^2 <= 1;
    case 'L-shape'
        domain = @(x,y) (x<0)|(y>0);
        closed_domain = @(x,y) (x<=0)|(y>=0);   
    case 'trapez'
        domain = @(x,y) ((x<0)|(y>0)).*(y<-x);
        closed_domain = @(x,y) ((x<=0)|(y>=0)).*(y<=-x);   
    case 'rectangle'
        domain = @(x,y) abs(y)<0.5;
        closed_domain = @(x,y) abs(y)<=0.5;
    case 'stadium'
        domain = @(x,y) ((abs(y)<0.5).*(abs(x)<0.5))|((x-0.5).^2+y.^2<0.5^2)|((x+0.5).^2+y.^2<0.5^2);
        closed_domain = @(x,y) ((abs(y)<=0.5).*(abs(x)<=0.5))|((x-0.5).^2+y.^2<=0.5^2)|((x+0.5).^2+y.^2<=0.5^2);
    case 'triangle'
        domain = @(x,y) (x<y);
        closed_domain = @(x,y) (x<=y);
    case 'triangle2'
        domain = @(x,y) (abs(x)<y);
        closed_domain = @(x,y) (abs(x)<=y);        
    case 'dumbell'
        domain = @(x,y) ((x+1).^2+y.^2<1)|((x-1).^2+y.^2<1)|(abs(y)<0.5);  
        closed_domain = @(x,y) ((x+1).^2+y.^2<=1)|((x-1).^2+y.^2<=1)|(abs(y)<=0.5);  
    case 'two-disks'
        domain = @(x,y) ((x+.5).^2+y.^2<.5^2)|((x-.5).^2+y.^2<.5^2)|((abs(y)<0.2).*(abs(x)<0.5)); 
        closed_domain = @(x,y) ((x+.5).^2+y.^2<=.5^2)|((x-.5).^2+y.^2<=.5^2)|((abs(y)<=0.2).*(abs(x)<=0.5)); 
    case 'heart1'
        domain = @(x,y) (((y+0.5).^2+(x).^2<0.5^2).*(x<=0.5))|(((y-0.5).^2+(x).^2<0.5^2).*(x<=0.5))|((x < 1-0.75*abs(y)).*(x>0).*(-0.75<=y).*(y<0.75));
        closed_domain = @(x,y) (((y+0.5).^2+(x).^2<=0.5^2).*(x<=0.5))|(((y-0.5).^2+(x).^2<=0.5^2).*(x<=0.5))|((x <= 1-0.75*abs(y)).*(x>=0).*(-0.75<=y).*(y<=0.75));
    case 'heart2'
        domain = @(x,y) (((y+0.5).^2+(x+0.5).^2<0.48^2).*(x<=-0.25))|(((y-0.5).^2+(x+0.5).^2<0.48^2).*(x<=-0.25))|((x < 1-1.35*abs(y)).*(x>-0.25))|((x<0).*(-0.5<x).*(abs(y)<0.5));
        closed_domain = @(x,y) (((y+0.5).^2+(x+0.5).^2<=0.48^2).*(x<=-0.25))|(((y-0.5).^2+(x+0.5).^2<=0.48^2).*(x<=-0.25))|((x <= 1-1.35*abs(y)).*(x>-0.25))|((x<=0).*(-0.5<=x).*(abs(y)<0.5));        
end     


if ~exist('closed_domain','var')
    closed_domain = domain;
end

% set boundary values by defining function on whole domain and erasing
% inner part then
p = zeros(size(x));
p([1 end],:)=-1;
p(:,[1 end])=-1;
p(find(~domain(x,y)))=-1;


distance_function = double(bwdist(p < 0)) * step_size;
dist_inner = distance_function(2:end-1, 2:end-1);
[max_dist, ~] = max(dist_inner(:));
lambda = 1/max_dist;

max_idx = find(dist_inner == max_dist);

max_dist_tmp = max_dist;
max_dist = 1;

center_idx = (length(max_idx)-1)/2+1;
max_idx = max_idx(center_idx);

% max_idx = max_idx(1);
% max_idx = [max_idx(1),  max_idx(end)];
%  max_idx = 1105; 
% max_idx = 4513;

% initialization
switch init
    case 'random'
        phi = rand(size(distance_function));        
    case 'distance'    
        phi = distance_function / max_dist_tmp;
    case 'zero'
        phi = zeros(size(distance_function));
end

u = phi(2:end-1,2:end-1);
u = u(:);

% building the distance function for a specified stencil
distance = generate_stencil(neighborhood_size, step_size, stencil_shape);
radius = floor(neighborhood_size/2);

% initialize boundary conditions on larger grid
bcinner = f(x,y);
bc = padarray(bcinner, [radius-1, radius-1], NaN);

% compute masks for the domain
xinn = x(2:end-1,2:end-1); yinn = y(2:end-1,2:end-1);
mask = domain(xinn, yinn);
mask = mask(:);
mask_vis = double(closed_domain(x, y));
mask_vis(mask_vis==0) = nan;

% initialize matrix of values
values = zeros(floor(neighborhood_size^2 /2)*neighborhood_size^2, numel(u));

% initialize matrix for values of ustar
ustar = zeros(size(values,2),1);

% initialize stopping criterions
rel_change = inf;
scheme_accuracy = inf;
iteration = 0;

% create output folder
if save2disk
    foldername = [problem, '_', shape, '_', 'high_ridge_', num2str(high_ridge),...
                 '_norm_', num2str(normalization), '_init_', init,... 
                 '_nbrs_', num2str(neighborhood_size)];
    if alpha < 1
        foldername = [foldername, '_alpha_', num2str(alpha)];
    end
             
    outputfolder = ['results','/',foldername];         

    if ~exist('results','dir')
        mkdir('results')
    end

    if ~exist(outputfolder, 'dir')
        mkdir(outputfolder)
    end
end

% initialize
if high_ridge
    u(max_idx) = max_dist;
end
if normalization
    u = max(u,0);        
    u = max_dist * u / norm(u,'inf');
end
B = speye(length(u));
umm = u;
um = u;
fm = zeros(size(u));
fmm = zeros(size(u));
%% iterate until convergence or maximum number of iterations is reached
fig = figure;
fig.PaperOrientation = 'landscape';
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',20);
% start time measurement
tic
while max(rel_change,scheme_accuracy) > TOL && iteration <= max_iterations    
    % generate patch rows
    pad_u = bc;
    pad_u(radius+1:end-radius,radius+1:end-radius) = reshape(u, size(xinn));
    patches = im2col(pad_u, [neighborhood_size neighborhood_size], 'sliding');
    
    % compute values abs( (v(k) - v(l) )/( d(k)+d(l) )
    for i = 1:floor(neighborhood_size^2 /2)
        values( 1+(i-1)*neighborhood_size^2:i*neighborhood_size^2, :) = abs(patches - circshift(patches,[i,0])) ./ ...
            repmat(distance' + circshift(distance',[i,0]), [1 size(patches,2)] );
    end
    
    
    % compute maximum values in each column
    [~, indices] = max(values, [], 1);

    % determine corresponding index pair within neighborhood
    index_pairs = zeros(2,size(indices,2));
    index_pair1 = mod(indices,neighborhood_size^2);
    index_pair1(index_pair1 == 0) = neighborhood_size^2;
    
    index_pair2 = mod(index_pair1 - ceil(indices / neighborhood_size^2) , neighborhood_size^2);
    index_pair2(index_pair2 == 0) = neighborhood_size^2;
    
    index_pairs(1,:) = index_pair1;
    index_pairs(2,:) = index_pair2;
            
    % compute new values of p1 as hh = (v(k1)*d(k2)+ v(k2)*d(k1))/(d(k1)+d(k2))
    % -> this corresponds to u* in the Oberman discretization
    for i = 1:size(patches,2)
        ustar(i) = (patches(index_pairs(1,i),i) * distance(index_pairs(2,i)) + ...
            patches(index_pairs(2,i),i) * distance(index_pairs(1,i)) ) / ...
            (distance(index_pairs(1,i)) + distance(index_pairs(2,i)));
    end
    
    F2 = u - ustar;
    if strcmpi(problem, 'infinity-harmonic')
        scheme = F2;
    elseif strcmpi(problem, 'eigenfunction')       
        % upwind gradients with general stencil
        slopes = ( u - patches' ) ./ distance;
        [max_slopes, local_idx] = max(slopes, [], 2);
        distances = distance(local_idx)';
        [local_subs_i, local_subs_j] = ind2sub([neighborhood_size, neighborhood_size],local_idx);
        [center_subs_i, center_subs_j] = ind2sub(size(xinn), 1:numel(u));
        global_subs_i = center_subs_i' + local_subs_i - 1;
        global_subs_j = center_subs_j' + local_subs_j - 1;
        global_idx = sub2ind(size(xinn) + [2*radius 2*radius] , global_subs_i, global_subs_j);
        
        val = pad_u(global_idx);
                        
        F1 = ((u - val ) - lambda.*distances.*sign(u).*abs(u).^alpha) ;
        

        scheme = min( F1 , F2);
    else
        error('Unknown problem!')
    end

    rho = 0.9;
    if and(iteration >= 2, iteration <= max_broyden)
        du = um-umm;
        df = fm-fmm;
        if strcmpi(broyden, 'good')
            B = B + ((df - B * du) * du')/norm(du)^2;
            unew = u - B \ scheme;
        elseif strcmpi(broyden, 'bad')
            B = B + ((du - B * df) * df')/norm(df)^2;
            unew = (u -  B * scheme); 
        else
            error('Unkown Broyden method')
        end
    else
        unew = (u - rho .* scheme); 
    end
    
    fmm = fm;
    fm = scheme;
    umm = um;
    um = u;
    
    if high_ridge
       unew(max_idx) = max_dist;
    end
    if normalization
        unew = max(unew,0);        
        unew = max_dist * unew / norm(unew,'inf');
    end
    % compute difference between two iterations
    diff = abs(u-unew).*mask;
  
         
    % update u and pad boundary values
    if strcmpi(problem, 'infinity-harmonic')
        scheme(max_idx) = 0;
    end
    u = unew.*mask;
    pad_u = bcinner;
    pad_u(2:end-1, 2:end-1) = reshape(um, size(xinn));
    solution = pad_u;
    scheme_exp = padarray(reshape(scheme.*mask, size(xinn)),[1,1]);
    
    % compute stopping criteria
    rel_change = norm(diff, 'inf') / norm(unew.*mask, 'inf');
    scheme_accuracy = norm(scheme.*mask, 'inf');
   
    % save to disk
    if save2disk && mod(iteration, savingFreq) == 0
        filename = ['solution_', 'grid_size_', strrep(num2str(step_size),'.','-'),... 
            '_iteration_', num2str(iteration)];
        save([outputfolder,'/',filename,'.mat'],'solution','scheme_exp')     
    end
     
    % visualize current solution
    if visualize && mod(iteration, visFreq) == 0
        fig;
        subplot(1,2,1);
        surf(x,y,solution.*mask_vis); title(['Solution at iteration ',num2str(iteration)]); 
        axis equal;
        view(-30,20); drawnow;
        subplot(1,2,2);
        surf(x,y,scheme_exp.*mask_vis); title(['Scheme at iteration ',num2str(iteration)]); 
        view(-30,20); 
        drawnow;
    end
    
    % give some output
    disp(['Current iteration: ' num2str(iteration) ... 
            ', relative change ' num2str(rel_change) ', accuracy ' num2str(scheme_accuracy)]);
    
    % increase iteration counter
    iteration = iteration + 1;   
end
toc
%%
if save2disk
        filename = ['solution_', 'grid_size_', strrep(num2str(step_size),'.','-'),... 
            '_iteration_', num2str(iteration-1)];      
        save([outputfolder,'/',filename,'.mat'],'solution','scheme_exp')       
end

if visualize 
        solution = max_dist_tmp * solution;
        scheme_exp = max_dist_tmp * scheme_exp;

        subplot(1,2,1);
        surf(x,y,solution.*mask_vis); title(['Solution at iteration ',num2str(iteration)]); 
        axis equal;
        view(-30,20); drawnow;
        subplot(1,2,2);
        surf(x,y,scheme_exp.*mask_vis); title(['Scheme at iteration ',num2str(iteration)]);
        view(-30,20); drawnow;
        if save2disk
            print([outputfolder,'/',filename],'-dpdf','-r300','-fillpage')
        end
        
        fig2 = figure(2);
        fig2.PaperOrientation = 'landscape';
        surf(x,y,solution.*mask_vis);  
        axis equal;
%         zlim([0, max_dist]);
%         axis([-1, 1, -1, 1, 0, max_dist]);
        view(-30,20); drawnow;
        if save2disk
            print([outputfolder,'/','final_surf_',filename],'-dpdf','-r300','-fillpage')
        end
        
        fig3 = figure(3);
        fig3.PaperOrientation = 'landscape';
        contour(x.*closed_domain(x,y),y.*closed_domain(x,y),solution);  
        axis equal;
        drawnow;
        if save2disk
            print([outputfolder,'/','final_contour_',filename],'-dpdf','-r300','-fillpage')
        end
end









