addpath 'data'
rng(0)

% The name of the test scene
test                = 'adam';
% The used inlier-outlier threshold in pixels
threshold           = 3.0;
% The used inlier-outlier threshold in pixels
confidence          = 0.9999;
% The maximum number of iterations in RANSAC
iteration_limit     = 5e3;
% The adaptively updated maximum iteration number in RANSAC
max_iterations      = iteration_limit;
% The truncated threshold for MSAC scoring
truncated_threshold = threshold * 3 / 2;

%% Load data
% Loading the source image
img1 = imread(strcat(test, 'A.png'));
% Loading the destination image
img2 = imread(strcat(test, 'B.png'));
% Loading the point correspondences
M = load(strcat(test, '.pts'));
% The number of points
N = size(M, 1);

% The points in the first image
source_points = M(:,1:3)';
% The points in the second image
destination_points = M(:,4:6)';

%% RANSAC homography fitting 
best_inliers = []; % The final set of inliers
best_homography = []; % The final homography
best_score = 0; % The score of the final homography

disp('RANSAC with 2SIFT solver started...')
for iter = 1 : iteration_limit
    % The randomly selected minimal sample
    indices = randperm(N, 2);

    % The points of the minimal sample in the images
    sample_src_points = M(indices, 1:3);
    sample_dst_points = M(indices, 4:6);

    % Calculating the normalizing transformation and normalizing
    % the point coordinates for estimating the homography.
    [normalized_src_points, normalized_dst_points, T1, T2] = ...
        normalizePoints(sample_src_points, sample_dst_points);
    
    % Scale difference between the two normalizing transformations.
    % This is requires to normalize the angles and scales
    normalizing_scale = T2(1,1) / T1(1,1);
    
    if isnan(normalized_src_points(1)) || isnan(normalized_dst_points(1))
        continue;
    end

    % The coordinates of the first point in the first image      
    u11 = normalized_src_points(1, 1);
    v11 = normalized_src_points(1, 2);
    % The coordinates of the first point in the second image
    u21 = normalized_dst_points(1, 1);
    v21 = normalized_dst_points(1, 2);
    % The coordinates of the second point in the first image    
    u12 = normalized_src_points(2, 1);
    v12 = normalized_src_points(2, 2);
    % The coordinates of the second point in the second image
    u22 = normalized_dst_points(2, 1);
    v22 = normalized_dst_points(2, 2);

    % The SIFT scale of the first point in the first image     
    q11 = M(indices(1), 7);
    % The SIFT scale of the first point in the second image
    q21 = M(indices(1), 8);
    % The SIFT rotation angle of the first point in the first image
    a11 = M(indices(1), 9);
    % The SIFT rotation angle of the first point in the second image
    a21 = M(indices(1), 10);
    % The SIFT scale of the second point in the first image
    q12 = M(indices(2), 7);
    % The SIFT scale of the second point in the second image
    q22 = M(indices(2), 8);
    % The SIFT rotation angle of the second point in the first image
    a12 = M(indices(2), 9);
    % The SIFT rotation angle of the second point in the second image
    a22 = M(indices(2), 10);

    % The sines and cosines of the rotation angles
    s11 = sin(a11);
    c11 = cos(a11);
    s21 = sin(a21);
    c21 = cos(a21);

    s12 = sin(a12);
    c12 = cos(a12);
    s22 = sin(a22);
    c22 = cos(a22);

    %% Application of the proposed solver
    % Compute the null space of the coefficient matrix formed from the
    % linear equations
    A = zeros(6, 9);
    A(1, :)     = [0, 0, 0, u11, v11, 1, -u11 * v21, -v11 * v21, -v21];
    A(2, :)     = [u11, v11, 1, 0, 0, 0, -u11 * u21, -v11 * u21, -u21];
    A(3, :)     = [0, 0, 0, u12, v12, 1, -u12 * v22, -v12 * v22, -v22];
    A(4, :)     = [u12, v12, 1, 0, 0, 0, -u12 * u22, -v12 * u22, -u22];
    A(5, :)     = [-s21 * c11, -s11 * s21, 0, c11 * c21, s11 * c21, 0, u21 * s21 * c11 - v21 * c11 * c21, u21 * s11 * s21 - v21 * s11 * c21, 0];
    A(6, :)     = [-s22 * c12, -s12 * s22, 0, c12 * c22, s12 * c22, 0, u22 * s22 * c12 - v22 * c12 * c22, u22 * s12 * s22 - v22 * s12 * c22, 0];

    n = null(A);

    data = zeros(47, 1);

    data(1:9) = n(:,1);
    data(10:18) = n(:,2);
    data(19:27) = n(:,3);

    % Coordinates of the first point correspondence
    data(28) = u11;
    data(29) = v11;
    data(30) = u21;
    data(31) = v21;

    % Coordinates of the second point correspondence 
    data(32) = u12;
    data(33) = v12;
    data(34) = u22;
    data(35) = v22;

    % SIFT parameters of the first correspondence
    data(36) = q11;
    data(37) = q21;
    data(38) = s11;
    data(39) = c11;
    data(40) = s21;
    data(41) = c21;

    % SIFT parameters of the second correspondence  
    data(42) = q12;
    data(43) = q22;
    data(44) = s12;
    data(45) = c12;
    data(46) = s22;
    data(47) = c22;

    % The normalizing scale
    k1 = normalizing_scale;
    k2 = normalizing_scale;
    data(48) = 1 / k1;
    data(49) = 1 / k2;

    % Apply the proposed solver operating on the null space
    Hs = solverHomography2SIFT(data);

    % Iterating through the possible solutions and selecting the one
    % with the maximum score
    min_error = 1e10;
    sol_number = 0;
    best_homography = [];
    for i = 1 : size(Hs, 2)
        alpha = Hs(1,i);
        beta = Hs(2,i);
        if abs(imag(alpha)) > 0 || abs(imag(beta)) > 0 % Keep just the real solutions
            continue;
        end

        % Recovering the homography from alpha and beta and
        % the null-vectors.
        h = alpha * n(:,1) + beta * n(:,2) + n(:,3);
        Hi = inv(T2) * reshape(h, 3, 3)' * T1;

        % Calculate the score of the homography
        pts2_t = Hi * source_points;
        pts2_t = rdivide(pts2_t, pts2_t(3,:));
        residuals = vecnorm(destination_points - pts2_t);
        
        % The inliers
        inliers = find(residuals < truncated_threshold);        
        % The MSAC score
        score = sum(1 - residuals(inliers) / truncated_threshold);
        
        % Update the so-far-the-best model if needed
        if score > best_score
            best_score = score;
            best_inliers = inliers;
            best_homography = Hi;
            
            % Updating the RANSAC iterationg number
            inlier_number = length(best_inliers);
            max_iterations = ...
                log(1 - confidence) / log(1 - (inlier_number / N)^2);
        end
    end 
    
    % Break if the iteration number has been exceeded
    if iter > max_iterations
        break;
    end
end

% Least-squares fitting on all the recovered inliers by the normalized DLT
% algorithm
H = normalizedDLT(M(best_inliers, 1:3), M(best_inliers, 4:6));

% Selecting the inliers and calculating the average error on them 
pts2_t = H * source_points;
pts2_t = rdivide(pts2_t, pts2_t(3,:));
residuals = vecnorm(destination_points - pts2_t);
best_inliers = find(residuals < truncated_threshold);    

disp(['Average re-reprojection error = ', num2str(mean(residuals(best_inliers))), ' px'])
disp(['Inlier number = ', int2str(length(best_inliers))])
disp(['Iteration number required for RANSAC = ', int2str(iter)])


% Draw points
close all;
image = [img1 img2];

imshow(image)
hold on;

for i = 1 : length(best_inliers)
    color = rand(1,3);
    
    plot([M(best_inliers(i), 1), size(img1, 2) + M(best_inliers(i), 4)], [M(best_inliers(i), 2) M(best_inliers(i), 5)], 'Color', color)
    
    scatter(M(best_inliers(i), 1), M(best_inliers(i), 2), 40, color, 'filled');
    scatter(size(img1, 2) + M(best_inliers(i), 4), M(best_inliers(i), 5), 40, color,'filled');
end
hold off;







