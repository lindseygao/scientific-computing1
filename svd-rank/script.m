clear all; close all; clc
% Homework 7

% Motion of a particle
load('particle_position.mat'); % particle position data
% 'A' matrix (x,y,z) positions of the particle from t = 0 to 100: 
    % [x(0) x(1) ... x(100)
    %  y(0) y(1) ... y(100)
    %  z(0) z(1) ... z(100)]
for row = 1:3 % ^change A to above format
    A(row,:) = A(row,:) - sum(A(row,:))/101;
end
% Compute the SVD of the matrix A
[U,S,V] = svd(A);
ans1 = diag(S);
Arank1 = S(1,1) * U(:,1) * V(:,1)'; % rank-1 approximation
error1 = norm(A-Arank1);
ans2 = error1
Arank2 = U(:,1:2)*S(1:2,1:2)*V(:,1:2)'; % rank-2 approximation
error2 = norm(A-Arank2);
ans3 = error2

%Problem 2 
% a) Plot the positions of the particles & rank-1 approximation
plot3(A(1,:),A(2,:),A(3,:), 'k.', 'Linewidth', [4]);, hold on
plot3(Arank1(1,:),Arank1(2,:),Arank1(3,:), 'r.', 'Linewidth', [4]);
legend('Position of Particles', 'Rank-1 Approximation', 'Location',...
       'Best', 'Fontsize', [12])
xlabel('x');
ylabel('y');
zlabel('z');
print('Problem2a','-dpng')
hold off
% b) Plot the positions of the particles & rank-2 approximation
plot3(A(1,:),A(2,:),A(3,:), 'k.', 'Linewidth', [3]), hold on
plot3(Arank2(1,:),Arank2(2,:),Arank2(3,:), 'r.')
xlabel('x');
ylabel('y');
zlabel('z');
legend('Position of Particles', 'Rank-2 Approximation', 'Location',...
       'Best', 'Fontsize', [12])
print('Problem2b','-dpng')

%% Problem 3 Image Compression using Singular Value Decomposition
% a) Perform SVD
A = imread('llama.jpg'); % MATLAB built-in image from the image processing toolbox
A = im2double(rgb2gray(A));
[U,S,V] = svd(A, 'econ');
Sing = diag(S);
ans4 = Sing(1:20,1);
% b) Calculate the proportion of total energy contained in rank-1
% approximation
ans5 = ans4(1)/sum(Sing);
% c) Proportion of total energy contained in rank-20 approx
ans6 = sum(ans4)/sum(Sing); % energy of rank-20 approx
% d) Find the smallest value of r such that the energy of the rank-r
% approximation is >= 0.9
for r = 1:length(Sing)
    energy = sum(Sing(1:r,1))/sum(Sing);
    if energy >= 0.9
        break
    end
end
ans7 = r;

% Plotting images from matrix above ('A') - original, rank-1, rank-20,
% rank-r images
subplot(2,2,1), imshow(A);
title('Original Image')
r = 1;
Approx = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
subplot(2,2,2), imshow(Approx);
title('Rank-1 Approximation')
hold on 
r = 20;
Approx = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
subplot(2,2,3), imshow(Approx);
title('Rank-20 Approximation')
r = ans7
Approx = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
subplot(2,2,4), imshow(Approx);
title('Rank-206 Approximation')
print 'problem4' -dpng
% b) calculate the total number of pixels 
pixels = 876*1314 % for full image
rank_r = (876*r) + r + (1314*r) % for rank-r approximation

%% Image denoising: goal = denoise the noisy image using a low rank approximation
load('NoisyImage.mat') % noisy image to denoise; in var 'A_noise'
% SVD on noisy image
[U,S,V] = svd(A_noise);
sing_vals = diag(S);
% Computing energy of rank-2 approximation
energy_2 = sum(sing_vals(1:2,1))/sum(sing_vals);
ans8 = energy_2
r = 2;
% Computing the rank-2 approximation of A_noise 
Approx = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
ans9 = norm(A_noise-Approx) % error of rank 2 approx & A_noise
ans10 = norm(A-Approx) % error between rank2 & A

% Problem 6
% a) Plotting singular values of noisy matrix (A_noise) on log-y scale
semilogy(sing_vals, 'ob')
title('Singular Values of Noisy Image on Logarithmic Scale')
ylabel('Singular Value (log scale)')
print problem6a.png -dpng
% b) Plotting the true image, noisy image, & rank-2 approximation
subplot(1,3,1), imshow(A)
title('True Image')
subplot(1,3,2), imshow(A_noise)
title('Noisy Image')
subplot(1,3,3), imshow(Approx)
title('Rank-2 Approx')
print problem6b.png -dpng
r = 1;
Approx = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
imshow(Approx)
title('Rank-1 Approximation')
print problem6c.png -dpng


