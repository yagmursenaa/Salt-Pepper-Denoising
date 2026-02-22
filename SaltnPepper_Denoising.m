A = imread('peppers.png');

class(A)
size(A)
ndims(A)
A = double(A);

if ndims(A) == 3
    R = A(:,:,1);
    G = A(:,:,2);
    B = A(:,:,3);
    I = 0.299*R+0.587*G+0.114*B;
else
    I =A;
end

figure;
imshow(uint8(I))
title('grayscale image')

snp_density=0.1;

[m, n]= size(I);
snp_img =I;                     
snp_noise= rand(m, n); % uniform distribution

snp_img(snp_noise <snp_density/2)= 0; 
snp_img(snp_noise>=snp_density/2 & snp_noise< snp_density)= 255; 

figure;
imshow(uint8(snp_img)) % display the image
title('image with noise ')

extended_img = zeros(m+2 , n+2);
extended_img(2:m+1,2:n+1) = snp_img;


%Median filter
med_filteredimg=snp_img;   

for i = 2:m+1
    for j =2:n+1

        window = [ extended_img(i-1,j-1), extended_img(i-1, j), extended_img(i-1,j+1);
                   extended_img(i,j-1), extended_img(i, j), extended_img(i , j+1);
                   extended_img(i+1, j-1), extended_img(i+1,j), extended_img(i+1, j+1)];

        windowvec = window(:);   % 9x1 vectör
        med_filteredimg(i-1, j-1) = median(windowvec); 
    end
end
figure;
imshow(uint8(med_filteredimg))
title('with median filter')


% Min filter(applied after the median filter)
min_filteredimg =zeros(m,n);
extended_med= zeros(m+2,n+2);
extended_med(2:m+1,2:n+1) =med_filteredimg;

for i =2:m+1
    for j = 2:n+1
        window = [extended_med(i-1,j-1), extended_med(i-1, j), extended_med(i-1,j+1);
                  extended_med(i,j-1), extended_med(i, j), extended_med(i , j+1);
                  extended_med(i+1, j-1), extended_med(i+1,j), extended_med(i+1, j+1) ];
       
        windowvec= window(:);
        min_filteredimg(i-1, j-1) =min(windowvec);
    end
end
figure;
imshow(uint8(min_filteredimg))
title('median+ min filter')


% PSNR 
mse_med = mean((I(:) -med_filteredimg(:)).^2);
psnr_med = 10* log10(255^2/ mse_med);
 
mse_med_min = mean((I(:) -min_filteredimg(:)).^2);
psnr_med_min = 10* log10(255^2/ mse_med_min);

disp(['PSNR (Median Filter): ', num2str(psnr_med)])
disp(['PSNR (Median+ Min Filter): ', num2str(psnr_med_min)])

% Global SSIM
L = 255; % dynamic range
C1 = (0.01 *L)^2;
C2 = (0.03 *L)^2;

mean_I =mean(I(:));
var_I = var(I(:));

mean_med = mean(med_filteredimg(:));
var_med = var(med_filteredimg(:));

cov_mat1 = cov(I(:), med_filteredimg(:));
cov1 = cov_mat1(1,2);

ssim1 =(2*mean_I*mean_med+C1)*(2*cov1 + C2) /((mean_I^2 +mean_med^2 + C1)*(var_I+var_med+C2));

mean_med_min =mean(min_filteredimg(:));
var_med_min =var(min_filteredimg(:));

cov_mat2 = cov(I(:), min_filteredimg(:));
cov2 = cov_mat2(1,2);

ssim2 = (2*mean_I*mean_med_min + C1)*(2*cov2 + C2) / ((mean_I^2 + mean_med_min^2 + C1) * (var_I + var_med_min + C2));

disp(['SSIM (Median Filter): ', num2str(ssim1)])
disp(['SSIM (Median + Min Filter): ', num2str(ssim2)])

