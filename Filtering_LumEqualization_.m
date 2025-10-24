close all
clear all

%Stimulus and Save Paths
Stimpath = '\exemple_visages\'; %path of original images
Savepath = '\exemple_visages\'; %path for saving filtered images
fmt = '.png'; %images format 


% The cut-off frequency is specified in cycle per image, so it is necessary to know the size of the
% image and the desired cut-off frequency in cycle per degree 
% (e.g., for images of  1024 pixels = 24 degrees 
% -> 1cpd = 1* 24 degrees per image = 24 cycles per image
% -> 2cpd = 2* 24 degrees per image = 48 cycles per image

fc1h=4.6;   % LSF cut-off frequency in cycles per image
fc2h=27.6;   % HSF cut-off frequency in cycles per image

% Definition of mean luminance between 0 and 1
lummoyenne = 0.5;


%list image files
files = dir([Stimpath '/*' fmt]);


for i=1:length(files)
    % if images are in .png and contain a transparent channel (i.e. oval
    % mask)
     if strcmp(fmt,'.png')
        [im map transparency] = imread([Stimpath '/' files(i).name]); %read image
        im = double(im)/255; %convert it to color levels between 0 and 1
        transparency = double(transparency)/255; %get oval mask
        graybckg = (ones(size(transparency))-transparency)/2; %define area for gray background
        im = im+transparency/2; %put the face on a gray bakcground
                
    else %if images are not in png 
    im=double(imread([Stimpath '/' files(i).name]))/255;
     end
    
    %convert to graysacale
    [H W P] = size(im);
    if size(im,3)==3
        im = sum(im,3)/3; 
    end
    
    %check mean luminance and RMS contrast of image
    mean(im(:))
    std(im(:))
        
    fc1w=W*fc1h/H;
    fc2w=W*fc2h/H;
    
    % attenuation of signal (e.g., 3dB)
    A=1/2;

    % Calculate sigma for each frequency cut-off 
    sigma1h=fc1h/sqrt(-log(A));
    sigma2h=fc2h/sqrt(-log(A));
    sigma1w=fc1w/sqrt(-log(A));
    sigma2w=fc2w/sqrt(-log(A));

    % create frequency grid
    [fx,fy]=meshgrid(0:W-1,0:H-1);fx=fx-(W-1)/2;fy=fy-(H-1)/2;

    % build gaussian filters for low-pass and high-pass filters 
    g1=exp(-((fx/sigma1w).^2+(fy/sigma1h).^2));
    g2=1 - exp(-((fx/sigma2w).^2+(fy/sigma2h).^2));
    
    % multiply the gaussian filtees with the amplitude spectrum 
    Im_bf=real(ifft2(fftshift(fftshift(fft2(im)).*g1)));
    Im_hf=real(ifft2(fftshift(fftshift(fft2(im)).*g2)));

    % normalization if LSF and HSF images 
    im=(im-min(im(:)))/(max(im(:))-min(im(:)));
    Im_bf=(Im_bf-min(Im_bf(:)))/(max(Im_bf(:))-min(Im_bf(:)));
    Im_hf=(Im_hf-min(Im_hf(:)))/(max(Im_hf(:))-min(Im_hf(:)));
    
    fprintf('Min luminance of original image : %d\n',min(im(:)));
    fprintf('Min luminance of LSF image : %d\n',min(Im_bf(:)));
    fprintf('Min luminance of HSF image : %d\n',min(Im_hf(:)));
    fprintf('Max luminance of original image: %d\n',max(im(:)));
    fprintf('Max luminance of LSF image  : %d\n',max(Im_bf(:)));
    fprintf('Max luminance of HSF image : %d\n',max(Im_hf(:)));


% stardardization of luminance values

Im_c = (im - mean(im(:)));
Im_bf_c =(Im_bf - mean(Im_bf(:)));
Im_hf_c =(Im_hf - mean(Im_hf(:)));

fprintf('Mean luminance of original image : %d\n',mean(Im_c(:)));
fprintf('Mean luminance of LSF image : %d\n',mean(Im_bf_c(:)));
fprintf('Mean luminance of HSF image : %d\n',mean(Im_hf_c(:)));
fprintf('SD of original image : %d\n',std(Im_c(:)));
fprintf('SD of LSF image : %d\n',std(Im_bf_c(:)));
fprintf('SD of HSF image : %d\n',std(Im_hf_c(:)));

%Luminance equalization
Im_n = Im_c + lummoyenne;
Im_bf_n = Im_bf_c + lummoyenne;
Im_hf_n = Im_hf_c + lummoyenne;



if ~exist(Savepath,'dir')
     mkdir('.',Savepath)
end

%show and save filtered images
%figure
%imshow(uint8(255*Im_bf_n));colormap(gray(256))
imwrite(uint8(255*Im_bf_n),[Savepath '/' files(i).name(1:(length(files(i).name)-4)) '_BFS_' num2str(fc1h) 'cpi_passe_bas_LUM' fmt],'Alpha',transparency);

%figure
%imshow(uint8(255*Im_hf_n));colormap(gray(256))
imwrite(uint8(255*Im_hf_n),[Savepath '/' files(i).name(1:(length(files(i).name)-4)) '_HFS_' num2str(fc2h) 'cpi_passe_haut_LUM' fmt],'Alpha',transparency);


end