close all
clear all
%spécification des répertoires
Stimpath = 'C:\Users\louis\Nextcloud\PROJECTS\Current\Collab\B37T_PaulineFavre\exemple_visages\'; %chemin des images originales
Savepath = 'C:\Users\louis\Nextcloud\PROJECTS\Current\Collab\B37T_PaulineFavre\exemple_visages\'; %chemin où enregistrer les images filtrées
fmt = '.png' %format des images

% On spécifie la fréquence coupures en cycle par image il faut donc connaitre la taille de
% l'image et les fréquences de coupures voulues en cycle par degrés 
% (e.g., dans mes manips: images de  1024 pixels = 24 degrés 
% -> 1cpd = 1* 24 degrés par images = 24 cycles par image
% -> 2cpd = 2* 24 degrés par images = 48 cycles par image

%à changer
fc1h=4.6;   % fréquence de coupure BF en cycle par image
fc2h=27.6;   % fréquence de coupure HF en cycle par image

% définition de la luminance moyenne pour des niveaux de gris de 0 (noir) à
% 1 (blanc)
lummoyenne = 0.5;


%lecture des noms de fichiers des images
files = dir([Stimpath '/*' fmt]);


for i=1%:length(files)
    %si les images sont en png et contiennent un canal transparent (= masque
    %oval)
    if strcmp(fmt,'.png')
        [im map transparency] = imread([Stimpath '/' files(i).name]);
        im = double(im)/255;
        transparency = double(transparency)/255;
        graybckg = (ones(size(transparency))-transparency)/2;
        im = im+transparency/2;
                
    else %si les images ne sont pas au format png avec un canal transparent pour le masque il faut utiliser des images où levisage est collé sur un fond gris
    im=double(imread([Stimpath '/' files(i).name]))/255;
    end
    
    [H W P] = size(im);
    if size(im,3)==3
        im = sum(im,3)/3;  %conversion couleur nb
    end

    mean(im(:))
    std(im(:))
    
    
    fc1w=W*fc1h/H;
    fc2w=W*fc2h/H;
    % choix de l'atténuation A du signal (par exemple, 3dB)
    A=1/2;

    % calcul de sigma pour une fréquence de coupure passe-bas fc1 et une
    % fréquence de coupure passe-haut fc2
    sigma1h=fc1h/sqrt(-log(A));
    sigma2h=fc2h/sqrt(-log(A));
    sigma1w=fc1w/sqrt(-log(A));
    sigma2w=fc2w/sqrt(-log(A));

    % création de la grille des fréquences
    [fx,fy]=meshgrid(0:W-1,0:H-1);fx=fx-(W-1)/2;fy=fy-(H-1)/2;

    % constuction du filtre gaussien g1 pour le filtrage passe-bas et g2 
    %pour le filtrage passe-haut
    g1=exp(-((fx/sigma1w).^2+(fy/sigma1h).^2));
    g2=1 - exp(-((fx/sigma2w).^2+(fy/sigma2h).^2));
    
    % multiplication du filtre au spectre de fréquences spatiales de l'image
    Im_bf=real(ifft2(fftshift(fftshift(fft2(im)).*g1)));
    Im_hf=real(ifft2(fftshift(fftshift(fft2(im)).*g2)));

    % normalisation de l'image hf et de l'image bf
    im=(im-min(im(:)))/(max(im(:))-min(im(:)));
    Im_bf=(Im_bf-min(Im_bf(:)))/(max(Im_bf(:))-min(Im_bf(:)));
    Im_hf=(Im_hf-min(Im_hf(:)))/(max(Im_hf(:))-min(Im_hf(:)));
    
    fprintf('Min de l''image originale : %d\n',min(im(:)));
    fprintf('Min de l''image BF : %d\n',min(Im_bf(:)));
    fprintf('Min de l''image HF : %d\n',min(Im_hf(:)));
    fprintf('Max de l''image originale : %d\n',max(im(:)));
    fprintf('Max de l''image BF : %d\n',max(Im_bf(:)));
    fprintf('Max de l''image HF : %d\n',max(Im_hf(:)));


% Egalisation de la luminance moyenne à 0
% et de la variance à 1
% Autrement dit centrage réduction

Im_c = (im - mean(im(:)));
Im_bf_c =(Im_bf - mean(Im_bf(:)));
Im_hf_c =(Im_hf - mean(Im_hf(:)));

fprintf('Moy de l''image originale : %d\n',mean(Im_c(:)));
fprintf('Moy de l''image BF : %d\n',mean(Im_bf_c(:)));
fprintf('Moy de l''image HF : %d\n',mean(Im_hf_c(:)));
fprintf('Std de l''image originale : %d\n',std(Im_c(:)));
fprintf('Std de l''image BF : %d\n',std(Im_bf_c(:)));
fprintf('Std de l''image HF : %d\n',std(Im_hf_c(:)));

%Egalisation de la luminance


Im_n = Im_c + lummoyenne;
Im_bf_n = Im_bf_c + lummoyenne;
Im_hf_n = Im_hf_c + lummoyenne;


% affichage des résultats pour vérification
fprintf('Moy de l''image NFn : %d\n',mean(Im_n(:)));
fprintf('Moy de l''image BFn : %d\n',mean(Im_bf_n(:)));
fprintf('Moy de l''image HFn : %d\n',mean(Im_hf_n(:)));
fprintf('Std de l''image NFn : %d\n',std(Im_n(:)));
fprintf('Std de l''image BFn : %d\n',std(Im_bf_n(:)));
fprintf('Std de l''image HFn : %d\n',std(Im_hf_n(:)));


if ~exist(Savepath,'dir')
     mkdir('.',Savepath)
end


%figure
%imshow(uint8(255*Im_bf_n));colormap(gray(256))
imwrite(uint8(255*Im_bf_n),[Savepath '/' files(i).name(1:(length(files(i).name)-4)) '_BFS_' num2str(fc1h) 'cpi_passe_bas_LUM' fmt],'Alpha',transparency);

%figure
%imshow(uint8(255*Im_hf_n));colormap(gray(256))
imwrite(uint8(255*Im_hf_n),[Savepath '/' files(i).name(1:(length(files(i).name)-4)) '_HFS_' num2str(fc2h) 'cpi_passe_haut_LUM' fmt],'Alpha',transparency);


end