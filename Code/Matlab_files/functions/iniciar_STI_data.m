function [principalFolder,imagesFolder] = iniciar_STI_data(path_to_folders)
    
principalFolder = pwd;
% f = filesep;
% path_to_folders = ['..',f,'..',f,'Imagenes',f,'STI'];
cd(path_to_folders)
imagesFolder = pwd;
cd(principalFolder)
end