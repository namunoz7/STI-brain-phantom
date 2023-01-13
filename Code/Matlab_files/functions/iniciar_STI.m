function [principalFolder,imagesFolder] = iniciar_STI()
    
principalFolder = pwd;
f = filesep;
tmp = ['..',f,'..',f,'Imagenes',f,'STI'];
cd(tmp)
imagesFolder = pwd;
cd(principalFolder)
end