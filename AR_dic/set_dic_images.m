%set dic images
%===========================
% images
bname = 'AR';

wdir = strcat('/home/ga288/crack_healing/AR_dic/')


imagepath = strcat('/home/ga288/crack_healing/',bname)
imagefiles = dir(fullfile(imagepath,'*.jpg'));% .bmp or .tif 
imagefiles = {imagefiles.name}';
images = strcat(imagepath,filesep,imagefiles);

save([wdir,strcat('DIC_images_',bname,'.mat')],'images')
