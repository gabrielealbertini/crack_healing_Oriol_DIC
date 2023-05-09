%set dic images
%===========================
% images
bname = 'Healed';

wdir = strcat('/home/ga288/crack_healing/Healed_dic/')


imagepath = strcat('/home/ga288/crack_healing/',bname)
imagefiles = dir(fullfile(imagepath,'*.jpg'));% .bmp or .tif 
imagefiles = {imagefiles.name}';
images = strcat(imagepath,filesep,imagefiles);

save([wdir,strcat('DIC_images_',bname,'.mat')],'images')
