clear all
close all
clc
%
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',18);
set(0,'defaultPatchEdgeColor',[0.5,0.5,0.5]);
set(0,'defaultPatchFaceColor','none');
set(0,'defaultFigureColormap',parula);
%
flag_mes_cur = 2;
file_save = ['DI_results_mesh',num2str(flag_mes_cur),'_D.mat'];
eval(['load(''',file_save,''')'])
radix = strrep(file_save,'_D.mat','');
file_strain = [radix,'_strains.mat'];
%
Udic = Mesh.Udic;
Rdic = Mesh.R;
Xdic = Mesh.coordinates(:,1);
Ydic = Mesh.coordinates(:,2);
Odic = ones(size(Ydic));
%
imagepath = '../../AR';
imagefiles = dir(fullfile(imagepath,'*.jpg'));% .bmp or .tif //uncompressed
imagefiles = {imagefiles.name}';
images = strcat(imagepath,filesep,imagefiles);

frgb = imread(images{1});
[n,m,~] = size(frgb);

% frgb(:,:,1) = fz;
% frgb(:,:,2) = fz;
% frgb(:,:,3) = fz;
% fz = single(fz);

dec  = 1;
imax = size(Rdic,2);
%
if ~exist(file_strain,'file')
    Edic = mesh_strain(Mesh,Mesh.Udic,'StrainDefinition','GreenLagrange');
    eval(['save(''',radix,'_strains.mat'',''Edic'')'])
else
    eval(['load(''',file_strain,''')'])
end
%
%h = figure;
co = 0;
for idx = 1:dec:imax
    co = co+1;
    h = figure;
    box on
    set(h,'Position',[350 150 1000 600])
%     maxe = 5*std(Edic{idx}(:,2));
%     ind = find(abs(Edic{idx}(:,2))>maxe);
%     Edic{idx}(ind,2) = NaN;

    grayImage = rgb2gray(imread(images{co}));
    frgb  = cat(3, grayImage, grayImage, grayImage);%ind2rgb(grayImage, 'gray');

    imagesc(frgb)
    axis image
    hold on
    correli_fig(Mesh,Udic{idx},Edic{idx}(:,2),m,n)
    hold off
    alpha(0.6)
    title(sprintf(' # %03d',idx))
    %pause(1)
    a=colorbar;
    ylabel(a,'$\varepsilon_{xx}$','Interpreter','latex','FontSize',26);%,'Rotation',270);
    caxis([-0.1, 0.1]);
    
    frame = getframe(h);
    img{co} = frame2im(frame);
    close(h)
end
%
file_print = strrep(file_save,'.mat','_exx_ref.gif');
for idx = 1:co
    [A,map] = rgb2ind(img{idx},256);
    if idx == 1
        imwrite(A,map,file_print,'gif','LoopCount',Inf,'DelayTime',.25);
    else
        imwrite(A,map,file_print,'gif','WriteMode','append','DelayTime',.25);
    end
end