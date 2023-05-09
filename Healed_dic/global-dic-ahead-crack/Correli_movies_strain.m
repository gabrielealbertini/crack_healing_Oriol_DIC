clear all
close all
clc
% specify the path to the Correli folder
Correlipath = '/home/ga288/sources/correli/'; %fullfile(pwd,'..','..','..');

% this adds the correct folder to the matlab path
addpath(fullfile(Correlipath,'ml_lib_core'));

%

set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',10);
set(0,'DefaultAxesLabelFontSize',1.7);
%set(0,'DefaultLabelFontName','Arial','DefaultLabelFontSize',12);
set(0, 'defaultTextInterpreter', 'latex'); 
%set(0, 'defaultLabelTextInterpreter', 'latex'); 

set(0,'defaultPatchEdgeColor',[0.5,0.5,0.5]);
set(0,'defaultPatchFaceColor','none');
set(0,'defaultFigureColormap',parula);
%
flag_mes_cur = 4;
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
imagepath = '../../Healed/';
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
istart=1;
imax=28;
for idx = istart:dec:imax
    co = co+1;
    h = figure;
   
    set(h,'Position',[350 150 1000 600])    
    box on

    %     maxe = 5*std(Edic{idx}(:,2));
    %     ind = find(abs(Edic{idx}(:,2))>maxe);
    %     Edic{idx}(ind,2) = NaN;

    grayImage = rgb2gray(imread(images{co}));
    frgb  = cat(3, grayImage, grayImage, grayImage);%ind2rgb(grayImage, 'gray');

    imagesc(frgb)
    axis image
    hold on
    correli_fig(Mesh,Udic{idx},Edic{idx}(:,2),m,n)
    xlim([800,1800])
    ylim([300,650])
    hold off
    %alpha(0.75)
    title(sprintf('Healed # %03d',idx),'fontsize',16)
    %pause(1)
    set(gca,"Position",[0.1 0.2 0.7 0.7]); 
    a=colorbar;
    ylabel(a,'$\varepsilon_{yy}$','Interpreter','latex','fontsize',18);

    
    caxis([0, 0.05]);
    figdir=['DI_results_mesh',num2str(flag_mes_cur,'%01.f'),'_D'];
    mkdir(figdir);
    
    fname = [figdir,'/Healed_eyy',num2str(idx,'%03.f'),'.png'];
    saveas(h,fname)
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
