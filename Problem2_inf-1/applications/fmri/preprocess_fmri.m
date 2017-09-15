clc;
clear;
close all;

load('/home/gustavo/Downloads/data-science-P1.mat')

fid=fopen('/home/gustavo/Gdrive/Papers/Mixed norm/lista_palabras');
tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

num_voxels = length(data{1});
num_words = length(tlines);
mris = zeros(num_voxels,num_words);

for ww = 1: length(tlines)
    num_tomas(ww) = 0;
    for ii=1:(length(data))
        if strcmp(info(ii).word,char(tlines{ww}))
            num_tomas(ww)=num_tomas(ww)+1;
            mris(:,ww)=mris(:,ww)+data{ii}';
        end
    end
    mris(:,ww)=mris(:,ww)/num_tomas(ww);
end

save('/home/gustavo/Gdrive/Papers/Mixed norm/mris.mat','mris')