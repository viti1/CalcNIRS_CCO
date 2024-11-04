mFolder = fileparts(mfilename('fullpath'));

recordsFolder = [ mFolder '\..\..\Records\UK_EXPERIMENTS_2024\Laser phantom data\']; 
recordsDir = dir([ recordsFolder 'Laser*']);
records = strcat(recordsFolder, {recordsDir.name})';

%%
prefix = '\results2\oxCCO_redCCO_';
figs = nan(size(numel(records)));
for ri = 1:numel(records)
    [conc, time_vector,figs(ri),substanceNames] = CalcNIRS([records{ri} '\results.mat'], [3,1],'adult_head',1:2);
    [~,recName ]  = fileparts(records{ri});
    savefig(figs(ri),[mFolder prefix recName '.fig']);
    saveas(figs(ri),[mFolder prefix recName '.png']);
    save([ mFolder prefix recName '.mat' ],'conc', 'time_vector','substanceNames');
end

%% save to PPT
import mlreportgen.ppt.*

% Define the directory containing your saved figures
figureDir = [ mFolder fileparts(prefix) '\' ]; % Current directory

% Start PowerPoint
ppt = Presentation([ figureDir '/NIRSwithYeast_oxCCO_redCCO.pptx' ]);
open(ppt);

% List of figure files
figureFiles = dir(fullfile(figureDir, '*.png'));

% Loop through each figure file and add it to the presentation
for i = 1:length(figureFiles)
    % Add a new slide
    slide = add(ppt,'Title and Picture');
    % Insert the figure
    imagePath = fullfile(figureDir, figureFiles(i).name);
    %replace(pictureSlide,'Title','Histogram of Vector');
    replace(slide,'Picture',Picture(imagePath));    
end

% Close PowerPoint
close(ppt);
rptview(ppt);