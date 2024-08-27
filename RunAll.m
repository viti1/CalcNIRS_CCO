mFolder = fileparts(mfilename('fullpath'));

recordsFolder = [ mFolder '\..\..\Records\UK_EXPERIMENTS_2024\Laser phantom data\']; 
recordsDir = dir([ recordsFolder 'Laser*']);
records = strcat(recordsFolder, {recordsDir.name})';


figs = nan(size(numel(records)));
for ri = 1:numel(records)
    [conc, timeVec,figs(ri)] = CalcNIRS([records{ri} '\results.mat'], [3,1],'adult_head',1:2);
    [~,recName ]  = fileparts(records{ri});
    savefig(figs(ri),[mFolder '\figs\oxCCO_' recName '.fig']);
    saveas(figs(ri),[mFolder '\figs\oxCCO_' recName '.png']);
    save([ mFolder '\figs\oxCCO_' recName '.mat' ],'conc', 'timeVec');
end

%% save to PPT
import mlreportgen.ppt.*

% Define the directory containing your saved figures
figureDir = [ mFolder '/figs/' ]; % Current directory

% Start PowerPoint
ppt = Presentation([ figureDir '/NIRSwithYeast_1.pptx' ]);
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