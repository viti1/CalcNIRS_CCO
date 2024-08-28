function [ conc, time_vector, fig , substanceNames] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile )
%% CalcNIRS - calculate and plot HbR HbO 
% Input:
%   dataFile - .mat file with intensity data.
%              it can have two possible formats, that will contain two "channels" 
%              (a) 'intensity1' and 'intensity2' fields
%              (b) 'cam1' and 'cam2' fields , each one has 'intensity' subfield 
%                
%  SDS  - Sourse-Detector Separation distance in cm. Value for each channel ( 2x1 vector corresponding to cam1 and cam2 )
%  tissueType - one of the rows in DPFperTissueFile (for example 'adult_forearm' \ 'baby_head' \ 'adult_head' \ 'adult_leg' )
%               default = 'adult_head'
%  plotChannelIdx - vector with numbers in the range of [1-2] indicating channels to plot. If empty - none is plotted. (default = [])
%    
%  extinctionCoefficientsFile - .csv file with the following columns : wavelength,	Water,	HbO2,	HHb,	FatSoybean, diffCCO, ox_redCCO, oxCCO, redCCO 
%                               default = '.\ExtinctionCoefficientsData.csv'    (if not passed or empty)                            
%  DPFperTissueFile - .txt file with two columns: Tissue and DPF  (Tissue is tissue type, corresponding to tissueType input variable)
%                     measured at 807nm 
%                     default = '.\DPFperTissue.txt'       (if not passed or empty)
%  relDPFfile - relative DPF according to wavelength
%              default = '.\RelativeDPFCoefficients.csv'   (if not passed or empty)     
%
% Output :
%   conc - cell array (cell for each channel) with concentration matrix ( rows : HbO2 HbR oxCCO redCCO)
%   time_vector - accorgin to the first wavelength in data file , in [seconds]
%   fig  - handle to figure. Empty if plotChannelIdx==[].
%
%% Check User Input
if ~exist('tissueType','var') || isempty(tissueType)
    tissueType = 'adult_head';
end

if ~exist('extinctionCoefficientsFile','var') || isempty(extinctionCoefficientsFile)
    extinctionCoefficientsFile = 'ExtinctionCoefficientsData.csv';
end

if ~exist(extinctionCoefficientsFile,'file') 
    error([ extinctionCoefficientsFile ' does not exist' ]);
else
    [~,~,ext] = fileparts(extinctionCoefficientsFile);
    if ~strcmp(ext,'.csv')
        error('extinctionCoefficientsFile should be a .csv file');
    end
end

if ~exist('DPFperTissueFile','var') || isempty(DPFperTissueFile)
    DPFperTissueFile = 'DPFperTissue.txt';
end
if ~exist(DPFperTissueFile,'file') 
    error([ DPFperTissueFile ' does not exist' ]);
else
    [~,~,ext] = fileparts(DPFperTissueFile);
    if ~strcmp(ext,'.txt')
        error('DPFperTissueFile should be a .txt file');
    end
end

if ~exist('relDPFfile','var') || isempty(relDPFfile)
    relDPFfile = 'RelativeDPFCoefficients.csv';
end
if ~exist(relDPFfile,'file') 
    error([ relDPFfile ' does not exist' ]);
else
    [~,~,ext] = fileparts(relDPFfile);
    if ~strcmp(ext,'.csv')
        error('relDPFfile should be a .csv file');
    end
end

% get record name
[~, recordName ] = fileparts(fileparts(dataFile));
%% Load Data and Coefficients
% 1. loda data
Data = load(dataFile);
nChannels =  nnz(startsWith( fieldnames(Data) , 'cam' ));
dataInfoFile = [fileparts(dataFile) '\info.mat'];
if exist(dataInfoFile,'file')
    info = load(dataInfoFile);
end

intensity_field = false;
if nChannels < 1
    nChannels =  nnz(startsWith( fieldnames(Data) , 'intensity' ));
    intensity_field = true;
    if nChannels < 1
        error('in data file expected at least one cam field');
    end
end
if ~exist('plotChannelIdx','var') 
    plotChannelIdx = 1:nChannels;
end
wavelengths = [785,808,830,860,680];

if numel(SDS) ~= nChannels
    error('You should provide the SDS for each channel')
end

% 2. Base DPF according to tissue
if isstring(tissueType); tissueType = char(tissueType); end
tblDPF = readtable(DPFperTissueFile);
DPFbase = tblDPF.DPF(ismember(tblDPF.Tissue,tissueType));
if isempty(DPFbase)
    error(['Unknown tissue type : "' tissueType '"']);
end
fid = fopen(DPFperTissueFile); line = fgetl(fid); fclose(fid);
wavelength_of_DPFdata = sscanf(line,'%%%*[^0123456789]%d');

% 3. Relative coefficient - depending on the wavelength
tblDPFrel = readtable(relDPFfile);
DPFrel = nan(numel(wavelengths),1);
for wi = 1:numel(wavelengths)
    DPFrel(wi) = tblDPFrel.relDPFcoeff(tblDPFrel.wavelength==wavelengths(wi) )/ tblDPFrel.relDPFcoeff(tblDPFrel.wavelength==wavelength_of_DPFdata);
end

% 4. total effective pathlength
Leff = DPFbase * DPFrel * SDS; % mattix with Leff for each wavelength(rows) for each channel(columns)

% 5. extinction coefficients
tblExtCoeff = readtable(extinctionCoefficientsFile);

%% Calc
  % init
  time_vector   = Data.time(1,:);
  
  % Create extinction coefficients matrix
  wavelength_idx = nan(size(wavelengths));
  for wi = 1:numel(wavelengths)
      wavelength_idx(wi) = find(tblExtCoeff.wavelength==wavelengths(wi) );      
  end
  
  substanceNames = {'HbO2','HHb','oxCCO','redCCO'};
  colors        = {'r','b',[0 0.8 0],'y'};

  extCoeffMat = []; 
  for mi = 1:numel(substanceNames)
        extCoeffMat = [ extCoeffMat , tblExtCoeff.(substanceNames{mi})(wavelength_idx)]; %#ok<AGROW>
  end

  idx_reftime = 1;      

  % loop over each channel    
  nOfTimePoints = size(Data.time,2) - idx_reftime + 1;
  conc = cell(1,nChannels);
  for ch_i = 1:nChannels
      if intensity_field
          raw_intensity = Data.(['intensity' num2str(ch_i)])';
      else
          raw_intensity = Data.(['cam' num2str(ch_i)]).intensity';
      end

      % 1. Calculate attenuation using raw intensity, and devide
      atten_div_Leff = nan(numel(wavelengths),nOfTimePoints);
      for wave_i = 1:numel(wavelengths)
          atten_div_Leff(wave_i,:) = log10(raw_intensity(idx_reftime,wave_i)./raw_intensity(:,wave_i)) .* 1/Leff(wave_i,ch_i);
      end

      % Calculate concentration  
      conc{ch_i} =  pinv(extCoeffMat) * atten_div_Leff ;
  end
  zero_idx = time_vector == 0;
  time_vector(zero_idx) = [];
  for ch=1:numel(conc)
      conc{ch}(:,zero_idx) = [];
  end
  
  save('DebugVika.mat'); %,'atten_div_Leff','extCoeffMat','conc','time_vector')

%% Plot
if ~isempty(plotChannelIdx)
    fig = figure('name',recordName,'Units','normalized','Position',[0.25      0.3      0.6   0.5]); 
    plot_idx = 1;
    for ch_i = plotChannelIdx(:)'
        subplot(numel(plotChannelIdx),1,plot_idx); plot_idx = plot_idx + 1;
        for mi = 1:numel(substanceNames)
            plot(time_vector/60,conc{ch_i}(mi,:)*1e6,'-','color',colors{mi}); hold on;
        end   
        xlabel('time [min]');
        ylabel('\Delta[\muM]');
        legend(substanceNames,'interpreter','none');
        title(sprintf('%s - SDS %g cm',recordName,SDS(ch_i)));
        grid on;

        ylims = get(gca,'YLim');

        if exist('info','var') && isfield(info,['DU_per_Pixel' num2str(ch_i)])
            text(5,ylims(1)+diff(ylims)*0.2,[ '<I> = ' num2str(round( info.(['DU_per_Pixel' num2str(ch_i)]))) ' DU'] )
        end
    end

    % plot timing of the events    
    events = ParseTiming(dataFile);
    if ~isempty(fieldnames(events))
        plot_idx = 1;
        for ch_i = plotChannelIdx(:)'
            subplot(numel(plotChannelIdx),1,plot_idx); plot_idx = plot_idx + 1;
            ylims = get(gca,'YLim');
            for k=1:numel(events)
                plot(events(k).min*[1 1], ylims,'-k');
                text(events(k).min+1,ylims(2)*0.8,events(k).name,'Rotation',-90,'Clipping','on','FontSize',12)
            end
            set(gca,'YLim',ylims);
            legend(substanceNames,'interpreter','none');        
        end
    end
else
    fig = [];
end

function events = ParseTiming(dataFile)
events = struct();
if exist([fileparts(dataFile) ,'\timing.txt'],'file')
    txt = fileread([fileparts(dataFile) ,'\timing.txt']);
    txt_split = strsplit(txt,'/');
    
    idx = 1;
    for k = 1:numel(txt_split)
        txt_split2 = strsplit(txt_split{k});
        if ~endsWith(txt_split2{end},'mins')
            if ~strcmp(txt_split2{end},'-')
                warning(['error parsing timing: ' txt_split{k}]);
            end
            continue;
        end
        events(idx).name = strjoin(txt_split2(1:end-1));
        events(idx).min = str2double(txt_split2{end}(1:end-4));
        if isnan(events(idx).min)
            warning(['error parsing timing: ' txt_split{k}]);
        end
        idx = idx + 1;
    end
end


