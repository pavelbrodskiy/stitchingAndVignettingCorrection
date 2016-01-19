% This application takes a folder with raw z-stacks and converts them into
% stitched confocal images. The files should have the default metamorph
% naming conventions, for example:
%
% FILENAME 00001_w1Confocal 405.TIF
% FILENAME 00001_w2Confocal 488.TIF
% FILENAME 00001_w3Confocal 561.TIF
% FILENAME 00001_w4Confocal 640.TIF
% 
% Inputs required in order to read the file will be a string array of
% channels in order (for example, here it would include 405, 488, 561, and
% 640), the filename immediately prior to the numbers (which would be
% 'FILENAME ' here), and the number of digits in the filenumber, as well
% as the index of the starting image, the dimensions of the montages and 
% the number of tiles.
%
% The order of the channels in the colors array determines which channel is
% represented by which color in the RGB image, in the order red, green,
% blue, white.
%
% This script z-projects each file, and saves the z-projections in a new
% folder, with the index and channel number. Then, if selected, a high-pass
% filter is applied to each image and the images are normalized. Finally,
% images are stitched with calculated overlaps with the assumption that the
% different channels in each tile should align with each other.

clear all
close all

%% Inputs
% Stitching Inputs
inputDirectory = 'E:\3D Stitching with vignetting correction\en CycE DCP\'; % The full directory of raw data
tempDirectory = [inputDirectory 'temp\']; % Temp directory
zprojDirectory = [inputDirectory 'zproj\']; % z-proj directory
outputDirectory = [inputDirectory 'stitched\']; % output directory
normalizationDirectory = [inputDirectory 'normalized\']; 
fftDirectory = [inputDirectory 'fft\']; 

fileName = 'en cycE DCP '; % FILENAME of input file
numChars = 5; % Number of characters in the file index (ei 00001)
tiles = 6*6; % Total number of tiles in folder
colors = {'640', '488', '561', '405'}; % Channels in order of red->green->blue->white
channels = {'405', '488', '561', '640'}; % Channels in order of name-number (ie w1, w2)
M = '6'; % Number of columns in montage
N = '6'; % Number of rows in montage
numembs = 1; % This is important, number of stitched embryos (should be tiles/M/N)
overlap = '30'; % Percentage of overlap in he tiles
outputNumber = 1; % Starting index of the final stitched image
outputName = '2015.02.08 cycE DCP '; % Name of the final stitched image
doVignettingCorrection = true; % Should we do FFT and normalization?
gaussSigma = 5; % If > 0, then sigma for gaussian filtering in 3D, else no filtering


% Vignetting Correction Inputs
if doVignettingCorrection
    x = 5; % X dimension of down-scaled overlap
    y = 5; % Y dimension of down-scaled overlap
    z = 5; % Z dimension of down-scaled overlap
    CropPercent = 0.01; % Percentage of image to crop away off edges
    name = 'Modified recursion\';
    chan = [405 488 561 640];
    start = 1; % Index of first tile
    filename = 'wt dpERK rbt ';
    
    
    % Parameters for the global optimization function
    lb = [1 0];
    ub = [200 5000];
    x0 = [2 200];
    MFE = 1e5;
    MI = 1e4;
    TolCon = 1e-3;
    TolFun = [1e-5 1e-2];
    TolX = 1e-5;
    time = 600;
    ST1 = 500;
    Trial = 2000;
end 

%% Initialize
addpath('C:\Fiji.app\scripts');
Miji(false);
mkdir(tempDirectory);
mkdir(outputDirectory);
mkdir(zprojDirectory);
mkdir(normalizationDirectory);
mkdir(fftDirectory);

%% Z-project and FFT
% for k = 1:tiles
%     for j = 1:length(channels)
%         path = [inputDirectory fileName num2str(k,['%0' num2str(numChars) 'd']) '_w' num2str(j) 'Confocal ' channels{j} '.TIF'];
%         savepath = [zprojDirectory fileName num2str(k,['%0' num2str(numChars) 'd']) '_w' num2str(j) 'Confocal ' channels{j} '_MIP.TIF'];
%         savepath2 = [fftDirectory fileName num2str(k,['%0' num2str(numChars) 'd']) '_w' num2str(j) 'Confocal ' channels{j} '_MIP.TIF'];
%     
%         MIJ.run('Open...', ['path=[' strrep(path, '\', '\\') ']']);
%         
%         if(doVignettingCorrection)
%             MIJ.run('Bandpass Filter...', 'filter_large=120 filter_small=0 suppress=None tolerance=5 process');
%         end
%         
%         MIJ.run('Tiff...', ['path=[' strrep(savepath2, '\', '\\') ']']);
%         
%         MIJ.run('Z Project...', 'projection=[Max Intensity]');
%         
%         MIJ.run('Tiff...', ['path=[' strrep(savepath, '\', '\\') ']']);
%         
%         MIJ.run('Close All');
%     end
% end

%% Normalization
if(doVignettingCorrection)
    aaaa1 = str2double(M);
    aaaa2 = str2double(N);
    aaa3 = str2double(overlap)/100;
    
[timeTakenout] = Global_Stitching_Optimization(aaaa1, aaaa2, x, y, z, aaa3, CropPercent, fftDirectory, lb, ub, x0, MFE, ... 
    MI, TolCon, TolFun, TolX, time, ST1, Trial, name, chan, numembs, start, fileName, normalizationDirectory, numChars, gaussSigma);
end

%% Combine Channels
for k = 1:tiles
    
    for j = 1:length(channels)
        path = [zprojDirectory '\' fileName num2str(k,['%0' num2str(numChars) 'd']) '_w' num2str(j) 'Confocal ' channels{j} '_MIP.TIF'];
        if(doVignettingCorrection)
        path = [normalizationDirectory '\' fileName num2str(k,['%0' num2str(numChars) 'd']) '_w' num2str(j) 'Confocal ' channels{j} '_MIP.TIF'];
        end
            MIJ.run('Open...', ['path=[' strrep(path, '\', '\\') ']']);
        MIJ.run('Duplicate...', ['title=' channels{j}]);
    end
    
    if length(colors) == 3
        MIJ.run('Merge Channels...', ['c1=[' colors{1} ']'...
            ' c2=[' colors{2} '] '...
            'c3=[' colors{3} ']'...
            ' create']);
    elseif length(colors) == 4
        MIJ.run('Merge Channels...', ['c1=[' colors{1} ']'...
            ' c2=[' colors{2} '] '...
            'c3=[' colors{3} ']'...
            'c4=[' colors{4} ']'...
            ' create']);
    else
        error('add more color options asshole');
    end
    
    MIJ.run('Tiff...', ['path=[' strrep([tempDirectory 'Tile ' num2str(k,'%05d') '.TIF]'], '\', '\\')]);
    MIJ.run('Close All');
    
end

%% Stitch

for i = 1:(str2double(M)*str2double(N)):(tiles-1)

MIJ.run('Grid/Collection stitching', ['type=[Grid: row-by-row] order=[Right & Down     ' ...
    '           ] grid_size_x=' M ' grid_size_y=' N ' tile_overlap=33 first_file_index_i=' num2str(i) ...
    ' directory=[' strrep(tempDirectory, '\', '\\') '] file_names=[Tile {iiiii}.tif] ' ... 
    'output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] ' ... 
    'regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement' ... 
    '_threshold=3.50 compute_overlap computation_parameters=[Save memory (but be slower)]' ... 
    ' image_output=[Fuse and display]']);

MIJ.run('Tiff...', ['path=[' strrep([outputDirectory outputName num2str(outputNumber,'%03d')], '\', '\\') '.TIF]']);
%MIJ.run('Tiff...', ['path=[' strrep([tempDirectory '\' fileName aaa '_w' num2str(j) 'Confocal ' channels{j} '_MIP'], '\', '\\') '.TIF]']);
        
MIJ.run('Close All');
    
    outputNumber = outputNumber + 1;
end

%% Cleanup
MIJ.exit;
rmdir(tempDirectory);
