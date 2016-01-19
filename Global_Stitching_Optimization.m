function [timeToFinish] = Global_Stitching_Optimization(M, N, x, y, OverlapPercent, CropPercent, input, lb1, ub1, x01, MFE1, MI1, ...
    TolCon1, TolFun1, TolX1, time1, ST11, Trial1, name, chan, numembs, start, filename, output, numChars)

%% Initialize Global Normalization

tic
%[perfectTiles, damagedTiles, rawTiles] = testtesttest(M,N);

warning('off','all');
global lb ub x0 MFE MI TolCon TolFun TolX time ST1 Trial
    
%mkdir(input,name);
channum = max(size(chan));
embnumel = int2str(numel(int2str(numembs)));

lb = lb1;
ub = ub1;
x0 = x01;
MFE = MFE1;
MI = MI1;
TolCon = TolCon1;
TolFun = TolFun1;
TolX = TolX1;
time = time1;
ST1 = ST11;
Trial = Trial1;

%% Input images

for m = channum:-1:1;
    t = 1;
    imageNumber = start;
    
    for k = 1:numembs
        %% Read images
        
        disp(['Embryo: ' int2str(k) ' Channel: ' int2str(m)]);
        for j = 1:N
            for i = 1:M
                rawTiles(i,j,:,:) = imread([input filename num2str(imageNumber,['%0' num2str(numChars) 'd']) '_w' num2str(m) 'Confocal ' num2str(chan(m)) '_MIP.TIF']);
                imageNumber = imageNumber + 1;
            end
        end
        
% rawTiles(1:M,1:N,1:512,1:512) = 0;
% 
% for i = 1:512
%     for j = 1:512
%         vigTemplate(i,j) = cos(sqrt((256-i)^2+(256-j)^2)/300);
%     end
% end
% imshow(vigTemplate,[])
% MIJ.createImage('result', vigTemplate, true);
% MIJ.run('Bandpass Filter...', 'filter_large=120 filter_small=0 suppress=None tolerance=5');
% MIJ.getImage('result');
% MIJ.closeAllWindows;
% 
% for i = 1:M
%     for j = 1:N
%         tileTempGen = uint8(vigTemplate .* normrnd(200,40,512));
%         a = fft2(tileTempGen);
%         a(120:end, 120:end) = 0;
%         rawTiles(M,N,:,:) = gaussianbpf(tileTempGen,0,120);
%     end
% end

        %% Crop Images
        
        [~, ~, totalxpixels, totalypixels] = size(rawTiles);
        
        xCrop = round(totalxpixels * CropPercent);
        yCrop = round(totalypixels * CropPercent);
        xOverlapPixels = round(OverlapPercent * totalxpixels - xCrop * 2);
        yOverlapPixels = round(OverlapPercent * totalypixels - yCrop * 2);
        cropxpixels = round(totalxpixels * (1-CropPercent*2));
        cropypixels = round(totalypixels * (1-CropPercent*2));
        cropTiles = rawTiles(:,:,(xCrop+1):(totalxpixels - xCrop),(yCrop+1):(totalypixels-yCrop));
        
        %% Obtain slice average overlap intensity (P)
        
        E = zeros(M,N,x,y,4);
        
        for i = 1:M
            for j = 1:N
                
                clear tempE1 tempE2 tempE3 tempE4;
                
                tempE1(:,:) = cropTiles(i,j,1:xOverlapPixels,:);
                tempE2(:,:) = cropTiles(i,j,:,(cropypixels - yOverlapPixels):cropypixels);
                tempE3(:,:) = cropTiles(i,j,(cropxpixels - xOverlapPixels):cropxpixels,:);
                tempE4(:,:) = cropTiles(i,j,:,1:yOverlapPixels);
                
                tempE2 = tempE2';
                tempE4 = tempE4';                
                
                E(i,j,:,:,1) = imresize(tempE1,[x,y]);
                E(i,j,:,:,2) = imresize(tempE2,[x,y]);
                E(i,j,:,:,3) = imresize(tempE3,[x,y]);
                E(i,j,:,:,4) = imresize(tempE4,[x,y]);
                        
            end
        end
        
        E = double(E);
        
        %% Run through recursive algorithm to get preliminary values
        
        [a, b] = a_b_recursive(E);
        
        %% Apply modification
        
        I = double(ones(M,N,cropxpixels,cropypixels)) - 1;
        
        for j = 1:N
            for i = 1:M
                I(i,j,:,:) = cropTiles(i,j,:,:) .*a(i,j);
            end
        end
        
        I = I + 1000 + max(b(:));

        for j = 1:N
            for i = 1:M
                I(i,j,:,:) = I(i,j,:,:) - b(i,j);
            end
        end        
        
        I = I + 1 - min(I(:));
        
        mkdir(output);
        for j = 1:N
            for i = 1:M
                hat(:,:) = I(i,j,:,:);
                imageNumber = (j-1)*N+i;
                %imwrite(uint16(hat),[output filename ' Tile - ' int2str(t) '.tif']);
                imwrite(uint16(hat),[output filename num2str(imageNumber,['%0' num2str(numChars) 'd']) '_w' num2str(m) 'Confocal ' num2str(chan(m)) '_MIP.TIF']);
                t = t + 1;
            end
        end
    end
end

warning('on','all');

timeToFinish = toc;

end