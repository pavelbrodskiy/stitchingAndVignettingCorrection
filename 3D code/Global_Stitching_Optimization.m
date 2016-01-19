function [timeToFinish] = Global_Stitching_Optimization(M, N, x, y, z, OverlapPercent, CropPercent, input, lb1, ub1, x01, MFE1, MI1, ...
    TolCon1, TolFun1, TolX1, time1, ST11, Trial1, name, chan, numembs, start, filename, output, numChars, gauss)

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
        disp('Reading Images')
        toc
        disp(['Embryo: ' int2str(k) ' Channel: ' int2str(m)]);
        for j = 1:N
            for i = 1:M
                
                fname = [input filename num2str(imageNumber,['%0' num2str(numChars) 'd']) '_w' num2str(m) 'Confocal ' num2str(chan(m)) '_MIP.TIF'];
                info = imfinfo(fname);
                numberOfImages = length(info);
                for kk = numberOfImages:-1:1
                    currentImage = imread(fname, kk, 'Info', info);
                    tempImage(:,:,kk) = currentImage;
                end
            
                if gauss > 0
                    tempImage = double(tempImage);
                    backgroundImage = imgaussfilt3(tempImage, gauss, 'FilterSize', 5);
                    newImage = tempImage - backgroundImage;
                    tempImage = newImage - min(newImage(:)) + 1;
                end
                
                rawTiles(i,j,:,:,:) = tempImage;
                
                imageNumber = imageNumber + 1;
                
            end
        end
        
        %% Crop Images
        
        [~, ~, totalxpixels, totalypixels, totalzpixels] = size(rawTiles);
        
        xCrop = round(totalxpixels * CropPercent);
        yCrop = round(totalypixels * CropPercent);
        xOverlapPixels = round(OverlapPercent * totalxpixels - xCrop * 2);
        yOverlapPixels = round(OverlapPercent * totalypixels - yCrop * 2);
        cropxpixels = round(totalxpixels * (1-CropPercent*2));
        cropypixels = round(totalypixels * (1-CropPercent*2));
        cropTiles = rawTiles(:,:,(xCrop+1):(totalxpixels - xCrop),(yCrop+1):(totalypixels-yCrop),:);
        
        %% Obtain slice average overlap intensity (P)
        disp('Cutting edges')
        toc
        E = zeros(M,N,x,y,z,4);
        
        for i = 1:M
            for j = 1:N
                
                clear tempE1 tempE2 tempE3 tempE4;
                
                tempE1(:,:,:) = cropTiles(i,j,1:xOverlapPixels,:,:);
                tempE2(:,:,:) = cropTiles(i,j,:,(cropypixels - yOverlapPixels):cropypixels,:);
                tempE3(:,:,:) = cropTiles(i,j,(cropxpixels - xOverlapPixels):cropxpixels,:,:);
                tempE4(:,:,:) = cropTiles(i,j,:,1:yOverlapPixels,:);
                
                tempE2 = permute(tempE2,[2,1,3]);
                tempE4 = permute(tempE4,[2,1,3]);                
                
                [xxx, yyy, zzz] = size(tempE1);
                
                for mm = xxx:-1:1
                    clear tempE111 tempE211 tempE311 tempE411
                    tempE111(:,:) = tempE1(mm,:,:);
                    tempE211(:,:) = tempE2(mm,:,:);
                    tempE311(:,:) = tempE3(mm,:,:);
                    tempE411(:,:) = tempE4(mm,:,:);
                    
                    tempE11(mm,:,:) = imresize(tempE111,[x,z]);
                    tempE21(mm,:,:) = imresize(tempE211,[x,z]);
                    tempE31(mm,:,:) = imresize(tempE311,[x,z]);
                    tempE41(mm,:,:) = imresize(tempE411,[x,z]);
                end
                
                for mm = z:-1:1
                    clear tempE111 tempE211 tempE311 tempE411
                    tempE111(:,:) = tempE11(:,:,mm);
                    tempE211(:,:) = tempE21(:,:,mm);
                    tempE311(:,:) = tempE31(:,:,mm);
                    tempE411(:,:) = tempE41(:,:,mm);
                    
                    E(i,j,:,:,mm,1) = imresize(tempE111,[x,y]);
                    E(i,j,:,:,mm,2) = imresize(tempE211,[x,y]);
                    E(i,j,:,:,mm,3) = imresize(tempE311,[x,y]);
                    E(i,j,:,:,mm,4) = imresize(tempE411,[x,y]);
                end
                
                
            end
        end
        
        E = double(E);
        
        %% Run through recursive algorithm to get preliminary values
        
        disp('Starting Recursion');
        toc
        [a, b] = a_b_recursive(E);
        
        %% Apply modification
        
        I = double(ones(M,N,cropxpixels,cropypixels,totalzpixels)) - 1;
        
        for j = 1:N
            for i = 1:M
                I(i,j,:,:,:) = cropTiles(i,j,:,:,:) .*a(i,j);
            end
        end
        
        I = I + 1000 + max(b(:));

        for j = 1:N
            for i = 1:M
                I(i,j,:,:,:) = I(i,j,:,:,:) - b(i,j);
            end
        end        
        
        I = I + 1 - min(I(:));
        
        mkdir(output);
        for j = 1:N
            for i = 1:M
                hat(:,:,:) = I(i,j,:,:,:);
                imageNumber = (j-1)*N+i;
                %imwrite(uint16(hat),[output filename ' Tile - ' int2str(t) '.tif']);
                outputFileName = [output filename num2str(imageNumber,['%0' num2str(numChars) 'd']) '_w' num2str(m) 'Confocal ' num2str(chan(m)) '_MIP.TIF'];
                for K=1:length(hat(1, 1, :))
                    imwrite(hat(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
                K
                end
                
                t = t + 1;
            end
        end
    end
end

warning('on','all');

timeToFinish = toc;

end