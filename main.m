%% USER INPUT

cam_param = 0.7;
% Represents the ratio "how much" the camera blurres the image. This is
% used to control how much the speed of an object impacts the length of the
% PSFs.

noise_var = 0.0080;
% The variance of the noise. Used to estimate the NSR of the iamges.

data_to_display = 1;
% 1: final output
% 2: objects original and deblurred + corresponding PSFs
% 3: segmentation and tracking information visualized

video_name = 's1.AVI';
% The filename of the video to precess.

skip_frames = 25;
% Frames of the beginning of the video to skip.


%% ALGORITHM

iptsetpref('ImshowBorder','tight');

VID = VideoReader(video_name);
for skip = 1:skip_frames
    readFrame(VID);
end

non_black = 0;
centroids_old = 0;
numObj_old = 0;

frame_nr = 1;

while(hasFrame(VID))
    
    % read image
    originalImage = readFrame(VID);

    % grayscale image
    grayImage = rgb2gray(originalImage);

    % threshold image
    thresh = 0.1; % (0-255)
    binaryImage = imbinarize(grayImage, thresh);

    % pre-process image
    binaryImage = imfill(binaryImage, 'holes');
    binaryImage = imerode(binaryImage, strel('cube', 5));
    binaryImage = imdilate(binaryImage, strel('cube', 5));
    
    % get centers, bounding boxes
    props = regionprops(binaryImage, 'Centroid', 'BoundingBox');
   
    % count objects
    numObj = size(props, 1); % get size of first element in props

    % make struct with relevent data
    objects = struct('blob', numObj, 'distances', numObj, 'match', numObj, 'vector', numObj, 'PSF', numObj, 'debl', numObj);
    
    % crop out the objects for deblurring
     for i=1 : numObj
        objects(i).blob  = imcrop(originalImage, [floor(props(i).BoundingBox(1)) floor(props(i).BoundingBox(2)) ceil(props(i).BoundingBox(3)) ceil(props(i).BoundingBox(4))] );
     end
     
    % track objects
    if (non_black == 1)
        for i=1 : numObj % calculate distance from each objects on the current frame to each of the old frame
            
            objects(i).distances = zeros(1, numObj_old); % space for distance of each current object to each of the old
            
            for j=1 : numObj_old % get the distances into the strcut as arrays
                objects(i).distances(j) = sqrt( (props(i).Centroid(1) - centroids_old(j).Centroid(1))^2 + (props(i).Centroid(2) - centroids_old(j).Centroid(2))^2);
            end
        
        end
        
        % now search for the minimum and index, the number indicates with what
        % image from the old frame the new one matches.
        for i=1 : numObj
             objects(i).match = find(objects(i).distances == min(objects(i).distances));
        end
        
        % are there new elements? if yes they will generate double matches,
        % sort them out by using the better match
        for i=1 : numObj
            cur_match = objects(i).match;
            best_dist = 10000000;
            best_match = 0;
            for j=1 : numObj
                if (objects(j).match == cur_match)
                   if (min(objects(j).distances) < best_dist) 
                      best_dist = min(objects(j).distances);
                      best_match = j;
                   end
                end
            end
            % now we have the best match numer
            for j=1 : numObj
                if (objects(j).match == cur_match)
                   if (best_match ~= j)
                       objects(j).match = 0;
                   end
                end
            end
        end
       
    end
    
    % get the motion vector of each matched object (x, y)
    % objects with no match will have no vector!
    for i=1 : numObj
        if (objects(i).match > 0 && numObj_old > 0)
            objects(i).vector = zeros (1, 2);
            objects(i).vector(1) = props(i).Centroid(1) - centroids_old(objects(i).match).Centroid(1);
            objects(i).vector(2) = props(i).Centroid(2) - centroids_old(objects(i).match).Centroid(2);
        end
    end
    
    % now we have the bounduing boxes and motion vectors (direction and speed)
    % of each individual object!
    % now calculate the PSFs from that and the user input variable that
    % relates to the internal camera settings like exposure time etc.
    for i=1 : numObj
        if (objects(i).match > 0 && numObj_old > 0)
            b = [1 0];
            objects(i).PSF = fspecial('motion', sqrt( ((objects(i).vector(1) - objects(i).vector(2))^2) * cam_param), -(acos(min(1,max(-1, objects(i).vector(:).' * b(:) / norm(objects(i).vector) / norm(b) ))) * 57.2957795131) );
        end
    end
    clear b; 
    
    % deblur the blobs using wiener filter
    for i=1 : numObj
        if (objects(i).match > 0 && numObj_old > 0)
            sizeBLB = size(objects(i).blob);
            sizePSF = size(objects(i).PSF);
            if (sizeBLB(1) > sizePSF(1) && sizeBLB(2) > sizePSF(2))
                Id = im2double(objects(i).blob);
                SNR = noise_var / var(Id(:)); % noise variance / image variance
                objects(i).debl = deconvwnr(objects(i).blob, objects(i).PSF, SNR);
            end
        end
    end
    clear sizeBLB sizePSF Id;

    % the display
    if (data_to_display == 0)
    % show nothing
    elseif (data_to_display == 1)
        sz = size(binaryImage);
        output = zeros(sz(1), sz(2));
        
        % standatd display
        imshow(output);
        hold on;

        for k = 1 : numObj
            blob_sz = size(objects(k).debl);
            if (blob_sz(1) > 1) && (blob_sz(2) > 1)
               image(props(k).BoundingBox(1), props(k).BoundingBox(2), objects(k).debl);
            end

        end

        hold off;
    
    elseif (data_to_display == 2)
        sz = size(binaryImage);
        output = zeros(sz(1), sz(2));
        
        % show all objects and their PSFs and beblurred versions
        for i=1 : numObj
            p1 = subplot(3, numObj, i);           imshow(objects(i).blob);
            p2 = subplot(3, numObj, i+numObj*2);  imshow(uint8(rescale(objects(i).PSF, 0, 255)));
            p3 = subplot(3, numObj, i+numObj);    imshow(objects(i).debl);
            linkaxes([p3, p2, p1]);
        end
        
    elseif (data_to_display == 3)
        sz = size(binaryImage);
        output = zeros(sz(1), sz(2));
        
        % show the original images with segmentation and tracking information
        imshow(originalImage); 
        hold on;

         for i=1 : numObj
            if (objects(i).match > 0 && numObj_old > 0)
                rectangle('Position', props(i).BoundingBox, 'LineWidth', 1, 'EdgeColor', 'r');
                quiver(props(i).Centroid(1), props(i).Centroid(2), objects(i).vector(1), objects(i).vector(2), 0, 'w', 'MaxHeadSize', 5, 'LineWidth', 2)
            end
         end

         hold off;
            
    end

    % pausing for frame by frame stepping
    disp('Press a key for next frame nr. ')
    disp(num2str(frame_nr));
    pause;
    
    frame_nr = frame_nr +1;
    
    % prepare for next iteration
    if (size(props, 1) ~= 0)
        non_black = 1;
        centroids_old = struct('Centroid', {props.Centroid});
        numObj_old = numObj;
    end

end