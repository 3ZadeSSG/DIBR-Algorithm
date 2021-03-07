%{
    ------------------------------------
    Objective 3 : Improve synthesized view
    ------------------------------------
    Description: This script read the generated virutal image, and tries to
                 remove the holes

    Step 1. Normal median filter has been used to remove tiny holes which are
    like noise.
    Step 2. Uses approacch from "An Image Inpainting Technique Based on the Fast Marching
    Method" for impainting. The Python implementation in Open CV library has been tested to
    see the results, this MATLAB implementation is similar implementation
    of Telea Algorithm since MATLAB 2018a does not have impainting functions built in
%}

input_image_path = 'Output Results\Output Task 1 - Generate Virtual Image\Output_Virtual_Image.png';
input_mask_path = 'Output Results\Output Task 3 - Remove Cracks, Missing Regions\Virtual_Bin_Image.png';
step_1_output_image_path = 'Output Results\Output Task 3 - Remove Cracks, Missing Regions\1_Cleaned_After_Masking_Bin_Image.png';
step_1_output_mask_path = 'Output Results\Output Task 3 - Remove Cracks, Missing Regions\1_Cleaned_After_Masking.png';
step_2_output_image_path = 'Output Results\Output Task 3 - Remove Cracks, Missing Regions\2_Cleaned_After_Impainting.png';

generated_image = imread(input_image_path);
generated_binary_image = imread(input_mask_path);

%{
----------------------------------------------------------------
  Step 1. Median Filter when applied on image with Kernel/Mask of 3x3 it can
  remove some but not all holes which are result of  disocclusion
----------------------------------------------------------------
%}
% padding by 'same' since 3x3 kernel will lead to indexing error
padImage=uint16(padarray(generated_image,[1,1],'replicate','both'));
newImage=padImage;
[height,width,RGB]=size(padImage);
for i=2:height-1
    for j=2:width-1
        for k=1:3
            mask=padImage(i-1:i+1,j-1:j+1,k);
            newImage(i,j,k)=median(median(mask));
        end
    end
end
result=uint8(newImage(2:height-1,2:width-1,:)); % Since orignal image was padded, select the appropriate size image
height = height-2; width = width-2;
mask_img = generated_binary_image;
for i=1:height
    for j=1:width
        if(mask_img(i,j)==255) && (mean(result(i,j,:)~=255))
            mask_img(i,j) = 0;
        end
    end
end

imwrite(mask_img,step_1_output_image_path);
imwrite(result,step_1_output_mask_path);

% Compare result of before and after processing
%subplot(1,2,1),imshow(generated_image),title("Before Applying Masking");
%subplot(1,2,2),imshow(result),title("After Applying Masking");

disp('Starting Step 2');
%{
----------------------------------------------------------------
   Step 2. Impainting approach as per the algorithm (Telea 2004)
----------------------------------------------------------------
%}
band = []; %empty priority queue
img = result;
mask = zeros(size(mask_img));
mask(mask_img>0)=1;
output_image = inpaint(img,mask,5); % call impainting with radius 5

imwrite(output_image,step_2_output_image_path); % save final image

% Plot results for comparision
subplot(2,2,1),imshow(generated_image),title('Input Image');
subplot(2,2,2),imshow(result),title('After Median Filter');
subplot(2,2,3),imshow(output_image),title('Final Image');


%{
----------------------------------------------------------------
    Compute X Y gradient of pixel
----------------------------------------------------------------
%}
function [grad_y,grad_x] = pixel_gradient(y, x, height, width, vals, flags)
    KNOWN = 0;BAND = 1; UNKNOWN = 2;INF = 1e6;EPS = 1e-6;
    val = vals(y,x);
    prev_y = y - 1;
    next_y = y + 1;
    if prev_y < 1 || next_y >= height
        grad_y = INF;
    else
        flag_prev_y = flags(prev_y, x);
        flag_next_y = flags(next_y, x);

        if flag_prev_y ~= UNKNOWN && flag_next_y ~= UNKNOWN
            grad_y = (vals(next_y, x) - vals(prev_y, x)) / 2.0;
        elseif flag_prev_y ~= UNKNOWN
            grad_y = val - vals(prev_y, x);
        elseif flag_next_y ~= UNKNOWN
            grad_y = vals(next_y, x) - val;
        else
            grad_y = 0.0;
        end
    end

    prev_x = x - 1;
    next_x = x + 1;
    if prev_x < 1 || next_x >= width
        grad_x = INF;
    else
        flag_prev_x = flags(y, prev_x);
        flag_next_x = flags(y, next_x);

        if flag_prev_x ~= UNKNOWN && flag_next_x ~= UNKNOWN
            grad_x = (vals(y, next_x) - vals(y, prev_x)) / 2.0;
        elseif flag_prev_x ~= UNKNOWN
            grad_x = val - vals(y, prev_x);
        elseif flag_next_x ~= UNKNOWN
            grad_x = vals(y, next_x) - val;
        else
            grad_x = 0.0;
        end
    end
end


%{
----------------------------------------------------------------
    Based on neighbors return RGB value of pixel to be impainted
----------------------------------------------------------------
%}
function pixel_sum = inpaint_pixel(y, x, img, height, width, dists, flags, radius)
    KNOWN = 0;BAND = 1; UNKNOWN = 2;INF = 1e6;EPS = 1e-6;
    dist = dists(y,x);
    [dist_grad_y,dist_grad_x] = pixel_gradient(y, x, height, width, dists, flags);
    pixel_sum = zeros(3,'double');
    weight_sum = 0.0;
    for nb_y=y-radius:y+radius
        if nb_y < 1 || nb_y >= height
            continue;
        end
        for nb_x=x-radius:x+radius
            if nb_x < 1 || nb_x >= width
                continue;
            end
            if flags(nb_y, nb_x) == UNKNOWN
                continue;
            end
            dir_y = y - nb_y;
            dir_x = x - nb_x;
            dir_length_square = dir_y^2 + dir_x^2;
            dir_length = sqrt(dir_length_square);
            if dir_length > radius
                continue;
            end
            dir_factor = abs(dir_y * dist_grad_y + dir_x * dist_grad_x);
            if dir_factor == 0.0
                dir_factor = EPS;
            end
            nb_dist = dists(nb_y, nb_x);
            level_factor = 1.0 / (1.0 + abs(nb_dist - dist));
            dist_factor = 1.0 / (dir_length * dir_length_square);
            weight = abs(dir_factor * dist_factor * level_factor);

            pixel_sum(1) = pixel_sum(1) + weight * img(nb_y, nb_x, 1);
            pixel_sum(2) = pixel_sum(2) + weight * img(nb_y, nb_x, 2);
            pixel_sum(3) = pixel_sum(3) + weight * img(nb_y, nb_x, 3);
            weight_sum = weight_sum+weight;
        end
    end
    pixel_sum = pixel_sum./weight_sum;
end




%{
----------------------------------------------------------------
   eikonal equation as per expression (3) in original paper
----------------------------------------------------------------
%}
function dist = solve_eikonal(y1, x1, y2, x2, height, width, dists, flags)
    KNOWN = 0;BAND = 1; UNKNOWN = 2;INF = 1e6;EPS = 1e-6;
    dist = INF;
    if y1 < 1 || y1 >= height || x1 < 1 || x1 >= width
        dist = INF;
        return;
    end
    if y2 < 1 || y2 >= height || x2 < 1 || x2 >= width
        dist = INF;
        return;
    end

    flag1 = flags(y1, x1);
    flag2 = flags(y2, x2);

    if flag1 == KNOWN && flag2 == KNOWN
        dist1 = dists(y1, x1);
        dist2 = dists(y2, x2);
        d = 2.0 - (dist1 - dist2)^2;
        if d>0.0
            r = sqrt(d);
            s = (dist1 + dist2 - r) / 2.0;
            if s >= dist1 && s >= dist2
                dist = s;
                return;
            end
            s = s+r;
            if s >= dist1 && s >= dist2
                dist = s;
                return;
            end
            dist = INF;
            return;
        end
    end
    if flag1 == KNOWN
        dist1 = dists(y1, x1);
        dist = 1.0 + dist1;
        return;
    end
    if flag2 == KNOWN
        dist2 = dists(y2, x2);
        dist = 1.0 + dist2;
        return;
    end
    dist = INF;
    return;
end


%{
----------------------------------------------------------------
returns priority queue with distance value between mask and pixels outside
mask with Fast Marching Method
----------------------------------------------------------------
%}
function [band] = compute_outside_dists(height,width,dists,flags,radius,band)
    KNOWN = 0;BAND = 1; UNKNOWN = 2;INF = 1e6;EPS = 1e-6;
    orig_flags = flags;
    [flag_H,flag_W] = size(flags);

    for i=1:flag_H
        for j=1:flag_W
            if orig_flags(i,j) == KNOWN
                flags(i,j) = UNKNOWN;
            end
            if orig_flags(i,j) == UNKNOWN
                flags(i,j) = KNOWN;
            end
        end
    end
    last_dist = 0.0;
    double_radius = radius*2;


    [band_rows,~] = size(band);
    disp(band_rows);
    while band_rows>0
        if last_dist >= double_radius
            break
        end
        [band,top] = dequeue(band);
        y = top(2);
        x = top(3);
        flags(y,x) = KNOWN;

        neighbors = [y-1 x;
                     y x-1;
                     y+1  x;
                     y x+1];
        [n_rows,~] = size(neighbors);
        for n_idx = 1:n_rows
            nb_y = neighbors(n_idx,1);
            nb_x = neighbors(n_idx,2);
            if nb_y < 1 || nb_y >= height || nb_x < 1 || nb_x >= width
                continue;
            end
            if flags(nb_y, nb_x) ~= UNKNOWN
                continue;
            end
            last_dist = min([
                solve_eikonal(nb_y-1, nb_x, nb_y, nb_x-1, height, width, dists, flags)
                solve_eikonal(nb_y+1, nb_x, nb_y, nb_x+1, height, width, dists, flags)
                solve_eikonal(nb_y-1, nb_x, nb_y, nb_x+1, height, width, dists, flags)
                solve_eikonal(nb_y+1, nb_x, nb_y, nb_x-1, height, width, dists, flags)
            ]);
            dists(nb_y,nb_x) = last_dist;
            flags(nb_y,nb_x) = BAND;
            band = enqueue(band,[last_dist,nb_y,nb_x]);
        end
        [band_rows,~] = size(band);
    end
    dists = dists.*(-1.0);
    disp('Compute Outside Dist Last Line');
end

%{
----------------------------------------------------------------
Returns initial mask, flags, narrow band priority queue
----------------------------------------------------------------
%}
function [dists, flags, band] = init(height, width, mask, radius)
    KNOWN = 0;BAND = 1;UNKNOWN = 2;INF = 1e6;EPS = 1e-6;
    disp('Init was called');
    dists = zeros(height,width,'double')+INF;
    flags = mask.*UNKNOWN;
    band = [];
    [mask_H,mask_W] = size(mask);
    [non_zero_size,~] = size(nonzeros(mask));
    index = 0;
    for i=1:mask_H
        for j=1:mask_W
            if mask(i,j) ~= 0
                y = i;
                x = j;
                neighbors = [y-1 x;
                             y x-1;
                             y+1 x;
                             y x+1];

                [n_rows,~] = size(neighbors);
                for n_idx = 1:n_rows
                    nb_y = neighbors(n_idx,1);
                    nb_x = neighbors(n_idx,2);
                    if nb_y < 1 || nb_y >= height || nb_x < 1 || nb_x >= width
                        continue;
                    end
                    if flags(nb_y, nb_x) == BAND
                        continue;
                    end
                    if mask(nb_y,nb_x) == 0
                        flags(nb_y,nb_x) = BAND;
                        dists(nb_y,nb_x) = 0.0;
                        band = enqueue(band,[0.0 nb_y nb_x]);
                    end
                end
                index = index+1;
            end
        end
    end
    disp(size(band));
    disp('Calling compute outside dists');
    [band] = compute_outside_dists(height,width,dists,flags,radius,band);
end



%{
----------------------------------------------------------------
for each region, iteratively impaints it with closest neighbor
----------------------------------------------------------------
%}
function img = inpaint(img,mask,radius)
    KNOWN = 0;BAND = 1;UNKNOWN = 2;INF = 1e6;EPS = 1e-6;
    [height,width,RGB] = size(img);
    [dists, flags, band] = init(height, width, mask, radius);
    disp(size(dists));
    disp(size(flags));
    disp(size(band));

    iteration = 0;
    [p_rows,~] = size(band);
    while p_rows>0
        iteration = iteration + 1;
        [band,top] = dequeue(band);
        y = top(2);
        x = top(3);
        flags(y,x) = KNOWN;
        neighbors = [y-1 x;
                     y x-1;
                     y+1  x;
                     y x+1];
        [n_rows,~] = size(neighbors);
        for n_idx = 1:n_rows
            nb_y = neighbors(n_idx,1);
            nb_x = neighbors(n_idx,2);
            if nb_y < 1 || nb_y >= height || nb_x < 1 || nb_x >= width
                continue
            end
            if flags(nb_y, nb_x) ~= UNKNOWN
                continue
            end
            nb_dist = min([
                solve_eikonal(nb_y - 1, nb_x, nb_y, nb_x - 1, height, width, dists, flags)
                solve_eikonal(nb_y + 1, nb_x, nb_y, nb_x + 1, height, width, dists, flags)
                solve_eikonal(nb_y - 1, nb_x, nb_y, nb_x + 1, height, width, dists, flags)
                solve_eikonal(nb_y + 1, nb_x, nb_y, nb_x - 1, height, width, dists, flags)
            ]);
            dists(nb_y,nb_x) = nb_dist;
            pixel_vals = inpaint_pixel(nb_y, nb_x, img, height, width, dists, flags, radius);

            img(nb_y, nb_x, 1) = pixel_vals(1);
            img(nb_y, nb_x, 2) = pixel_vals(2);
            img(nb_y, nb_x, 3) = pixel_vals(3);

            flags(nb_y,nb_x) = BAND;
            band = enqueue(band,[nb_dist,nb_y,nb_x]);
        end
        [p_rows,~] = size(band);
    end
    disp(iteration);
end



%{
----------------------------------------------------------------
Returns queue with added new value
----------------------------------------------------------------
%}
function band = enqueue(band,y)
    band = [band;y];
end
%{
----------------------------------------------------------------
Pops and returns top element as per priority queue rule
----------------------------------------------------------------
%}
function [band,top] = dequeue(band)
    [val,idx] = min(band(:,3));
    top = band(idx,:);
    band(idx,:) = [];
end
