%{
    ------------------------------------
    Objective 1 : Generate 2D Virtual Image with DIBR
    ------------------------------------
    Description : This script uses C_original and C_virtual camera positions
                 to generate N numebr of virutal images based on virtual
                 camera positioned at N points connecting C_original and
                 C_virtual
%}

%{
    file_depth = Location for image file containing depth info (H,W,3)
    file_texture = Location for image file containing color info (H,W,3)   
    save_output = set it to "false" in case saving generated image is not
    required
%}
file_depth = 'D_original.png';
file_texture = 'V_original.png';
save_output = true;

%{
    ---------------------------
    Original Camera Parameters
    ---------------------------
    K_original  [3 x 3] = Intrinsic Parameters
    Rt_original [3 x 4] = Extrinsic Parameters
    Znear = Nearest Point from Original Camera
    Zfar = Farthese Point from Original Camera
%}
K_original = [1732.87   0.0     943.23;
              0.0       1729.90 548.845040;
              0         0       1];
Rt_original = [1.0 0.0 0.0 0.0;
               0.0 1.0 0.0 0.0;
               0.0 0.0 1.0 0.0];
Zfar = 2760.510889;
Znear = 34.506386;

%{
    ---------------------------
    Virtual Camera Parameters
    ---------------------------
    K_virtual  [3 x 3] = Intrinsic Parameters
    Rt_virtual [3 x 4] = Extrinsic Parameters
%}
K_virtual = [1732.87    0.0     943.23;
             0.0        1729.90 548.845040;
             0          0       1];
Rt_virtual = [1.0 0.0 0.0 1.5924;
              0.0 1.0 0.0 0.0;
              0.0 0.0 1.0 0.0];

image_original = imread(file_texture); % read RGB image
image_depth = imread(file_depth);      % read Depth image
Z_Map = getDepthMap(image_depth,Znear, Zfar); % create actual depth map for DIBR function

image_output = DIBR(image_original,Z_Map,K_original, Rt_original,K_virtual, Rt_virtual); % Invoke DIBR

% If saving output is required then save the image
if save_output
    imwrite(image_output,'Output Results\Output Task 1 - Generate Virtual Image\Output_Virtual_Image.png');
end

%{
    Plot Input vs Output Image
%}
subplot(1,2,1),imshow(image_original),title('Input Image');
subplot(1,2,2),imshow(image_output),title('Output Image');

