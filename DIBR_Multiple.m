%{
    ------------------------------------
    Objective 2 : Generate N numebr of Virtual Images
    ------------------------------------
    Description : This script uses C_original and C_virtual camera positions
                 to generate N numebr of virutal images based on virtual
                 camera positioned at N points connecting C_original and
                 C_virtual

    Drawbacks: Since as per task description, DIBR function needs to be
               used, the script will be slow because each time it calls
               the function where common step like 2D to 3D reprojection 
               is done.
%}

N = 5;
file_depth = 'D_original.png';
file_texture = 'V_original.png';

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

%{
    Calculate the C1 and C2, then generate [R|t] i.e., t element for each
    point i=1:N between C1 and C2 at equal distance.
%}
P_o = K_original*Rt_original;       % Original Camera Projection Matrix
P_v = K_virtual*Rt_virtual;         % Virtual Camera Projection Matrix

C_o = -inv(P_o(1:3,1:3))*P_o(:,4);  % Projection center C for original camera
C_v = -inv(P_v(1:3,1:3))*P_v(:,4);  % Projection center C for virtual camera

% Calculate distance and then distnace ratio t
d = sqrt(((C_v(1,1)-C_o(1,1))^2)+((C_v(2,1)-C_o(2,1))^2)+((C_v(3,1)-C_o(3,1))^2));
dt = d/(N+2);
t = dt/d;

% for each distance ration step, go from original to virtual points
x_o = C_o(1,1); y_o = C_o(2,1); z_o = C_o(3,1);
x_v = C_v(1,1); y_v = C_v(2,1); z_v = C_v(3,1);

% calculate Rt_virtual for each i=1 to N number of points and generate image 
for i=1:N
   x_t = ((1-t)*x_o)+(t*x_v);
   y_t = ((1-t)*y_o)+(t*y_v);
   z_t = ((1-t)*z_o)+(t*z_v);
   
   x_o = x_t; y_o = y_t; z_o = z_t;  % new virutal point will be t step from original point
   Rt_virtual(:,4) = -[x_t;y_t;z_t]; % Since P = K[R|t], and P = KR[I|-C], so this calculated 't' for each Rt matrix
   virtual_image = DIBR(image_original,Z_Map,K_original, Rt_original,K_virtual, Rt_virtual); % call DIBR and save the image for i
   
   file_output = sprintf('%s_%d_%s','Output Results\Output Task 2 - Generate N Virtual Images\N',i,'Image.png');
   imwrite(virtual_image,file_output);
end