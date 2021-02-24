function Z_Map = getDepthMap(Depth_Image,Znear, Zfar)
%{
    -------------------------------------------
    HELPER FUNCTION : For creating depth map
    -------------------------------------------
    Given depth image (H,W,3) where each channel pixel has equal value,
    Znear, and Zfar vlalue the function generates the depth map
%}


Depth_Image = rgb2gray(Depth_Image);
[H,W] = size(Depth_Image);              
Z_Map = zeros(H,W,'double');            
Depth_Image = cast(Depth_Image,'double');

%{
    ------------------------------------------------
     For this task inverse mapping is used for depth
    ------------------------------------------------
                           1
    Z = _____________________________________________
        ((I(ij)/255)x (1/Znear - 1/Zfar)) + (1/Zfar)
%}
for i = 1:H
    for j = 1:W
        t1 = Depth_Image(i,j)/255;
        t2 = (1/Znear) - (1/Zfar);
        t3 = t1*t2;
        t4 = t3+(1/Zfar);
        Z_Map(i,j) = 1/t4;
    end
end