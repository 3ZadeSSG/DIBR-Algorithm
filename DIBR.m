function V_v = DIBR(V_o,D_o,K_o,Rt_o,K_v,Rt_v)

%{
    ------------------------------------------------------------
    HELPER FUNCTION : To be used by scripts for objective 1 & 2
    ------------------------------------------------------------
%}

[H,W,RGB] = size(V_o);        % H and W information about image is required, color channel is not necessary at this point
V_v = zeros(H,W,3);         % Result image that needs to be returned

fx = K_o(1,1);
fy = K_o(2,2);
u0 = K_o(1,3);
v0 = K_o(2,3);

P_o = K_o*Rt_o;
P_v = K_v*Rt_v;
%C_o = -inv(P_o(1:3,1:3))*P_o(:,4);
%C_v = -inv(P_v(1:3,1:3))*P_v(:,4); 

Cam_XYZ = zeros(H*W,3);
index = 1;
for v = 1:H
    for u = 1:W
        x = (u-u0)*D_o(v,u)/fx;
        y = (v-v0)*D_o(v,u)/fy;
        z = D_o(v,u);
        Cam_XYZ(index,:) = [x y z];
        index = index+1;
    end
end
Cam_XYZ = Cam_XYZ';

%{
---------------------------------------------------------------------------
        Optional Info: Visualization of points with Depth information
---------------------------------------------------------------------------
Set target distance 250,so this way additonal points with larger distance 
will be discarded thus making matlab's scatter3 function little faster and
easier to navigate
%}
%visualize3D(Cam_XYZ,250);

M = [Cam_XYZ;ones(1,H*W)];
m_V = P_v*M;

v_points = zeros(H*W,2);
index = 1;
for v = 1:H
    for u = 1:W
        x = m_V(1,index);
        y = m_V(2,index);
        u_P = x/D_o(v,u);
        v_P = y/D_o(v,u);
        v_points(index,:) = [v_P u_P];
        index = index+1;
    end
end
v_points = cast(v_points,'int64');

index = 1;
myMap = zeros(H,W); % Keeping track of already seen points, a better way would have been use of HashMap for saving space
output_image = zeros(H,W,3);
for i = 1:H
    for j = 1:W
        x_o = v_points(index,1);
        y_o = v_points(index,2);
        if(x_o>0 && x_o<=H && y_o>0 && y_o<=W)
            if (myMap(x_o,y_o)==0)
                output_image(i,j,1) = V_o(x_o,y_o,1);
                output_image(i,j,2) = V_o(x_o,y_o,2);
                output_image(i,j,3) = V_o(x_o,y_o,3);
                myMap(x_o,y_o) = 1;
            else
                output_image(i,j,1) = 255;
                output_image(i,j,2) = 255;
                output_image(i,j,3) = 255;
            end
        end
        index = index+1;
    end
end
output_image = cast(output_image,'uint8');
V_v = output_image;


