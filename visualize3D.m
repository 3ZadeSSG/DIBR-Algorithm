function visualize3D(Cam_XYZ,Distance_Filter)
%{
    ---------------------------------------------
    HELPER FUNCTION : for 3D points visualization
    ---------------------------------------------
    The function takes points projected in XYZ according to depth
    information and generates a 3D scatter plot. For faster rendering
    keep Distance_Filter value close to 250, which will discard any points
    having Z greated than 250, thus eliminating points which are far, the 
    plot process becomes much faster
%}

% Select and remove unecessary points which are far
toDelete = Cam_XYZ(3,:)>Distance_Filter;
Cam_XYZ(:,toDelete) = [];

% Plot with jet color so blue points will be closest
zscaled = Cam_XYZ(3,:);
cn = ceil(max(zscaled));
cm = colormap(jet(cn));
scatter3(Cam_XYZ(1,:),Cam_XYZ(2,:),Cam_XYZ(3,:),[], cm(ceil(zscaled),:),'.');