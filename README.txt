Running the Scripts:-

Deliverable 1 >> main.m  : Running the script will call implemented DIBR function and output virtual image, the image will be plotted and also saved in “Output Results\Output Task 1 - Generate Virtual Image” directory, or output location can be changed. Thus completing first deliverable task.

Deliverable 1 >> DIBR_Multiple.m : Running the script after setting value of N will result generation of views between original and virtual camera with equal distance apart. 
Note since original DIBR function also saves the binary image, in case user wants to test the third deliverable task, it’s better to change the path in DIBR.m function for binary image to some temporary location when running this function because the third task will try to read binary image associated to generated virtual image (m_v)

OR 

Run this second deliverable task script first, and then run first task script so the binary image gets overwritten by m_v’s rather than images generated between m_o and m_v
Output of this function can be seen inside directory “Output Results\Output Task 2 - Generate N Virtual Images”. Thus completing second deliverable task.

Deliverable 1 >> removeArtifacts.m : Running the script will read the virtual image file and binary (mask) image file for trying to reduce the artifacts. It will read those images from:-

Virtual Image: “Output Results\Output Task 1 - Generate Virtual Image\”
Mask Image: “Output Results\Output Task 3 - Remove Cracks, Missing Regions\”

The function will plot the results for comparison. Thus completing third deliverable task.





Other Description related to design choices:-

1. To synthesize the virtual view, the depth information and original camera matrix is used to backproject 2D image points to 3D world coordinates (A visualization function is provided for optional visualization). The virtual camera matrix and depth is used to project and convert it to virtual 2D image.

2. Since task is more related to testing functionality, there is no edge cases checking and no error handling has been implemented.

3. For DIBR, the basic inverse projection can be done in matrix operation, but to show steps the iterative approach has been used.

4. All 3 tasks can be binded in one executable script, but to show individual functionality as per task lists, separate scripts have been created

5. For third task, employing any method is allowed, hence impainting algorithm has been used. Additonally generated image contains lots of smaller holds, which are like noise, hence a median filter has been applied before impainting to get rid of the additional small holes. 
Used MATLAB 2018a version does not has any build in library for TELA Impainting algorithm (Telea, A., 2004. An image inpainting technique based on the fast marching method. Journal of graphics tools, 9(1), pp.23-34.), hence a code has been added to execute that task, comparing it with it’s python implementation in openCV, the execution speed is slow since optimization is not done for making it execute faster.

 
