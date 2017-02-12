%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written for Course :- Computational Electromagnetics, Fall 2011
%                       Department of Electrical Engineering
%                       Indian Institute of Technology Madras
%                       Chennai - 600036, India
%
% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumanthra Chaudhuri, M.S. Electrical Engg.
%                       Projesh Basu, M.S. Electrical Engg.
%                       Nikhil Kumar CS, M.S. Electrical Engg.
%
% Instructor :- Ananth Krishnan
%               Assistant Professor
%               Department of Electrical Engineering
%               Indian Institute of Technology Madras
%
% Any correspondance regarding this program may be addressed to
% Prof. Ananth Krishnan at 'computational.em.at.iit.madras@gmail.com'
%
% Copyright/Licensing :- For educational and research purposes only. No
% part of this program may be used for any financial benefit of any kind
% without the consent of the instructor of the course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Finite Element Method (FEM) solution to Laplacian"
% 
%         Objective of the program is to solve for the steady state voltage
% distribution in a region 0<x<30 units, 0<y<30 units, given that one of the 
% sides of square is excited with a voltage of 45*(x/xmax)*((xmax-x)/xmax) 
% Volts (xmax=30) and all other sides are maintained at 0 Volts. This voltage 
% at the boundary is symmetrical with its maximum value at centre of the 
% boundary namely at x=15. 
%
%         The gridding for the Finite Element Method is very similar to that 
% of Finite Difference Method with the exception that the squares in the
% FDM are further divided into two triangles using a right diagonal to form
% the triangular elements of the FEM grid. The points are given a bilinear 
% fit in x and y for voltage distribution. The parameters for the fit are 
% obtained using standard Ritz FEM formulation.
%
%         The tolerance in error between iterations is kept at 0.01 V. This 
% may be tweaked to a higher or lower value for lower or higher accuracy 
% respectively. Imagesc command by default uses image axis settings, which
% are different from normal plot command and hence x and y axis may look 
% flipped. Read Matlab documentation on imagesc for more details.
%
%         Program stops when iteration number in plot does not change or can 
% be closed anytime by just closing the plot window. On normal completion,the
% program plots the electric field in a quiver plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Clearing variables in memory and Matlab command screen
clear all;
clc;

%Dimensions of the simulation grid in x (xdim) and y (ydim) directions
xdim=2;
ydim=2;

%Area of square elements in the grid
area_sq_elem=1;

%Step size in the grid in cartesian directions viz. side of the square
%element
step_size=area_sq_elem^0.5;

%No of square elements
sq_elems=(xdim*ydim)/area_sq_elem;

%Dimension in x direction in steps
dim2=xdim/step_size;

%Dimension in y direction in steps
dim1=ydim/step_size;

%No of triangular elements
no_of_elems=sq_elems*2;

%Total no of points for voltage distribution
total_no_pts=(dim1+1)*(dim2+1);

%grid_pts matrix consisting of numbering of all the points in the domain in
%a page read fashion from left to right and then top to bottom along with
%the x and y coordinates of the points in the adjacent columns of the
%matrix
grid_pts=zeros(total_no_pts,3);

%look_up_table matrix consisting of numbering of all the triangular
%elements in a page read fashion just like before along with the cardinal
%numbers of the three points making the triangular elements (obtained from
%the grid_pts matrix) in the adjacent columns of the matrix
look_up_table=zeros(no_of_elems,4);

%Populating the first columns of the above two matrices with serial numbers
i1=1:1:no_of_elems;
look_up_table(:,1)=i1';
i1=1:1:total_no_pts;
grid_pts(:,1)=i1';

% Polpulating the x coordinate and y coordinate columns of the grid_pts
% matrix with the appropriate values from the knowledge of the grid
% geometry
for i1=1:1:total_no_pts
    if rem(i1,dim1+1)==0
        grid_pts(i1,3)=((i1/(dim1+1))-1)*step_size;
    else
        grid_pts(i1,3)=(floor(i1/(dim1+1)))*step_size;
    end
    if rem(i1,dim1+1)==0
        grid_pts(i1,2)=dim1*step_size;
    else
        grid_pts(i1,2)=(i1-(dim1+1)*floor(i1/(dim1+1))-1)*step_size;
    end
end


% Polpulating the second, third and fourth columns of the look_up_table 
% matrix with cardinal numbers (obtained from the grid_pts matrix) of three 
% points forming the corresponding triangular element in the first column.
for i1=1:1:no_of_elems
    k=ceil(i1/(2*dim1))-1;
    if rem(i1,2)==1
        look_up_table(i1,2)=(i1+1)/2+dim1+1+k;
        look_up_table(i1,3)=(i1+1)/2+k;
        look_up_table(i1,4)=(i1+1)/2+1+k;
    else
        look_up_table(i1,2)=i1/2+dim1+1+k;
        look_up_table(i1,3)=i1/2+1+k;
        look_up_table(i1,4)=i1/2+dim1+2+k;
    end
end
        
% Area of all triangular elements are equal to half the area of all square
% elements as they are obtained by dividing the square elements by 2.
area_of_elem=area_sq_elem/2;

% Initializing and populating the P and Q matrices. P and Q matrices are 
% defined for every element as P=(y2-y3,y3-y1,y1-y2) & Q=(x3-x2,x1-x3,x2-x1)
% where (x1,y1),(x2,y2),(x3,y3) are the coordinates of the three points
% making up the triangular element starting from the bottom left point and
% moving anti-clockwise. Also, C matrices are iitialized for local eleements
% which is then accumulated into a global matrix
P=zeros(no_of_elems,3);
Q=zeros(no_of_elems,3);
C=zeros(no_of_elems,3,3);
for i1=1:1:no_of_elems
    P(i1,1)=grid_pts(look_up_table(i1,3),3)-grid_pts(look_up_table(i1,4),3);
    P(i1,2)=grid_pts(look_up_table(i1,4),3)-grid_pts(look_up_table(i1,2),3);
    P(i1,3)=grid_pts(look_up_table(i1,2),3)-grid_pts(look_up_table(i1,3),3);
    Q(i1,1)=grid_pts(look_up_table(i1,4),2)-grid_pts(look_up_table(i1,3),2);
    Q(i1,2)=grid_pts(look_up_table(i1,2),2)-grid_pts(look_up_table(i1,4),2);
    Q(i1,3)=grid_pts(look_up_table(i1,3),2)-grid_pts(look_up_table(i1,2),2);
    for j=1:1:3
        for k=1:1:3
            C(i1,j,k)=(1/(4*area_of_elem))*(P(i1,j)*P(i1,k)+Q(i1,j)*Q(i1,k));
        end
    end
end

% Initializing matrix C_global from which propagation voltage distribution 
% is calculated using iterative analysis
C_global=zeros(total_no_pts,total_no_pts);

% Calculating matrix C1_global by adding up local C matrices
for i=1:1:no_of_elems
    for j=1:1:3
        for k=1:1:3
            C_global(look_up_table(i,j+1),look_up_table(i,k+1))=C_global(look_up_table(i,j+1),look_up_table(i,k+1))+C(i,j,k);
        end
    end
end

%Initializing previous (V_prev) voltage distribution matrix viz the voltage
%distribution matrix before any iteration
V_prev=zeros(dim1+1,dim2+1);

% Giving Dirichlet boundary cnditions for previous voltage distribution matrix
i=0:1:dim1;
xmax=xdim/step_size;
V_prev(1:dim1+1,1)=(1/(xmax^2))*45*i.*(xmax-i);

%Initializing present (V_now) voltage distribution matrix viz the voltage
%distribution matrix after any iteration and the matrix is given a critical
%difference from V_prev matrix at (2,2) to enter loop of update.
V_now=V_prev;
V_now(2,2)=1;

%Iteration counter iter initialized
iter=0;

%Update loop begins with condition as difference between two matrics being
%greater than 0.01
while(max(max(abs(V_prev-V_now)))>1e-2)

    %Iteration counter incremented
    iter=iter+1;
    
    %Initial critical difference between matrices removed in first
    %iteration
    if iter==1
        V_now(2,2)=0;
    end
    
    %Updating V_prev of present iteration with V_now of previous iteration 
    V_prev=V_now;
    %Clearing V_now for update
    V_now(:,:)=0;
    
    %Updating V_now using standard update procedure obtained from FEM
    %analysis
    for k1=1:1:dim2+1
        for k2=1:1:dim1+1
            for i1=1:1:dim2+1
                for i2=1:1:dim1+1
                    if(C_global((k1-1)*(dim1+1)+k2,(i1-1)*(dim1+1)+i2)>0 || C_global((k1-1)*(dim1+1)+k2,(i1-1)*(dim1+1)+i2)<0)
                        if(i1<k1 || i2<k2 || i1>k1 || i2>k2)
                            V_now(k2,k1)=V_now(k2,k1)-(1/C_global((k1-1)*(dim1+1)+k2,(k1-1)*(dim1+1)+k2))*V_prev(i2,i1)*C_global((k1-1)*(dim1+1)+k2,(i1-1)*(dim1+1)+i2);
                        end
                    end
                end
            end
        end
    end
    
    % Giving Dirichlet boundary conditions for present voltage distribution matrix
    i=0:1:dim1;
    xmax=xdim/step_size;
    V_now(1:dim1+1,1)=(1/(xmax^2))*45*i.*(xmax-i);
    V_now(1:dim1+1,dim2+1)=0;
    V_now(1,1:dim2+1)=0;
    V_now(dim1+1,1:dim2+1)=0;
    
    %Plotting the voltage distribution using color scaled plot after every
    %iteration
%     imagesc(V_now);colorbar;
%     title(['\fontsize{20}Voltage distribution on a ',int2str(xdim),' x ',int2str(ydim),' grid at iteration no ',int2str(iter)],'Color','k');
%     xlabel('x-axis','FontSize',20);
%     ylabel('y-axis','FontSize',20);
%     set(gca,'FontSize',20);
%     getframe;
end
%Update loop ends


% Electric field distribution obtained from voltage distribution using the
% gradient function
x=1:1:dim2+1;
y=1:1:dim1+1;
[Ex,Ey]=gradient(V_now);

%Plotting of electric field with directed arrows using Quiver plot command
%quiver(x,y,-Ex,-Ey);
% title(['\fontsize{20}Electric field lines distribution on the ',int2str(xdim),' x ',int2str(ydim),' grid'],'Color','k');
% xlabel('x-axis','FontSize',20);
% ylabel('y-axis','FontSize',20);
% set(gca,'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%