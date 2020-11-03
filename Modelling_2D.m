clear all;
clear;

%Set parameters
Vis = 1e16; 
Rho = 2800; 
dPdY = 31000;
Gy = 10; 

%Set parameters for the x and y grid
xsize = 100;
numx = 50;
dx = xsize/(numx-1);
ysize = 100;
numy = 50;
dy = ysize/(numy-1);
k = numx*numy;

%Settin the grid plus L and R
x = 0:dx:xsize;
y = 0:dy:ysize;
L = zeros(k, k);
R = zeros(k, 1);

%Filling in the grid L and R. K bound the outcomes to one vector. 
for i= 1:1:numy
    for j = 1:1:numx
        k = (i-1)*numy+j;
        
        if (j == 1 || j == numx || i == 1 || i == numy )
            L(k,k) = 1;
            R(k,1) = 0;
            
        else  
          
            L(k,k-numy) = 1/(dx^2);                
            L(k,k-1) = 1/(dy^2);
            L(k,k) = -2/(dx^2)-2/(dy^2);
            L(k,k+1) = 1/(dy^2);
            L(k,k+numy) = 1/(dx^2);
            R(k,1) = (1/Vis)*(dPdY-Rho*Gy); 
            
        end
        
    end
  
 end
   
%Computing S
S = L\R;

%Translating S back to a numx*numy matrix
for i = 1:1:numy
    for j = 1:1:numx
        k = (j-1)*numy+i;
        v_y(i,j) = -S(k);
        
    end
end

%Plotting the figure
figure(2)
pcolor(x,y,v_y);
title('2D magma flow pattern for a dike')
xlabel('m')
ylabel('m')
