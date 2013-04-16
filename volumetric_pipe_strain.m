% Function to calculate volumetric surface strain from a
% prolate source at depth

%--------------------------------------------------

%variables:

%c1 = top depth of pipe
%c2 = bottom of the pipe
%r = sqrt(x^2 + y^2)
%R1 = sqrt(r^2 + c1^2)
%R2 = sqrt(r^2 + c2^2) 
% a = source radius [m]
% P = source pressure difference from surround [Pa]
% nu = Poisson Ratio
% G = Shear Modulus (Rigidity) [Pa]
% E = Young's Modulus [Pa]
% x = radial distance on the surface [m]
% y = ?
% ur = radial displacement

%----------------------------------------------------

clear all
close all
%Inputs

c1 = 1000;
c2 = 50000;
x = -4000:1:4000;
y = -4000:1:4000;
%r = -50000:50:50000;
r = sqrt(x.^2 + y.^2);
R1 = sqrt(r.^2 + c1.^2);
R2 = sqrt(r.^2 + c2.^2); 
a = 500;
P = 10e6;
nu = 0.25;
G = 8e9;


%----------------------------------------------------
%calculations for ur

ur = (a^2*P)/(4*G);
j = (c1.^3)./R1.^3; 
k = (2*c1*(-3+5*nu))./R1; 
l = (5*c2.^3*(1-2*nu)-2*c2*r.^2*(-3+5*nu))./R2.^3; 
ur = ur*(j+k+l);
ur = ur.*(x./r.^2);


%------------------------------------------------------
%calculations for tangential pipe strain

%E00 = ur/x
Normal= (a^2*P)./(4*G*c1.^2); 

tang_strain=ur./r;

norm_tang_strain=tang_strain./Normal;

%--------------------------------------------------------
%calculations for radial pipe strain

%Err = diff ur/ diff x
a = diff(ur);
b = diff(x);

rad_strain=a./b;

norm_rad_strain=rad_strain./Normal;

%--------------------------------------------------------
%calculations for volumetric pipe strain

%EV = (1-2*nu/1-nu)*(Err+E00)
volumetric_pipe=(1-2*nu./1-nu)*(rad_strain + tang_strain(:,1:max(size(tang_strain))-1));

norm_volumetric_pipe=((1-(2*nu))./(1-nu))*(norm_rad_strain + abs(norm_tang_strain(1:8000)));
%norm_volumetric_pipe=volumetric_pipe./Normal;

top_depth = x./c1;

%----------------------------------------------------------------
%Plot Graph

plot (top_depth(1:8000),norm_volumetric_pipe, 'LineWidth', 2)


%----------------------------------------------------------------------------------
% Tidy Graph
%----------------------------------------------------------------------------------

% Set Graph Title in fontsize
title('Volumetric Pipe Strain (Lisowski)', 'FontSize', 12, 'FontName', 'Arial');

% Set Axis
xlabel('Distance (Top Depths)', 'FontSize', 12, 'FontName', 'Arial')
ylabel('Volumetric Strain (strain units)', 'FontSize', 12)

% Set Graph Background Color
set(gca,'Color',[1 1 1]);

% Show Grid Lines
%grid minor
%grid
grid off

% Set Graph Limits
%xMin xMax yMin yMax
%axis([-1.e4 1.e4 -2e-7 10.5e-6])
waitforbuttonpress()

%legend command
%line_1_name = 'Analytical Model';
%line_2_name = 'Numerical Model';
%legend(line_1_name, line_2_name, 'Location','NorthEast')

%----------------------------------------------------------------------------------
% Save high resolution version of graph to working directory
%----------------------------------------------------------------------------------

% Save figure 1 to jpeg,called output_test.jpg
% at a resolution of 500 dots per inch
% text is (for commercial print) 300
% images are 2000
print(1,'-djpeg','lisowski_pipe_volumetric','-r500')

%----------------------------------------------------------------------------------

