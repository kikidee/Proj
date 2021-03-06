% Function to calculate radial surface strain from a
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
%calculations for radial pipe strain


%Err = diff ur/ diff x
Normal= (a^2*P)./(4*G*c1.^2); 

a = diff(ur);
b = diff(x);

rad_strain=a./b;

norm_rad_strain=rad_strain./Normal;

top_depth = x./c1;

%plot(x,ur)

figure(1)

plot(top_depth(1:8000), norm_rad_strain, 'LineWidth', 2)

%tang_strain=ur./x;

%norm_tang_strain=tang_strain./Normal;

%plot(top_depth,norm_tang_strain)

% for i=1:max(size(ur))
%  if ((i < 4000) || (i > 6000))
%    plot((ur(i)/r(i))/(a^2*P)/(4*G*c1.^2), top_depth(i))
% end
%end

%----------------------------------------------------------------------------------
% Tidy Graph
%----------------------------------------------------------------------------------

% Set Graph Title in fontsize
title('Radial Pipe Strain (Lisowski)', 'FontSize', 12, 'FontName', 'Arial');

% Set Axis
xlabel('Distance (Top Depths)', 'FontSize', 12, 'FontName', 'Arial')
ylabel('Radial Strain (strain units)', 'FontSize', 12)

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
print(1,'-djpeg','lisowski_pipe_rad','-r500')

%----------------------------------------------------------------------------------