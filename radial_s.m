% Function to calculate the radial surface strain 
% from a spherical source at depth

%--------------------------------------------------

%variables:

% d = source depth [m]
% a = source radius [m]
% P = source pressure difference from surround [Pa]
% nu = Poisson Ratio
% G = Shear Modulus (Rigidity) [Pa]
% E = Young's Modulus [Pa]
% x = radial distance on the surface [m]
% ur = radial displacement

%----------------------------------------------------

clear all
close all
%Inputs

d = 5000;
a = 1000;    
P = 10e6;
nu = 0.25;
E = 20e9;
G = 8e9;
r = -15000:5:15000;
%r = sqrt(x.^2 + d.^2);


%calcuate ur from mogi equation
ur = P.*(a.^3).*(1-nu)./G.*(r./(r.^2+d.^2).^1.5);



% Calculations

% Err = radial displacement/radial distance
Rstrain = diff(ur)./diff(r);

% normalising factor
Normal = (3*P*(a^3))/(4*G*d^3);

norm_r_strain = Rstrain/((3*P*a^3)/(4*G*d^3));

source_depth = r/d;

plot(source_depth(1:6000), norm_r_strain, 'LineWidth', 2)

%-----------------------------------------------------------------------------------------

% Set Graph Title in fontsize
title('Radial Strain (Lisowski)', 'FontSize', 12, 'FontName', 'Arial');

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
print(1,'-djpeg','liswoski_sphere_radial','-r500')

%----------------------------------------------------------------------------------



