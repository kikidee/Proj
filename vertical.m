% Function to calculate the vertical surface strain 
% from a spherical source at depth

%--------------------------------------------------

%variables:

% d = source depth [m]
% a = source radius [m]
% P = source pressure difference from surround [Pa]
% nu = Poisson Ratio
% G = Shear Modulus (Rigidity) [Pa]
% E = Young's Modulus [Pa]
% r = radial distance on the surface [m]
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

%calcuate ur from mogi equation
ur = P.*(a.^3).*(1-nu)./G.*(r./(r.^2+d.^2).^1.5);

% Calculations

% Err = radial displacement/radial distance
Rstrain = diff(ur)./diff(r);

%         v
% ezz= - ----
%        1-v
%-----------------
% E00 = radial displacement/radial distance
Tstrain = ur./r;

% normalising factor
%Normal = (3*P*(a^3))/(4*G*d^3);

%norm_r_strain = Rstrain/((3*P*a^3)/(4*G*d^3));

%norm_tang_strain = Tstrain/((3*P*a^3)/(4*G*d^3)); %ad

%Ezz = vertical strain
Vert_strain = (nu/(1-nu))*(Rstrain+Tstrain(:,1:max(size(Tstrain))-1));

%Normalised Ezz
norm_vert_strain = ((nu/(1-nu))*(Rstrain+Tstrain(:,1:max(size(Tstrain))-1)))/((3*P*a^3)/(4*G*d^3));
Vert_strain =  Vert_strain*-1;
norm_vert_strain = norm_vert_strain*-1;

source_depth = r/d;

plot (source_depth(1:6000), norm_vert_strain, 'LineWidth', 2)

%----------------------------------------------------------------------------------------------------

% Set Graph Title in fontsize
title('Vertical Strain (Lisowski)', 'FontSize', 12, 'FontName', 'Arial');

% Set Axis
xlabel('Distance (Top Depth)', 'FontSize', 12, 'FontName', 'Arial')
ylabel('Vertical Strain (strain units)', 'FontSize', 12)

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
print(1,'-djpeg','lisowski_sphere_vert_strain','-r500')

%----------------------------------------------------------------------------------