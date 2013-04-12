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

d = 10000;
a = 1000;    
P = 10e6;
nu = 0.25;
E = 20e9;
G = 8e9;
r = -25000:1:25000;

%calcuate ur from mogi equation
ur = P.*(a.^3).*(1-nu)./G.*(r./(r.^2+d.^2).^1.5);
uz = P.*(a.^3).*(1-nu)./G.*(d./(r.^2+d.^2).^1.5);
% Calculations

% Err = radial displacement/radial distance
Rstrain = diff(ur)./diff(r);

% E00 = radial displacement/radial distance
Tstrain = ur./r;

% normalising factor
%Normal = (3*P*(a^3))/(4*G*d^3);

norm_r_strain = Rstrain/((3*P*a^3)/(4*G*d^3));

norm_tang_strain = Tstrain/((3*P*a^3)/(4*G*d^3)); %ad

%Ezz = vertical strain
Vert_strain = (nu/(1-nu))*(Rstrain+Tstrain(:,1:max(size(Tstrain))-1));

%Normalised Ezz
norm_vert_strain = ((nu/(1-nu))*(Rstrain+Tstrain(:,1:max(size(Tstrain))-1)))/((3*P*a^3)/(4*G*d^3));
Vert_strain =  Vert_strain*-1;
norm_vert_strain = norm_vert_strain*-1;


%Volumetric strain = Err + E00 + Ezz

Vol_strain = (Vert_strain + Tstrain(:,1:max(size(Tstrain))-1) + Rstrain);
norm_vol_strain= (Vert_strain + Tstrain(:,1:max(size(Tstrain))-1) + Rstrain)/((3*P*a^3)/(4*G*d^3));

source_depth = r/d;

%plot (r(:,1:max(size(r))-1), Vol_strain);

%legend('Volumetric Strain')
%xlabel('Distance (source depth)')
%ylabel('Volumetric Strain')

%------------------------------------
%Call comsol benchmark data

spherray=load('sphere_bench2.txt', '-ascii');
dist=spherray(:,1);
udisp=spherray(:,3);
%wdisp=spherray(:,3);
volstrain=spherray(:,2);
tangstrain=udisp./dist;
radstrain=diff(udisp)./diff(dist);
vertstrain= (nu/(1-nu))*(tangstrain(1:100)+radstrain);
vstrain = (tangstrain(1:100)+radstrain + vertstrain(1:100));

%------------------------------------------------------------------------------
%Plot Data
%------------------------------------------------------------------------------

figure(1)
hold on
plot(r(1:50000),Vol_strain,'r','LineWidth',2)
for i=1:max(size(volstrain))
  if (mod(i,2.5) == 0)
    scatter(dist(i),volstrain(i),60,'g', 'fill')
  end
end


% Set Graph Title in fontsize
title('Sphere Benchmark at 10km depth', 'FontSize', 12, 'FontName', 'Arial');

% Set Axis
xlabel('Distance (meters)', 'FontSize', 12, 'FontName', 'Arial')
ylabel('Strain (strain units)', 'FontSize', 12)

% Set Graph Background Color
set(gca,'Color',[1 1 1]);


% Show Grid Lines
%grid minor
%grid
grid off

% Set Graph Limits
%xMin xMax yMin yMax
axis([-1e4 1e4 -0.5e-7 13e-7])
waitforbuttonpress()

%legend command
line_1_name = 'Analytical Model';
line_2_name = 'Numerical Model';
legend(line_1_name, line_2_name, 'Location','NorthEast')

% Save figure 1 to jpeg,called output_test.jpg
% at a resolution of 500 dots per inch
% text is (for commercial print) 300
% images are 2000
print(1,'-djpeg','output_sphere_10','-r500')


