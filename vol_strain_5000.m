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
Vert_strain = Vert_strain*-1;
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

spherray=load('sphere_bench_5000m.txt', '-ascii');
dist=spherray(:,1);
udisp=spherray(:,3);
wdisp=spherray(:,4);
volstrain=spherray(:,2);
tangstrain=udisp./dist;
radstrain=diff(udisp)./diff(dist);
%vertstrain= (nu/(1-nu))*(tangstrain(1:100)+radstrain);
%vstrain = (tangstrain(1:100)+radstrain + vertstrain(1:100));

%------------------------------------------------------------------------------
%Plot Data
%------------------------------------------------------------------------------

figure(1)
hold on
plot(r,uz,'r','LineWidth',2)

x_fig = 0;
y_fig = 0;
spacing = 0.008;
for i=1:max(size(wdisp))
  x_new = dist(i)/max(dist);
  y_new = wdisp(i)/(1*max(wdisp));
  if (abs(x_new-x_fig)^2 + abs(y_new-y_fig)^2) > spacing
    x_fig = x_new;
y_fig = y_new;
    scatter(dist(i),wdisp(i),50,'g', 'fill')
  end
end

plot(dist,wdisp, 'g', 'LineWidth',2)


% Set Graph Title in fontsize
title('Sphere Benchmark for Vertical Displacement at 5km depth', 'FontSize', 10, 'FontName', 'Arial');

% Set Axis
xlabel('Distance (metres)', 'FontSize', 10, 'FontName', 'Arial')
ylabel('Displacement (metres)', 'FontSize', 10)

% Set Graph Background Color
set(gca,'Color',[1 1 1]);


% Show Grid Lines
%grid minor
%grid
grid off

% Set Graph Limits
%xMin xMax yMin yMax
%axis([-1.e4 1.e4 -0.015 0.015])
waitforbuttonpress()

%legend command
line_1_name = 'Analytical Model';
line_2_name = 'Numerical Model';
legend(line_1_name, line_2_name, 'Location','NorthWest')

% Save figure 1 to jpeg,called output_test.jpg
% at a resolution of 500 dots per inch
% text is (for commercial print) 300
% images are 2000
print(1,'-djpeg','output_sphere_5_vertdisp','-r500')

