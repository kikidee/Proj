close all
clear all
%---------------------------------------------
%Script to plot depth against author
%---------------------------------------------

label    = ['Berkerley et al (1998)''Mattioli et al(1998)''Shepherd et al (1998)''Mattioli & Herd (2003)''Widiwijayanti et al (2005)''Voight et al (2006)''Elsworth et al (2008)''Foroozan et al (2010)''Voight et al (2010b)''Mattioli et al (2010)''Hautmann et al (2010b)''Linde et al (2010)'];
spheresx = [2 3 4 5 6 12];
spheresy = [-6000 -750 -12000 -855 -6000 -6000];

figure(1)
hold on
scatter(spheresx,spheresy, 'filled', 'b', 'MarkerEdgeColor', [0 0 1])

petx     = [1];
pety     = [-5000];

scatter(petx, pety, 50,'d', 'blue', 'filled')

spheredxa = [7 7];
spheredya = [-12000 -5000];

plot(spheredxa, spheredya, '-oblue','MarkerFaceColor',[0 0 1],'MarkerEdgeColor', [0 0 1],'linewidth', 2 )

spheredxb = [8 8];
spheredyb = [-6000 -17000];

plot(spheredxb, spheredyb, '-oblue','MarkerFaceColor',[0 0 1],'MarkerEdgeColor', [0 0 1], 'linewidth', 2 )

prolx    = [9 10 11];
proly    = [-10400 -10400 -12250];

scatter(prolx, proly, 'blue', 'square', 'filled')

% Set Graph Title in fontsize
title('Changing opinions of Magma System Location', 'FontSize', 12, 'FontName', 'Arial');


% Set Y Axis
ylabel('Depth (meters)', 'FontSize', 12)

% Set Graph Background Color
set(gca,'Color',[1 1 1]);

% Set X Axis Number of Ticks
set(gca, 'XTick', 0:1:12)

% Define Rotation
text_rotation = 45;

% Define Labels
x_labels = {'';'Berkerley et al (1998)';'Mattioli et al(1998)';'Shepherd et al (1998)';'Mattioli & Herd (2003)';'Widiwijayanti et al (2005)';'Voight et al (2006)';'Elsworth et al (2008)';'Foroozan et al (2010)';'Voight et al (2010b)';'Mattioli et al (2010)';'Hautmann et al (2010b)';'Linde et al (2010)'};

% Get Current Tick Labels
tick_labels=get(gca,'XTickLabel');

% Remove Current Tick Labels
set(gca,'XTickLabel',[]);

% Get Tick Label Positions
x_label_positions=get(gca,'XTick');
c=get(gca,'YTick');

% Generate Y Positions
y_label_positions = repmat(  ...
	c(1)-.1*(c(2)-c(1)),       ... %{ Copy This Value (Y Position of Label}      }%
	length(x_label_positions), ... %{ Into An Array This Long (Number of Labels) }%
	1                          ... %{ With Width 1 (IE Array not Matrix)         }%
);
y_label_positions = y_label_positions - 500;

% Write New Tick Labels
text( ... 
  x_label_positions,                       ... %{ X Position of Text Label }%
  y_label_positions,                       ... %{ Y Position of Text Label }%
  x_labels,                                ... %{ Text To Be Printed       }%
  'HorizontalAlignment','right',           ... %{ Align Right              }%
  'rotation',text_rotation                 ... %{ Angle of Rotation        }%
);


% Show Grid Lines
%grid minor
%grid
grid off


% Set Graph Limits
%xMin xMax yMin yMax
axis([0 13 -18000 0])

% Set Whitespace around graph
% works from bottom left
% offset_y offset_x height_x height_y
% offset & height sums to 1 - we've left a 0.1 border at top right
set(gca, 'Position',[0.2 0.3 0.7 0.6])

waitforbuttonpress()

%legend command
%line_1_name = 'Petrological Constraint,';
%line_2_name = 'Single Sphere';
%line_3_name = 'Dual Sphere';
%line_4_name = 'Prolate';
%legend(line_2_name, line_1_name, line_3_name, line_4_name, 'Location','NorthEast')

%Save figure 1 to jpeg,called output_depths.jpg
% at a resolution of 500 dots per inch
% text is (for commercial print) 300
% images are 2000
print(1,'-djpeg','output_depth','-r500')