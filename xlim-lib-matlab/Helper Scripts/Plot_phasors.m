data = importdata("SIM2/pest.txt");
fileID= fopen('Simulation_Parameters.txt');
W = importdata('Simulation_Parameters.txt');
%      W = fscanf(fileID, '%s');
% end
fclose(fileID);
GD = data(1,1);
SD = data(1,2);
GDf = data(2,1);
SDf = data(2,2);
th = linspace(0, pi, 100);
R = 1;  %or whatever radius you want
x = R*cos(th);
y = R*sin(th);
plot(x,y);
hold on
    plot(GD,SD, '.', 'MarkerSize',20)
    plot(GDf,SDf, '.', 'MarkerSize',20)