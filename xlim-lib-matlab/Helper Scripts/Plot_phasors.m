data = importdata("C:\Users\qazisd\OneDrive - National Institutes of Health\qazisd\tfitz-tconn-presentation\l410\pest.txt");

GD = data(1,1);
SD = data(1,2);
GDf = data(2,1);
SDf = data(2,2);
th = linspace(0, pi, 100);
R = 0.5;  %or whatever radius you want
x = R*cos(th);
y = R*sin(th);
plot(x,y);
hold on
    plot(GD,SD, '.', 'MarkerSize',20)
    plot(GDf,SDf, '.', 'MarkerSize',20)