y = zeros(14,1);
x = zeros(14,1);
%load("sols_j.mat")

for i = 1:length(sols_j)
    x(i) = 9*cosd(sols_j(i,1)) + 6.5*cosd(sols_j(i,1)+sols_j(i,2)) + 3*cosd(sols_j(i,1)+sols_j(i,2)+sols_j(i,3));
    y(i) = 9*sind(sols_j(i,1)) + 6.5*sind(sols_j(i,1)+sols_j(i,2)) + 3*sind(sols_j(i,1)+sols_j(i,2)+sols_j(i,3));
end
x2 = linspace(-18.5,18.5,1000);
y2 = sqrt(18.5^2-x2.^2);

figure
hold on    
plot(x,y,"r")



plot(x2,y2,"b")
axis equal

hold off