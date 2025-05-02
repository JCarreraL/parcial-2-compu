load("coords.mat")
load("sols_j.mat")

x = coords(:,1);
y = coords(:,2);
time_vec = 1:length(x);

x_newton_coeffs = newton_coeff(time_vec,x);
y_newton_coeffs = newton_coeff(time_vec,y);


t = linspace(1,length(x),1000);
x_newton = zeros(1,length(t));
y_newton = zeros(1,length(t));
for i = 1:length(t)
    x_newton(i) = horner_eval(x_newton_coeffs,time_vec,t(i));
    y_newton(i) = horner_eval(y_newton_coeffs,time_vec,t(i));
end

x_spline_coeffs = spline_comp(time_vec,x');
y_spline_coeffs = spline_comp(time_vec,y');


x_spline = zeros(1,length(t));
y_spline = zeros(1,length(t));
for i = 1:length(t)
    x_spline(i) = evaluate_spline(time_vec,x_spline_coeffs,t(i));
    y_spline(i) = evaluate_spline(time_vec,y_spline_coeffs,t(i));
end

spline_interpol = figure;
plot(x_spline,y_spline)
axis("equal")
title('Spline')

newton_interpol = figure;
plot(x_newton,y_newton)
axis("equal")
xlim([-7 1])
ylim([7 17])
title('Newton')

x_linear = interp1(time_vec,x,t,'linear');
y_linear = interp1(time_vec,y,t,'linear');

linear_interpol = figure;
plot(x_linear,y_linear)
title('Linear')
axis("equal")
xlim([-7 1])
ylim([7 17])

%% Funciones
function yval = evaluate_spline(x, coeffs, xval)
    % Find which interval xval is in
    i = find(xval >= x, 1, 'last');
    if i >= length(x)
        i = length(x) - 1;  % clamp to last valid interval
    end

    dx = xval - x(i);
    a = coeffs(i,1);
    b = coeffs(i,2);
    c = coeffs(i,3);
    d = coeffs(i,4);

    yval = a + b*dx + c*dx^2 + d*dx^3;
end

function coeffs = spline_comp(x,y)
    
    n = length(x);
    h = diff(x);
    A = zeros(n,n);
    A(1,1) = 1;
    A(n,n) = 1;
    for i = 2:n-1
        A(i,i-1) = h(i-1);
        A(i,i)   = 2*(h(i-1) + h(i));
        A(i,i+1) = h(i);
    end
    s = zeros(n,1);
    s(2:n-1) = 3 * ( (y(3:n)-y(2:n-1))./h(2:n-1) - (y(2:n-1)-y(1:n-2))./h(1:n-2) );

    c = A\s;
    
    a = y(1:n-1)';
    b = zeros(n-1,1);
    d = zeros(n-1,1);

    for i = 1:n-1
        b(i) = (y(i+1) - y(i))/h(i) - h(i)/3*(2*c(i) + c(i+1));
        d(i) = (c(i+1) - c(i)) / (3*h(i));
    end
    c = c(1:n-1);
    coeffs = [a,b,c,d];

end


function [Q, z] = hermite(x, y, dy)
    n = length(x);
    z = repelem(x,2);      
    f = repelem(y,2);      

    Q = zeros(2*n, 2*n);   
    Q(:,1) = f(:);        

    
    for i = 1:2*n
        if mod(i,2) == 1
            Q(i,2) = dy((i+1)/2);     
        else
            Q(i,2) = (Q(i,1) - Q(i-1,1)) / (z(i) - z(i-1));
        end
    end

   
    for j = 3:2*n
        for i = 1:(2*n - j + 1)
            Q(i,j) = (Q(i+1,j-1) - Q(i,j-1)) / (z(i+j-1) - z(i));
        end
    end
end


function pol = lagrange_interpolation(x,y)
    syms xv
    n = length(x);
    L = ones(n,1);
    
    for i = 1:n
        for j = 1:n
            if j ~= i
                L(i) = L(i)*(xv-x(j))/(x(i)-x(j));
            end
        end
    end

    pol = dot(y,L);

end

function a = newton_coeff(x,y)
    n = length(x);
    Q = zeros(n,n);
    Q(:,1) = y;
    for j = 2:n
        for i = 1:n-j+1
            Q(i,j) = (Q(i+1,j-1)-Q(i,j-1))/(x(i+j-1)-x(i));
        end
    end
    a = Q(1,:);
end

function eval = horner_eval(a,x,t)
    n = length(a);

    eval = a(n);
    for i = n-1:-1:1
        eval = a(i) + (t-x(i))*eval;
    end

end
