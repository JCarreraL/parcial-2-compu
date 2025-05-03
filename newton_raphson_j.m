%% Soluciones (angulos)
clc, clear
load('coords.mat')

sols_j = [20,0,0;pseudo_newton(coords,1000)];


%% Funciones

% Newton-Raphson
function[sols] = pseudo_newton(coord,iteraciones)
    
    sols = zeros(length(coord),3);
    
    x0 = [20;0;0];
    
    for j = 1:length(coord)
        
        syms t1 t2 t3

        F = [(coord(j,1) - (9*(cosd(t1))+6.5*(cosd(t1+t2))+3*(cosd(t1+t2+t3))));
            (coord(j,2) - (9*(sind(t1))+6.5*(sind(t1+t2))+3*(sind(t1+t2+t3))))];

        x = num2cell(x0);


        J = jacobian(F);

        fx = matlabFunction(F);
        jx = matlabFunction(J);


        j_eval = jx(x{:});
        f_eval = fx(x{:});

        M = pinv(j_eval);
    
        x = num2cell(cell2mat(x)-(M*f_eval));
    
        for i = 1:iteraciones
            if (x{1}<0 || x{1}>180) || (x{2}<-120 || x{2}>120) || (x{3}<-120 || x{3}>120)
                x{1} = randi([0,180],1);
                x{2} = randi([-120,120],1);
                x{3} = randi([-120,120],1);
            end

            j_eval = jx(x{:});
            f_eval = fx(x{:});
    
            M = pinv(j_eval);
        
            x = num2cell(cell2mat(x)-(M*f_eval));
           
        end
        x = cell2mat(x);
        
        sols(j,:) = x;
        x0 = x;
    end
end
    