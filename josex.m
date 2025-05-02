%% Encontrar ángulos

% Parámetros del sistema
L1 = 9; 
L2 = 6.5;
L3 = 3;

% Sistema de ecuaciones
syms t1 t2 t3

F = [L1*cos(t1) + L2*cos(t1+t2) + L3*cos(t1+t2+t3);
     L1*sin(t1) + L2*sin(t1+t2) + L3*sin(t1+t2+t3)];

% Condición inicial
t0 = [20;0;0];
pos0 = double(subs(F, [t1, t2, t3], t0.'));

% Restricciones de angulos
% 0 < t1 < 180
% -120 < t2 < 120
% -120 < t3 < 120

% Restricciones de coordendas
% Lo que alcance el brazo con los angulos

%% Coordenadas letras J y L

% Coordenadas letra J

coordJ = [
          11.7 13.7;
          11.7 15.7;
          11.7 17.7;
          11.7 19.7;
          9.7 19.7;
          9.7 17.7;
          9.7 15.7;
          9.7 13.7;
          9.2 12.2;
          7.7 11.7;
          5.7 11.7;
          5.7 9.7;
          7.7 9.7;
          9.7 10.2;
          11.2 11.7;
         ];

% Coordenadas letra L

coordL = [
          9 11;
          9 13;
          9 15;
          9 17;
          7 17;
          7 15;
          7 13;
          7 11;
          7 9;
          9 9;
          11 9;
          11 11;
         ];

%% Encontrar ternas de angulos para J y L

% Ternas de angulos J

xJ = coordJ';
tJ = zeros(3,15);

for i = 1:15
    Fi = F - xJ(:,i);
    tJ(:,i) = NR(Fi, t0, 100);
end

tJ = tJ';

% Ternas de angulos L

xL = coordL';
tL = zeros(3,12);

for i = 1:12
    Fi = F - xL(:,i);
    tL(:,i) = NR(Fi, t0, 100);
end

tL = tL';

%% Funciones

function[x,dif,k] = NR(F,x0,N)
    dif = 1;
    eps = 10^-6;
    x = x0;
    k = 1;

    var = transpose(symvar(F));
    J = jacobian(F,var);
    fx = matlabFunction(F, 'Vars', {var});
    jx = matlabFunction(J, 'Vars', {var});

    while (k <= N) && (dif > eps)
        x0cell = {x0};
        delta_x = -pinv(jx(x0cell{:}))*fx(x0cell{:});
        x = x0 + delta_x;
        dif = abs(max(x0-x)) / abs(max(x));
        x0 = x;
        k = k + 1;
    end
end