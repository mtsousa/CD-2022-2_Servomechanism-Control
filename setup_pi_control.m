% Autores: Isabela Maria Pereira Cruzeiro
%          Matheus Teixeira de Sousa
%
% Este código recebe os parâmetros a e b, o período de amostragem e os
% requisitos de projeto, e retorna as matrizes calculadas

function [A, B, C, D, K1, K2, Ke, Gaa, Gab, Gba, Gbb, Ha, Hb] = setup_pi_control(a, b, T, ts, Mp)
    % Malha aberta
    K_ma = 420;
    A = [0 1; 0 -a];
    B = [0; b];
    C = [1 0];
    D = 0;

    % Espaço de estados discretizado
    [G, H] = c2d(A, B, T);
    H = K_ma*H;
    
    % Calcula função de transferência
    [num, den] = ss2tf(G, H, C, D);

    % Transforma para canônica controlável
    M = [H G*H];
    W = [den(2) 1; 1 0];
    P = M*W;

    G_bar = (P\G)*P;
    H_bar = P\H;
    C_bar = C*P;
    D_bar = D;

    % Calcula os polos desejados em S
    sigma = 4/ts;
    omega_d = -sigma*pi/log(Mp);

    modulo_z = exp(-sigma*T);
    fase_z = omega_d*T;
    polos = [modulo_z*cos(fase_z)+modulo_z*sin(fase_z)*1j; modulo_z*cos(fase_z)-modulo_z*sin(fase_z)*1j];

    z = tf('z');
    den_z = (z)*(z-polos(1))*(z-polos(2));
    aux = get(den_z);
    coefs_cell = aux.Numerator();
    coefs = coefs_cell{1};

    alpha1 = coefs(2);
    alpha2 = coefs(3);
    alpha3 = coefs(4);

    % Calcula o ganho K_chapeu
    G_chapeu = [G_bar H_bar; 0 0 0];
    H_chapeu = [0; 0; 1];

    phi_G = G_chapeu^3 + alpha1*G_chapeu^2 + alpha2*G_chapeu + alpha3*eye(3);
    K_chapeu = [0 0 1]*inv([H_chapeu G_chapeu*H_chapeu (G_chapeu^2)*H_chapeu])*phi_G;

    % Calcula os ganhos K2 e K1
    N = [G_bar - eye(2) H_bar; C_bar*G_bar C_bar*H_bar];
    K2K1 = (K_chapeu + [0 0 1])/(N);
    K2_bar = [K2K1(1) K2K1(2)];
    K1_bar = K2K1(3);

    K2 = K2_bar/P;
    K1 = K1_bar;

    % Projeto do observador ordem mínima
    polo_min = 0.1;
    Gaa = G(1, 1);
    Gab = G(1, 2);
    Gba = G(2, 1);
    Gbb = G(2, 2);
    Ha = H(1, 1);
    Hb = H(2, 1);
    Ke = (Gbb - polo_min)/Gab;

end