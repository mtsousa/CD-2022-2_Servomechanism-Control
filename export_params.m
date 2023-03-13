% Autores: Isabela Maria Pereira Cruzeiro
%          Matheus Teixeira de Sousa
%
% Este código salva os parâmetros de projeto em um arquivo .mat

% Estima os valores de a e b
[a, b] = fit_transfer_function();

% Define os requisitos de projeto
T = 0.01;
ts = 0.3;
Mp = 0.1;

% Calcula os ganhos do sistema
[A, B, C, D, K1, K2, Ke, Gaa, Gab, Gba, Gbb, Ha, Hb] = setup_pi_control(a, b, T, ts, Mp);

% Salva o workspace com as matrizes projetadas
save('params_ts03_mp10.mat')