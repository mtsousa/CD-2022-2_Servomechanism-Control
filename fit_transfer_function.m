% Autores: Isabela Maria Pereira Cruzeiro
%          Matheus Teixeira de Sousa
%
% Este código calcula os parâmetros a e b a partir dos dados coletados

function [a, b] = fit_transfer_function()    
    % Lê os valores salvos
    load('real_data.mat');
    out_real = mean_output;

    % Pega o valor máximo do output
    [out_max, idx_max] = max(out_real);
    t_max = time(idx_max);

    % Calcula os polos do sistema
    omega_d = pi/t_max;
    sigma = -log((out_max - out_real(1000))/out_real(1000))/t_max;

    % Calcula o formato do denominador
    z = tf('z');
    den = (z + sigma)^2 + (omega_d)^2;
    aux = get(den);
    coefs_cell = aux.Numerator();
    coefs = coefs_cell{1};

    % Calcula o valor de a e b
    a = coefs(2);
    b = coefs(3)/420;

end