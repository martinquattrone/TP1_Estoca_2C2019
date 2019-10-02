close all;
clear all;

%% Definicion de constantes
SNR = 0:30;
h = 0.1;
mu = 0;
sigma=1;
a = (10.^(SNR/10)*sigma/h.^2).^(1/2);
n_realizaciones = 10000;
n_etapas = 9;

%% Digital monte carlo 
for k = 1:length(a)
    x_n = (2*a(k)) * randi([0 1], n_realizaciones, 1) - a(k) ;
    x_n_sombrero = x_n;
    e_n = zeros(n_realizaciones,1);
    for i = 1:n_realizaciones
        for j = 1:n_etapas
            x_n_sombrero(i) = h*x_n_sombrero(i)+ normrnd(mu,sigma);
            if x_n_sombrero(i)>=0
                x_n_sombrero(i)=  a(k);
            else
                x_n_sombrero(i)= -a(k);
            end
        end
        if x_n(i) ~= x_n_sombrero(i)
            e_n(i) = 1;
        end
    end
    mean_e(k) = mean(e_n);
end

semilogy(SNR, mean_e, 'DisplayName', 'Monte Carlo', 'Linewidth', 1.8, 'Linestyle', '--')
hold on

%% Digital teorico
SNR_DB = linspace (0 , 30 , 1000) ;
SNR = 10.^(SNR_DB/10);
n = 9;

prob_digital = 1/2*(1-((1-2*qfunc(SNR.^(1/2))).^n));
leyenda = "Digital n=" + n;
semilogy(SNR_DB, prob_digital,'DisplayName', leyenda, 'Linewidth', 1.2)
hold on

%% Teorico analogico
SNR_DB = linspace (0 , 30 , 1000) ;
SNR = 10.^(SNR_DB/10);
n = 9;

k = 0 : (n-2);
SNR_n = (SNR* (SNR/(SNR+1))^(n-1) )/(1 + sum( (SNR/(SNR+1)).^(k+1) ));
prob_analogico = qfunc(SNR_n.^(1/2));
leyenda = "Analogic n=" + n;
semilogy(SNR_DB, prob_analogico,'DisplayName', 'Teorico analogico', 'Linewidth', 1.2)
hold on

%% Caso analogico
SNR = 0:30;
a = (10.^(SNR/10)*sigma/h.^2).^(1/2);
SNR_Ndb = 10.^(SNR./10);
G = (1/h).*(((SNR_Ndb.^(-1)+1).^(-1)).^(1/2));
mean_e = zeros(length(G),1);
x_n = x_n .* 0;

%% Monte carlo analogico
for k = 1:length(a)
    x_n = (2*a(k)) * randi([0 1], n_realizaciones, 1) - a(k) ;
    x_n_sombrero = x_n;
    e_n = zeros(n_realizaciones,1);
    for i = 1:n_realizaciones
        for j = 1:(n_etapas-1)
            x_n_sombrero(i) = G(k)*(h*x_n_sombrero(i)+normrnd(mu,sigma));
        end
        x_n_sombrero(i) = (h*x_n_sombrero(i)+normrnd(mu,sigma));
        %x_n_sombrero(i) = x_n_sombrero(i)/G(k); %Divido por G por que tengo n-1 ganancias
        if x_n_sombrero(i)>=0
                x_n_sombrero(i)=  a(k);
            else
                x_n_sombrero(i)= -a(k);
        end
     
        if x_n(i) ~= x_n_sombrero(i)
            e_n(i) = 1;
        end
    end
    mean_e(k) = mean(e_n);
end

semilogy(SNR, mean_e, 'DisplayName', 'Monte Carlo', 'Linewidth', 1.8, 'Linestyle', '--')
hold on



%% Configuraciones del grafico
grid on
lgd = legend;
legend('Location','southwest')
xlim([0 30])
ylim([10^-6 1]) 
xlabel('Signal Noise Ratio - SNR [dB]', 'fontweight','bold','fontsize',11)
ylabel('Error Probability', 'fontweight','bold','fontsize',11)

