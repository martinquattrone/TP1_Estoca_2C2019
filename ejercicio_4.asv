close all;
clear all;

%% Definicion de constantes
SNR = 0:30;
h = 0.1;
mu = 0;
sigma=1;
a = (10.^(SNR/10)*sigma/h.^2).^(1/2);
n_realizaciones = 100000;
n_etapas = 9;

figure(1)

%% Digital MONTE CARLO 
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

semilogy(SNR, mean_e, 'DisplayName', 'Monte Carlo Digital', 'Linewidth', 1.8, 'Linestyle', '--')
hold on

%% Digital teorico
SNR_DB = linspace (0 , 30 , 1000) ;
SNR = 10.^(SNR_DB/10);
n = 9;

prob_digital = 1/2*(1-((1-2*qfunc(SNR.^(1/2))).^n));
leyenda = "Teorico Digital n=" + n;
semilogy(SNR_DB, prob_digital,'DisplayName', leyenda, 'Linewidth', 1.2)
hold on



%% Analogico MONTE CARLO
SNR = 0:30;
a = (10.^(SNR/10)*sigma/h.^2).^(1/2);
SNR_Ndb = 10.^(SNR./10);
G = (1/h).*(((SNR_Ndb.^(-1)+1).^(-1)).^(1/2));
mean_e = zeros(length(G),1);
x_n = x_n .* 0;

for k = 1:length(a)
    x_n = (2*a(k)) * randi([0 1], n_realizaciones, 1) - a(k) ;
    x_n_sombrero = x_n;
    e_n = zeros(n_realizaciones,1);
    for i = 1:n_realizaciones
        for j = 1:(n_etapas-1)
            x_n_sombrero(i) = G(k)*(h*x_n_sombrero(i)+normrnd(mu,sigma));
        end
        x_n_sombrero(i) = (h*x_n_sombrero(i)+normrnd(mu,sigma));
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

semilogy(SNR, mean_e, 'DisplayName', 'Monte Carlo Analogico', 'Linewidth', 1.8, 'Linestyle', '--')
hold on

%% Teorico Analogico
SNR_DB = linspace (0 , 30 , 1000) ;
SNR = 10.^(SNR_DB/10);
n = 9;

k = 0 : (n-2);
SNR_n = (SNR* (SNR/(SNR+1))^(n-1) )/(1 + sum( (SNR/(SNR+1)).^(k+1) ));
prob_analogico = qfunc(SNR_n.^(1/2));
leyenda = "Teorico Analogico n=" + n;
semilogy(SNR_DB, prob_analogico,'DisplayName', leyenda, 'Linewidth', 1.2)
hold on

%% Configuraciones del grafico 1

grid on
lgd = legend;
legend('Location','southwest')
xlim([0 25])
ylim([10^-6 1]) 
xlabel('Signal Noise Ratio - SNR [dB]', 'fontweight','bold','fontsize',11)
ylabel('Error Probability', 'fontweight','bold','fontsize',11)


%% EJERCICIO 4.3
%% Calculo de funciones de densidad

figure(2)
h = 0.25;
n = 9;
a = 10;

var_w = 0.25;                    % se define la varianz ade los Wi
sigma_w = sqrt(var_w);

SNR = (power(a*h, 2)/var_w);    % se definen el valor de la ganancia
G = sqrt(power(SNR/(SNR+1),n-1));

mu_y = h * a * G;               % se define la media y la varianza de Y|X=A     
k = 0: (n-2);
sigma_y = sigma_w * sqrt(sum( (SNR/(SNR+1) ).^(k+1) ));

aux = linspace(-10,10,100000);     % vector auxiliar
size_l = size(aux);

pdfs_a = normpdf(aux,mu_y,sigma_y); % se obtienen las densidades 
pdfs_an = normpdf(aux,(-mu_y),sigma_y);

line_1 = plot(aux,pdfs_a,'DisplayName', leyenda, 'Linewidth', 2.2 );
hold on
line_2 = plot(aux,pdfs_an, 'DisplayName', leyenda, 'Linewidth', 2.2);
hold on

y_a  = normrnd( mu_y, sigma_y, size(aux));
y_an = normrnd(-mu_y, sigma_y, size(aux));

line_3 = histogram(y_a, 'Normalization', 'pdf');
line_4 = histogram(y_an, 'Normalization', 'pdf');
legend([line_1,line_2,line_3,line_4],{'f_{Y|X_1=A}','f_{Y|X_1=-A}','Histograma f_{Y|X_1=A}','Histograma f_{Y|X_1=-A}'});
hold on;

% Configruaciones del grafico
figure(2)
grid on
legend('Location','northwest')
xlim([-8 8])
ylim([0 0.35]) 
ylabel({'$f_{Y|X_1}$'},'Interpreter', 'latex', 'fontweight','bold','fontsize',16)





