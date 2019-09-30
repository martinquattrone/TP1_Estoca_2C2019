close all;
%% Definiciones

SNR_DB = linspace (-5 , 30 , 1000) ;
SNR = 10.^(SNR_DB/10);
n = 1 : 4 : 25;

%% Probabilidad de error digital
for i=1:7
    prob_analogico = 1/2*(1-((1-2*qfunc(SNR.^(1/2))).^n(i)));
    leyenda = "Digital n=" + n(i);
    semilogy(SNR_DB, prob_analogico,'DisplayName', leyenda, 'Linewidth', 1.2)
    hold on
end

%% Probabilidad de error analogico

for i=1:7
    k = 0 : (n(i)-2);
    SNR_n = (SNR* (SNR/(SNR+1))^(n(i)-1) )/(1 + sum( (SNR/(SNR+1)).^(k+1) ));
    prob_analogico = qfunc(SNR_n.^(1/2));
    leyenda = "Analogic n=" + n(i);
    semilogy(SNR_DB, prob_analogico,'DisplayName', leyenda, 'Linewidth', 1.2,'Linestyle', '-.')
    hold on
end


%% Configuraciones del grafico
grid on
lgd = legend;
legend('Location','southwest')
xlim([-5 30])
ylim([10^-6 1]) 
xlabel('Signal Noise Ratio - SNR [dB]', 'fontweight','bold','fontsize',11)
ylabel('Error Probability', 'fontweight','bold','fontsize',11)