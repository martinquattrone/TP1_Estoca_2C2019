SNR_DB = linspace (-5 , 30 , 100) ;
SNR = 10.^(SNR_DB/10);
n = 1 : 4 : 25;

for i=1:7
    prob_digital = 1/2*(1-((1-2*qfunc(SNR.^(1/2))).^n(i)));
    leyenda = "n=" + n(i);
    semilogy(SNR_DB, prob_digital,'DisplayName', leyenda, 'Linewidth', 1.2)
    hold on
end

grid on
lgd = legend;
legend('Location','southwest')
xlim([-5 15])
ylim([10^-6 1]) 
xlabel('SNR [dB]')
ylabel('log(P_{e})')