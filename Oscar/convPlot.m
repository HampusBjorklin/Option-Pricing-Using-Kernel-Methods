function [] = convPlot(list, error, figNum)
    figure(figNum)
    loglog(list,error, "r."); hold on, grid on
    pfit = polyfit(log(list)',log(error(:)),1);
    % pfit(2) = pfit(2)/ 10;
    pval = polyval(pfit, log10(list));
    plot(list, 10.^(pval), 'b-');

    xlim([0.8*min(list), 1.2*max(list)])
    legend("Simulated Data", "Polyfit, k = " + num2str(pfit(1)))
end

