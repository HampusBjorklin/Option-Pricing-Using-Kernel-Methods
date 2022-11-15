function [price] = basketSolverSingle(x, y, Strike)

    Settle = datetime(2009,5,1);
    Maturity  = datetime(2012,5,1);
    
    % Define RateSpec
    Rate = 0.005;
    Compounding = -1;
    RateSpec = intenvset('ValuationDate', Settle, 'StartDates',...
    Settle, 'EndDates', Maturity, 'Rates', Rate, 'Compounding', Compounding);
    
    % Define the Correlation matrix. Correlation matrices are symmetric,
    % and have ones along the main diagonal.
    Corr = [1 0; 0 1];

     % Define BasketStockSpec
        AssetPrice =  [x;y]; 
        Volatility = [0.12;0.15];
        Quantity = [0.5;0.5];
        BasketStockSpec = basketstockspec(Volatility, AssetPrice, Quantity, Corr);
        
        % Compute the price of the call basket option
        OptSpec = {'call'};
        AmericanOpt = 0; % American option
        
        price = basketbyls(RateSpec, BasketStockSpec, OptSpec, Strike, Settle, Maturity,...
                    'AmericanOpt',AmericanOpt);

    end
    
