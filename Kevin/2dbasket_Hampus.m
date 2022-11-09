Settle = datetime(2009,5,1);
Maturity  = datetime(2012,5,1);

% Define RateSpec
Rate = 0.05;
Compounding = -1;
RateSpec = intenvset('ValuationDate', Settle, 'StartDates',...
Settle, 'EndDates', Maturity, 'Rates', Rate, 'Compounding', Compounding);

% Define the Correlation matrix. Correlation matrices are symmetric,
% and have ones along the main diagonal.
Corr = [1 0.50; 0.50 1];

% Define BasketStockSpec
AssetPrice =  [35;45]; 
Volatility = [0.12;0.15];
Quantity = [0.5;0.5];
BasketStockSpec = basketstockspec(Volatility, AssetPrice, Quantity, Corr);

% Compute the price of the call basket option
OptSpec = {'call'};
Strike = 42;
AmericanOpt = 1; % American option

M1 = 5;
M2 = 5;
m1_stock_prices = 10*(1:M1)
m2_stock_prices = 10*(1:M2)
price_mesh = zeros(M1, M2);

for i = 1:M1
    for j = 1:M2
        AssetPrice =  [m1_stock_prices(i); m2_stock_prices(j)]; 
        BasketStockSpec = basketstockspec(Volatility, AssetPrice, Quantity, Corr);
        price_mesh(i, j) = basketbyls(RateSpec, BasketStockSpec, OptSpec, Strike, Settle, Maturity,...
            'EuropeanOpt',AmericanOpt)
    end
end

mesh(m1_stock_prices, m2_stock_prices, price_mesh)

