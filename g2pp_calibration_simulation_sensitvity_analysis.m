%% spot rate curve and ATM strikes for reference date 2023-12-29
zeroCurve = IRFunctionCurve('zero', today, @(t) nss(t));

% strike rates on 2023-12-29 
strikes = [2.5613	2.6332	2.5359	2.4048	2.2945; ...
    2.7157	2.5203	2.3373	2.2037	2.1067; ...
    2.2971	2.109	1.9874	1.9048	1.8292; ...
    1.8987	1.8058	1.7444	1.6795	1.6056; ...
    1.7034	1.6554	1.5917	1.5139	1.4486]./100;
%% 2023-12-29: Base swaption prices
marketPrices = [0.03246178 0.05811404 0.07602049 0.09012628 0.10177115;
    0.03496414 0.05984678 0.07717243 0.09102828 0.10255426;
    0.03266004 0.05468460 0.06973371 0.08257292 0.09264322;
    0.02920470 0.04841808 0.06123786 0.07301530 0.08177447;
    0.02562538 0.04215041 0.05372694 0.06400217 0.07157762];

% perform 10 G2++ calibrations with CRSP2 independent from each other
params_base = zeros(10, 5);
for i = 1:10
    params_base(i,:) = g2pp_calibration_global(strikes, marketPrices);
end
%% 2023-12-29: Stressed swaption prices (PCA AR(1) approach)
marketPrices = [0.06575717 0.11266907 0.14708596 0.17379744 0.19576285; ...
    0.05711917 0.09456879 0.11982600 0.13780037 0.15379605; ...
    0.04901588 0.07827618 0.09824440 0.11227347 0.12362716; ...
    0.04243872 0.06634451 0.08083884 0.09432507 0.10251506; ...
    0.03838713 0.05700220 0.06967552 0.07892543 0.08552567];

% perform 10 G2++ calibrations with CRSP2 independent from each other
params_scenario_ar1 = zeros(10, 5);
for i = 1:10
    params_scenario_ar1(i,:) = g2pp_calibration_global(strikes, marketPrices);
end
%% 2023-12-29: Stressed swaption prices (PCA Copula approach)
marketPrices = [0.06815114 0.11293252 0.14617452 0.17077131 0.19091385; ...
    0.05616029 0.08937374 0.11107814 0.12478307 0.13747674; ...
    0.04872270 0.07509529 0.08822732 0.09702614 0.10479228; ...
    0.04402125 0.06481151 0.07420473 0.08158876 0.08640673; ...
    0.03982943 0.05521803 0.06414588 0.06966038 0.07444727];

params_scenario_copula = zeros(10, 5);
for i = 1:10
    params_scenario_copula(i,:) = g2pp_calibration_global(strikes, marketPrices);
end
%% 2023-12-29: Stressed swaption prices (naive approach)
marketPrices = [0.06699010 0.10911751 0.13998694 0.16377518 0.18258327; ...
    0.05738687 0.09077289 0.11396113 0.13221598 0.14639228; ...
    0.04946029 0.07501207 0.09144963 0.10768927 0.11947116; ...
    0.04150511 0.06238198 0.07511011 0.08997146 0.09968488; ...
    0.03533866 0.05225467 0.06374968 0.07572664 0.08380653];

% perform 10 G2++ calibrations with CRSP2 independent from each other
params_scenario_naive = zeros(10, 5);
for i = 1:10
    params_scenario_naive(i,:) = g2pp_calibration_global(strikes, marketPrices);
end
%% RMSE computation for HW2 model parameters
marketPrices = [0.03246178 0.05811404 0.07602049 0.09012628 0.10177115;
    0.03496414 0.05984678 0.07717243 0.09102828 0.10255426;
    0.03266004 0.05468460 0.06973371 0.08257292 0.09264322;
    0.02920470 0.04841808 0.06123786 0.07301530 0.08177447;
    0.02562538 0.04215041 0.05372694 0.06400217 0.07157762];
optionPeriods = [5, 10, 15, 20, 25]; % Vector of option periods (1 year, 2 years, etc.)
swapTenors = [5, 10, 15, 20, 25]; % Vector of swap tenors (5 years, 10 years, etc.)

% compute RMSE for 10 G2++ parametrizations
for i = 1:10
    [error] = objectiveFunction(params_base(i,:), marketPrices, strikes, optionPeriods, swapTenors);
    params_base(i,6) = sqrt(error);
    disp(['RMSE: ', num2str(sqrt(error))]);
end
%% relative change: model vs observed swaption prices
marketPrices = [0.06688755 0.11392315 0.14851787 0.17517356 0.19707893; ...
    0.05747183 0.09452069 0.11939097 0.13677555 0.15234458; ...
    0.04936917 0.07837002 0.09731878 0.11053303 0.12134414; ...
    0.04306623 0.06661137 0.08031630 0.09285455 0.10050322; ...
    0.03897012 0.05714653 0.06925108 0.07790107 0.08422140];
optionPeriods = [5, 10, 15, 20, 25]; % Vector of option periods (e.g., 1 year, 2 years, etc.)
swapTenors = [5, 10, 15, 20, 25]; % Vector of swap tenors (e.g., 5 years, 10 years, etc.)
% G2++ parametrization (a,b,simga,eta,rho)
params = [0.699764813801145	0.0889206838827151	0.0521808414517144	0.0255999363951097	0.990000000000000];

for i = 1
    [error, matrix] = objectiveFunction(params, marketPrices, 0, optionPeriods, swapTenors);
    disp(['RMSE: ', num2str(sqrt(error))]);
end
%% G2++ interest rate term structure for different times to maturity
zeroCurve = IRFunctionCurve('zero', today, @(t) nss(t));
zeroCurve = toRateSpec(zeroCurve, today+1:today+465*50);


tmp = [0.8755    0.0461    0.0734    0.0243    0.9900];
a = tmp(1);
b = tmp(2);
sigma = tmp(3);
eta = tmp(4);
rho = tmp(5);

G2PP = LinearGaussian2F(zeroCurve,a,b,sigma,eta,rho);
nDates = 360;
DT = 1/12;
nTrials = 100;
tenor = (1:10)';
Paths = simTermStructs(G2PP, nDates,'ntrials',nTrials,'deltatime',DT,'tenor',tenor); 

trialIdx = 1;
figure
surf(tenor,0:nDates,Paths(:,:,trialIdx))
datetick y keepticks keeplimits
title(['Evolution of the Zero Curve for Trial:' num2str(trialIdx) ' of G2++ Model'])
xlabel('Tenor (Years)')
%% Monte Carlo simulation of BEL, guarantee and surplus for HW2/G2++ calibrations
K = 100000;
guaranteed_rate=0.01;
n_years = 30;
contract_value_0 = 10000;

simulated_data = zeros(8, 3, K);
simulated_data_summary = zeros(8, 6);
simulated_no_surplus_freq = zeros(8,30);

investment_strategies = ["1yr", "10yr"];
scenarios = ["base", "stress_ar1", "stress_copula", "stress_simple"];
counter_i = 1;
for investment_strategy = investment_strategies
  for scenario = scenarios
    scenario_params = [0 0 0 0 0];
    if strcmp(scenario, "base")
      scenario_params = [0.981500000000000	0.0490000000000000	0.0104000000000000	0.0123000000000000	0.914700000000000];
    elseif strcmp(scenario, "stress_ar1")
      scenario_params = [0.9725    0.0878    0.0831    0.0258    0.9900];
    elseif strcmp(scenario, "stress_copula")
      scenario_params = [0.1123    0.6351    0.0287    0.0480    0.9900];
    elseif strcmp(scenario, "stress_simple")
      scenario_params = [0.0824000000000000	0.451200000000000	0.0215000000000000	0.0386000000000000	0.990000000000000];
    end
    scenario_params
    
    [bel, guaranteed_benefit, surplus_benefit, simulated_no_surplus_frequencies] = pv_liabilities(scenario_params, n_years, contract_value_0, K, investment_strategy);
    simulated_data(counter_i,1,:) = bel;
    simulated_data(counter_i,2,:) = guaranteed_benefit;
    simulated_data(counter_i,3,:) = surplus_benefit;
    simulated_data_summary(counter_i,1) = mean(bel);
    simulated_data_summary(counter_i,2) = mean(guaranteed_benefit);
    simulated_data_summary(counter_i,3) = mean(surplus_benefit);
    simulated_data_summary(counter_i,4) = sqrt(var(bel));
    simulated_data_summary(counter_i,5) = sqrt(var(guaranteed_benefit));
    simulated_data_summary(counter_i,6) = sqrt(var(surplus_benefit));
    simulated_no_surplus_freq(counter_i,:) = simulated_no_surplus_frequencies; 

    counter_i = counter_i + 1;
  end
end
simulated_data_summary

%% sensitivity analysis of G2++ model parameters
hw2_a_sensitivity = hw2_sensitivity("a");
hw2_b_sensitivity = hw2_sensitivity("b");
hw2_sigma_sensitivity = hw2_sensitivity("sigma");
hw2_eta_sensitivity = hw2_sensitivity("eta");
hw2_rho_sensitivity = hw2_sensitivity("rho");

%% plot results of sensitivity analysis of G2++ model parameters
plot_hw2_change_of_params(hw2_a_sensitivity, ' a');
plot_hw2_change_of_params(hw2_rho_sensitivity, ' \rho');
plot_hw2_change_of_params(hw2_b_sensitivity, ' b');
plot_hw2_change_of_params(hw2_eta_sensitivity, ' \eta');
plot_hw2_change_of_params(hw2_rho_sensitivity, ' \rho');

%%
function [] = plot_hw2_change_of_params(hw2_sensitivity, param_string)
    scenario_identifiers = ["Base scenario", "AR(1) stress", "Copula stress", "Naive stress"];
    % plot the lines
    figure;
    hold on;
    for i = 1:4
        xaxis = -0.5:0.1:0.5;
        if param_string == " \rho"
            xaxis = -1:0.2:1;
        end
            
        plot(xaxis, hw2_sensitivity(:,i,1), 'DisplayName', scenario_identifiers(i));
    end
    hold off;

    % customize the plot
    xlabel(strcat('Relative change of parameter ', param_string));
    ylabel('BEL');
    title('BEL estimates for changes in HW2 model parameters');
    legend show; % Display the legend
    grid on;
end

function [result] = hw2_sensitivity(param_name)
    n_years = 30; % contract duration
    n_trials = 10000; % number of simulations
    premium_val = 10000; % single upfront premium
    investment_strategy = "10yr";

    % (# change values) x (# scenarios) x (# key values: bel, guarantee, surplus)
    result = zeros(5, 4, 8);
    
    j = 1;
    for scenario = ["base", "stress_ar1", "stress_copula", "stress_simple"]
        scenario_params = [0 0 0 0 0];
        if strcmp(scenario, "base")
          scenario_params = [0.981500000000000	0.0490000000000000	0.0104000000000000	0.0123000000000000	0.914700000000000];
        elseif strcmp(scenario, "stress_ar1")
          scenario_params = [0.699764813801145	0.0889206838827151	0.0521808414517144	0.0255999363951097	0.990000000000000];
        elseif strcmp(scenario, "stress_copula")
          scenario_params = [0.1123    0.6351    0.0287    0.0480    0.9900];
        elseif strcmp(scenario, "stress_simple")
          scenario_params = [0.0824000000000000	0.451200000000000	0.0215000000000000	0.0386000000000000	0.990000000000000];
        end

        i = 1;
        for change_rel = -0.5:0.1:0.5
            scenario_params_tmp = scenario_params;
            switch param_name
                case "a"
                    scenario_params_tmp = [scenario_params(1) * (1+change_rel) scenario_params(2:5)];
                case "b"
                    scenario_params_tmp = [scenario_params(1) scenario_params(2) * (1+change_rel) scenario_params(3:5)];
                case "sigma"
                    scenario_params_tmp = [scenario_params(1:2) scenario_params(3) * (1+change_rel) scenario_params(4:5)];
                case "eta"
                    scenario_params_tmp = [scenario_params(1:3) scenario_params(4) * (1+change_rel) scenario_params(5)];
                case "rho"
                    scenario_params_tmp = [scenario_params(1:4) change_rel*2];
            end
    
            [bel, guaranteed_benefit, surplus_benefit] = pv_liabilities(scenario_params_tmp, n_years, ...
                premium_val, n_trials, investment_strategy);
            result(i, j,:) = [mean(bel), median(bel), var(bel)];
            i = i + 1;
        end
        j = j + 1;
    end
end

%% stylistic insurance liability projection and valuation
function [val] = contract_value_increase(contract_value, guaranteed_rate, realized_rate)
    val = contract_value .* max(guaranteed_rate, 0.9 * realized_rate);
end

function [bel, guaranteed_benefit, surplus_benefit, simulated_no_surplus_frequencies] = pv_liabilities(params, n_years, ...
    premium_val, n_trials, investment_strategy)
  a = params(1);
  b = params(2);
  sigma = params(3);
  eta = params(4);
  rho = params(5);

  if nargin < 5
      investment_strategy = "1yr";
  end

  guaranteed_rate = 0.01;
  
  contract_value = premium_val;
  guaranteed_benefit = premium_val;
  surplus_benefit = 0;

  % number of simulations with surplus == 0
  simulated_no_surplus_frequencies = zeros(1,30);

  tenor = 1;
  if strcmp(investment_strategy, "10yr")
    tenor = 10;
  end

  rng(1) % set seed
  zeroCurve = IRFunctionCurve('zero', today, @(t) nss(t));
  G2PP = LinearGaussian2Fcustom(zeroCurve,a,b,sigma,eta,rho);
  [zero_rates, forward_rates, short_rates] = G2PP.simTermStructs(n_years*12,'ntrials',n_trials,'deltatime',1/12,'tenor',(1:10)'); 

  for i = 1:n_years
    realized_rate = 0;
    if strcmp(investment_strategy, "1yr")
      realized_rate = squeeze(zero_rates((i-1)*12+1,1,:));
      investment_sum = contract_value;
    elseif strcmp(investment_strategy, "10yr")
      realized_rate = squeeze(zero_rates(floor(i/tenor)*tenor*12+1,1,:));
      if mod(i-1, tenor) == 0
          investment_sum = contract_value;
      end
    end

    if realized_rate <= guaranteed_rate
        realized_rate
        simulated_no_surplus_frequencies(1,i) = simulated_no_surplus_frequencies(1,i) + 1;
    end
        
    cashflow_increase = contract_value_increase(investment_sum, guaranteed_rate, realized_rate);
    contract_value = contract_value + cashflow_increase;
  
    guaranteed_benefit = guaranteed_benefit + guaranteed_benefit * guaranteed_rate;
    surplus_benefit = surplus_benefit + (cashflow_increase - guaranteed_benefit * guaranteed_rate);
  end
  
  discount_factor = exp(- sum(short_rates)/12);
  bel = discount_factor(1,:) .* contract_value';%contract_value';
  guaranteed_benefit = discount_factor(1,:) .* guaranteed_benefit;
  surplus_benefit = discount_factor(1,:) .* surplus_benefit';
end

%% Nelson Siegel Svennson yield curve on 2023-12-29
% https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html
function [nss] = nss(t)
    beta0 =	2.466781;
    beta1 =	1.763840;
    beta2 =	-3.719427;
    beta3 =	2.955184;
    tau1 =	1.592413;
    tau2 =	15.593772;
    frac1 = (1-exp(-t./tau1)) ./ (t./tau1);
    frac2 = (1-exp(-t./tau2)) ./ (t./tau2);
    nss = (beta0 + beta1 .* (frac1) + beta2 .* (frac1 - exp(-t./tau1)) + beta3 .* (frac2 - exp(-t./tau2)))./100;
end
%%
function SwaptionPrice = swaptionbylg2f(inCurve,a,b,sigma,eta,rho,X,...
    ExerciseDate,Maturity,varargin)
    %SWAPTIONBYLG2F Compute swaption price using LG2F model
    %
    % Syntax:
    %
    %   Price = swaptionbylg2f(ZeroCurve,a,b,sigma,eta,rho,Strike,ExerciseDate,Maturity)
    %   Price = swaptionbylg2f(ZeroCurve,a,b,sigma,eta,rho,Strike,ExerciseDate,Maturity,...
    %                           'name1','val1')
    %
    % Description:
    %
    %   Compute swaption price for 2 factor additive Gaussian interest rate
    %   model given zero curve, a, b, sigma, eta and rho parameters for the following
    %   equation
    %
    %   r(t) = x(t) + y(t) + phi(t)
    %   dx(t) = -a*x(t)dt + sigma*dW_1(t), x(0) = 0
    %   dy(t) = -b*y(t)dt + eta*dW_2(t), y(0) = 0
    %   dW_1(t)*dW_2(t) = rho
    %
    % Inputs:
    %
    %   ZeroCurve - IRDataCurve or RateSpec. This is the zero curve that is
    %               used to evolve the path of future interest rates.
    %
    %   a - Mean reversion for 1st factor; specified as a scalar.
    %
    %   b - Mean reversion for 2nd factor; specified as a scalar.
    %
    %   sigma - Volatility for 1st factor; specified as a scalar.
    %
    %   eta - Volatility for 2nd factor; specified as a scalar.
    %
    %   rho - Scalar correlation of the factors.
    %
    %   Strike - NumSwaptions x 1 vector of Strike values.
    %
    %   ExerciseDate - NumSwaptions x 1 vector of serial date numbers or date strings
    %                  containing the swaption exercise dates.
    %
    %   Maturity - NumSwaptions x 1 vector of serial date numbers or date strings
    %              containing the swap maturity dates.
    %
    % Optional Inputs:
    %
    %   Reset - NumSwaptions x 1 vector of reset frequencies of swaption -- default is 1.
    %
    %   Notional - NumSwaptions x 1 vector of notional values of swaption -- default is 100.
    %
    %   OptSpec - NumSwaptions x 1 cell array of strings 'call' or 'put'. A call
    %             swaption entitles the buyer to pay the fixed rate. A put
    %             swaption entitles the buyer to receive the fixed rate.
    %             Default is call.
    %
    % Example:
    %
    %   Settle = datenum('15-Dec-2007');
    %
    %   ZeroTimes = [3/12 6/12 1 5 7 10 20 30]';
    %   ZeroRates = [0.033 0.034 0.035 0.040 0.042 0.044 0.048 0.0475]';
    %   CurveDates = daysadd(Settle,360*ZeroTimes,1);
    %
    %   irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
    %
    %   a = .07;
    %   b = .5;
    %   sigma = .01;
    %   eta = .006;
    %   rho = -.7;
    %
    %   Reset = 1;
    %   ExerciseDate = daysadd(Settle,360*5,1);
    %   Maturity = daysadd(ExerciseDate,360*[3;4],1);
    %   Strike = .05;
    %
    %   Price = swaptionbylg2f(irdc,a,b,sigma,eta,rho,Strike,ExerciseDate,Maturity,'Reset',Reset)
    %
    % Reference:
    %
    %   [1] Brigo, D and F. Mercurio. Interest Rate Models - Theory and
    %   Practice. Springer Finance, 2006.
    %
    % See also LINEARGAUSSIAN2F
    
    % Copyright 2012-2020 The MathWorks, Inc.
    
    narginchk(9, 15);
    
    if ~isafin(inCurve,'RateSpec') && ~isa(inCurve,'IRDataCurve') && ~isa(inCurve,'IRFunctionCurve') && ~isa(inCurve,'ratecurve')
        error(message('fininst:swaptionbylg2f:invalidCurve'));
    end
    
    if ~isscalar(a),error(message('fininst:swaptionbylg2f:invalidA')),end
    if ~isscalar(b),error(message('fininst:swaptionbylg2f:invalidB')),end
    if ~isscalar(sigma),error(message('fininst:swaptionbylg2f:invalidSigma')),end
    if ~isscalar(eta),error(message('fininst:swaptionbylg2f:invalidEta')),end
    if ~isscalar(rho),error(message('fininst:swaptionbylg2f:invalidRho')),end
    
    ExerciseDate = datenum(ExerciseDate);
    Maturity = datenum(Maturity);
    
    if any(ExerciseDate > Maturity)
        error(message('fininst:swaptionbylg2f:MaturityBeforeExercise'));
    end
    
    p = inputParser;
    
    p.addParameter('reset',2);
    p.addParameter('notional',1);
    p.addParameter('optspec',{'put'});
    
    try
        p.parse(varargin{:});
    catch ME
        newMsg = message('fininst:swaptionbylg2f:optionalInputError');
        newME = MException(newMsg.Identifier, getString(newMsg));
        newME = addCause(newME, ME);
        throw(newME)
    end
    
    Reset = p.Results.reset;
    Notional = p.Results.notional;
    
    if ischar(p.Results.optspec) || isstring(p.Results.optspec)
        OptSpec = cellstr(p.Results.optspec);
    elseif iscell(p.Results.optspec)
        OptSpec = p.Results.optspec;
    else
        error(message('fininst:swaptionbylg2f:invalidOptSpec'));
    end
    
    if ~all(ismember(lower(OptSpec),{'call','put'}))
        error(message('fininst:swaptionbylg2f:invalidOptSpec'));
    end
    
    try
        [X, ExerciseDate, Maturity, Reset, Notional, OptSpec] = finargsz(1, X(:), ExerciseDate(:),Maturity(:),...
            Reset(:), Notional(:), OptSpec(:));
    catch ME
        throwAsCaller(ME)
    end
    
    w = double(strcmpi(OptSpec,'call'));
    w(w == 0) = -1;
    
    if isafin(inCurve,'RateSpec')
        Settle = inCurve.ValuationDate;
        PM = @(t) intenvget(intenvset(inCurve,'EndDates',daysadd(Settle, round(365*t),3)),'Disc')';
    elseif isa(inCurve,'IRCurve') || isa(inCurve,'IRFunctionCurve')
        Settle = inCurve.Settle;
        PM = @(t) inCurve.getDiscountFactors(daysadd(Settle, round(365*t),3))';
    else
        Settle = inCurve.Settle;
        if any(inCurve.Basis == [0 2 3 8 9 10 12])
            PM = @(t) discountfactors(inCurve,daysadd(Settle, round(365*t),inCurve.Basis))';
        else
            PM = @(t) discountfactors(inCurve,daysadd(Settle, round(360*t),inCurve.Basis))';
        end
    end
    
    V = @(t,T) sigma^2/a^2*(T - t + 2/a*exp(-a*(T-t)) - 1/(2*a)*exp(-2*a*(T-t)) - 3/2/a) + ...
        eta^2/b^2*(T - t + 2/b*exp(-b*(T-t)) - 1/(2*b)*exp(-2*b*(T-t)) - 3/2/b) + ...
        2*rho*sigma*eta/(a*b)*(T - t + (exp(-a*(T-t)) - 1)/a + (exp(-b*(T-t)) - 1)/b - ...
        (exp(-(a + b)*(T-t)) - 1)/(a + b));
    
    A = @(t,T) PM(T)./PM(t) .*exp(1/2*(V(t,T) - V(0,T) + V(0,t)));
    
    B = @(z,t,T) (1 - exp(-bsxfun(@times,z,(T-t))))/z;
    
    nSwaptions = length(Maturity);
    SwaptionPrice = zeros(nSwaptions,1);
    
    optOptions = optimset('Jacobian','on','display','off');
    
    for swapidx=1:nSwaptions
        T = yearfrac(Settle,ExerciseDate(swapidx),1);
        Tenor = yearfrac(ExerciseDate(swapidx),Maturity(swapidx),1);
        
        ti = T:1/Reset(swapidx):(Tenor + T);
        tau = diff(ti);
        c = X(swapidx)*tau;
        c(end) = c(end) + 1;
        ti(1) = [];
        
        ux = -(sigma^2/a^2 + rho*sigma*eta/a/b)*(1 - exp(-a*T)) + ...
            sigma^2/(2*a^2)*(1 - exp(-2*a*T)) + ...
            rho*sigma*eta/(b*(a+b))*(1 - exp(-b*T - a*T));
        
        uy = -(eta^2/b^2 + rho*sigma*eta/a/b)*(1 - exp(-b*T)) + ...
            eta^2/(2*b^2)*(1 - exp(-2*b*T)) + ...
            rho*sigma*eta/(a*(a+b))*(1 - exp(-b*T - a*T));
        
        sigx = sigma*sqrt((1-exp(-2*a*T))/2/a);
        sigy = eta*sqrt((1-exp(-2*b*T))/2/b);
        rhoxy = rho*sigma*eta/((a+b)*sigx*sigy)*(1-exp(-(a+b)*T));
        
        x = linspace(ux - 10*sigx,ux + 10*sigx,1001)';
        
        cA = c.*A(T,ti);
        [ybar,~,exitflag] = fsolve(@(ybar) localObjFcn(ybar,x,cA,...
            B(a,T,ti),B(b,T,ti)),-x,optOptions);
        
        if exitflag <= 0
            %error(message('fininst:swaptionbylg2f:rootFailure'));
            SwaptionPrice(swapidx) = 99;
            break
        end
        
        h1 = (ybar - uy)./(sigy*sqrt(1 - rhoxy^2)) - rhoxy*(x - ux)./(sigx*sqrt(1 - rhoxy^2));
        h2 = bsxfun(@plus,B(b,T,ti).*sigy*sqrt(1 - rhoxy^2),h1);
        
        lambda = bsxfun(@times,A(T,ti).*c,exp(-bsxfun(@times,B(a,T,ti),x)));
        
        k = bsxfun(@times,-B(b,T,ti),bsxfun(@plus,uy - .5*(1 - rhoxy.^2)*sigy^2*B(b,T,ti),...
            rhoxy*sigy*(x-ux)/sigx));
        
        Y = exp(-1/2*((x - ux)./sigx).^2)./(sigx*sqrt(2*pi)) .* ...
            (normcdf(-w(swapidx)*h1) - sum(lambda.*exp(k).*normcdf(-w(swapidx)*h2),2));
        
        TempVal = trapz(x,Y);
        
        SwaptionPrice(swapidx) = w(swapidx)*Notional(swapidx)*TempVal*PM(T);
    end
end

function [f,g] = localObjFcn(y,x,cA,Ba,Bb)
    % LOCALOBJFUN Local function for solving roots of equation
    
    tmpSum = bsxfun(@times,cA,exp(-bsxfun(@times,Ba,x) - bsxfun(@times,Bb,y)));
    f = sum(tmpSum,2) - 1;
    g = diag(-sum(bsxfun(@times,Bb,tmpSum),2));
end

function [calibratedParams] = g2pp_calibration_global(strikes, marketPrices)
    optionPeriods = [5, 10, 15, 20, 25]; % vector of option periods (1 year, 2 years, etc.)
    swapTenors = [5, 10, 15, 20, 25]; % vector of swap tenors (5 years, 10 years, etc.)
    
    % initial guesses for the G2++ model parameters [a, b, sigma, eta, rho]
    initialParams = [rand(), rand(), 0.1*rand(), 0.05*rand(), rand()*2-1];
    
    % bounds for the G2++ parameters
    lb = [0.01, 00.01, 0, 0, -0.99]; % Lower bounds for a, b, sigma, eta, rho
    ub = [1, 1, 0.2, 0.2, 0.99];  % Upper bounds for a, b, sigma, eta, rho
    
    % objective function to minimize (sum of squared differences)
    objFun = @(params) objectiveFunction(params, marketPrices, strikes, optionPeriods, swapTenors);
    
    % Calibrate the G2++ parameters using simulated annealing
    opt.algorithm = NLOPT_GN_CRS2_LM;
    opt.lower_bounds = lb;
    opt.upper_bounds = ub;
    opt.min_objective = objFun;
    opt.maxeval = 10000;
    calibratedParams = nlopt_optimize(opt, initialParams);
    error = objFun(calibratedParams);
    
    % Display the calibrated parameters
    disp('Calibrated G2++ Parameters:');
    disp(['a: ', num2str(calibratedParams(1))]);
    disp(['b: ', num2str(calibratedParams(2))]);
    disp(['sigma: ', num2str(calibratedParams(3))]);
    disp(['eta: ', num2str(calibratedParams(4))]);
    disp(['rho: ', num2str(calibratedParams(5))]);
    disp(['error', num2str(error)])
end

function [calibratedParams] = g2pp_calibration_local(initialParams, strikes, marketPrices)
    optionPeriods = [5, 10, 15, 20, 25]; % Vector of option periods (e.g., 1 year, 2 years, etc.)
    swapTenors = [5, 10, 15, 20, 25]; % Vector of swap tenors (e.g., 5 years, 10 years, etc.)
    
    % Bounds for the G2++ parameters
    lb = [ ...
        max(initialParams(1)-0.2,0.01), ...
        max(initialParams(2)-0.2,0.01), ...
        max(0.001), ...
        max(0.001), ...
        -1]; % Lower bounds for a, b, sigma, eta, rho
    ub = [ ...
        min(initialParams(1)+0.1,1), ...
        min(initialParams(2)+0.1,1), ...
        min(0.2), ...
        min(0.2), ...
        1];  % Upper bounds for a, b, sigma, eta, rho
    
    % Objective function to minimize (sum of squared differences)
    objFun = @(params) objectiveFunction(params, marketPrices, strikes, optionPeriods, swapTenors);
    
    % Calibrate the G2++ parameters using simulated annealing
    calibratedParams = fminsearchbnd(objFun, initialParams, lb, ub);
end

function [error, matrix] = objectiveFunction(params, marketPrices, strikes, optionPeriods, swapTenors)
    % Unpack G2++ model parameters
    a = params(1);
    b = params(2);
    sigma = params(3);
    eta = params(4);
    rho = params(5);
    
    % Other required parameters for pricing (e.g., rate curve, dates)
    zeroCurve = IRFunctionCurve('zero',today, @(t) nss(t));
    
    % Initialize vector to hold model prices
    modelPrices = zeros(length(optionPeriods), length(swapTenors));
    
    % Loop over each swaption and calculate the model price using G2++ model
    for i = 1:length(optionPeriods)
        for j = 1:length(swapTenors)
            optionMaturity = addtodate(today, optionPeriods(i), 'year'); % Option period
            swapTenor = addtodate(optionMaturity, swapTenors(j), 'year');         % Swap tenor
            
               numerator = zeroCurve.getDiscountFactors(optionMaturity) - zeroCurve.getDiscountFactors(swapTenor);
               denominator = 0.5*sum(zeroCurve.getDiscountFactors((optionMaturity+0.5):183:swapTenor));

               swap_rate = numerator / denominator;
             
            %marketVolas(i,j) = marketPrices(i,j) * sqrt(2*pi/optionPeriods(i)) / (0.5*denominator);

            % Pricing the swaption using swaptionbylg2f() with G2++ parameters
            modelPrices(i,j) = swaptionbylg2f(zeroCurve, a, b, sigma, eta, rho, ...
                swap_rate, optionMaturity, swapTenor, 'notional', 1, 'optspec','put');
            %modelVolas(i,j) = modelPrices(i,j) * sqrt(2*pi/optionPeriods(i)) / (denominator);
        end
    end
    
    % Calculate the sum of squared differences between market and model prices
    error = sum(ones(1,5)*((modelPrices - marketPrices).^2)) ./ (size(optionPeriods,2) * size(swapTenors,2));

    if isnan(error)
        error = 99;
    end
end