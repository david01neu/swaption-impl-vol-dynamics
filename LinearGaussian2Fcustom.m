% Manually extended version of https://de.mathworks.com/help/fininst/lineargaussian2f.html
classdef LinearGaussian2Fcustom
%LINEARGAUSSIAN2F Create a 2 factor additive Gaussian interest rate model
%
% Syntax:
%
%   OBJ = LinearGaussian2F(ZeroCurve,a,b,sigma,eta,rho)
%
% Description:
%
%   Create a 2 factor additive linear Gaussian interest rate model by specifying
%   the zero curve, a, b, sigma, eta and rho parameters for the following
%   equation
%
%   r(t) = x(t) + y(t) + phi(t)
%   dx(t) = -a*x(t)dt + sigma*dW_1(t), x(0) = 0
%   dy(t) = -b*y(t)dt + eta*dW_2(t), y(0) = 0
%   dW_1(t)*dW_2(t) = rho
%
% Properties:
%
%   ZeroCurve - IRDataCurve or RateSpec. This is the zero curve that is
%               used to evolve the path of future interest rates.
%
%   a - Mean reversion for 1st factor; specified either as a scalar or as
%       a function handle which takes time as as input and returns a scalar
%       mean reversion value.
%
%   b - Mean reversion for 2nd factor; specified either as a scalar or as
%       a function handle which takes time as as input and returns a scalar
%       mean reversion value.
%
%   sigma - Volatility for 1st factor; specified either as a scalar or as
%           a function handle which takes time as as input and returns a
%           scalar mean volatility.
% 
%   eta - Volatility for 2nd factor; specified either as a scalar or as a
%         function handle which takes time as as input and returns a scalar
%         mean volatility.
%
%   rho - Scalar correlation of the factors.
%
% Methods:
%
%   [ZeroRates, ForwardRates] = simTermStructs(nPeriods)
%   [ZeroRates, ForwardRates] = simTermStructs(nPeriods,'name1','val1')
%
% Example:
%
%   Settle = datenum('15-Dec-2007');
%   CurveTimes = [1:5 7 10 20]';
%   ZeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
%   CurveDates = daysadd(Settle,360*CurveTimes,1);
%
%   irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
%   
%   a = .07;
%   b = .5;
%   sigma = .01;
%   eta = .006;
%   rho = -.7;
%
%   G2PP = LinearGaussian2F(irdc,a,b,sigma,eta,rho);
% 
%   SimPaths = G2PP.simTermStructs(10,'nTrials',100);
%
% Reference:
%
%   [1] Brigo, D and F. Mercurio. Interest Rate Models - Theory and
%   Practice. Springer Finance, 2006.
%
% See also LIBORMARKETMODEL, HULLWHITE1F

% Copyright 1999-2012 The MathWorks, Inc.
    
    properties
        ZeroCurve
        a
        b
        sigma
        eta
        rho
        
    end
    properties (Access = private)
        SDE
        PM
    end
    
    methods (Access = public)
        function obj = LinearGaussian2Fcustom(inCurve,ina,inb,insigma,ineta,inrho)
        %LINEARGAUSSIAN2F Create a 2 factor additive Gaussian interest rate model
            
            narginchk(6,6);
            
            % If RateSpec, convert to be an IRDataCurve
            if isafin(inCurve,'RateSpec')
                obj.ZeroCurve = IRDataCurve('Zero',inCurve.ValuationDate,inCurve.EndDates,...
                    inCurve.Rates,'Basis',inCurve.Basis,'Compounding',inCurve.Compounding);
            elseif isa(inCurve,'IRDataCurve') || isa(inCurve,'IRFunctionCurve')
                %obj.ZeroCurve = inCurve;
                obj.ZeroCurve = IRDataCurve('Zero',inCurve.Settle,inCurve.Settle+1:inCurve.Settle+365*30,...
                    inCurve.getZeroRates(inCurve.Settle+1:inCurve.Settle+365*30),'Basis',inCurve.Basis,'Compounding',inCurve.Compounding);
            else
                error(message('fininst:LinearGaussian2F:invalidCurve'));
            end
            obj.ZeroCurve = inCurve;

            
            % If scalar, convert to be a function handle
            if isa(ina,'function_handle')
                obj.a = @(t,V) ina(t);
            elseif isscalar(ina)
                obj.a = @(t,V) ina;
            else
                error(message('fininst:LinearGaussian2F:invalidA'));
            end
            
            % If scalar, convert to be a function handle
            if isa(inb,'function_handle')
                obj.b = @(t,V) inb(t);
            elseif isscalar(inb)
                obj.b = @(t,V) inb;
            else
                error(message('fininst:LinearGaussian2F:invalidB'));
            end
            
            % If scalar, convert to be a function handle
            if isa(insigma,'function_handle')
                obj.sigma = @(t,V) insigma(t);
            elseif isscalar(insigma)
                obj.sigma = @(t,V) insigma;
            else
                error(message('fininst:LinearGaussian2F:invalidSigma'));
            end
            
            % If scalar, convert to be a function handle
            if isa(ineta,'function_handle')
                obj.eta = @(t,V) ineta(t);
            elseif isscalar(ineta)
                obj.eta = @(t,V) ineta;
            else
                error(message('fininst:LinearGaussian2F:invalidEta'));
            end
            
            if isscalar(inrho)
                obj.rho = inrho;
            else
                error(message('fininst:LinearGaussian2F:invalidRho'))
            end
            
            obj.PM = @(t) obj.ZeroCurve.getDiscountFactors(daysadd(obj.ZeroCurve.Settle,round(365*t),3))';
            
            obj.SDE = hwv(@(t,X) diag([obj.a(t);obj.b(t)]),zeros(2,1), ...
                @(t,X) diag([obj.sigma(t);obj.eta(t)]),'Correlation',...
                [1 inrho;inrho 1],'StartState',[0;0]);
            
        end
        function [ZeroRates, ForwardRates, ShortRates, ZeroCouponBond] = simTermStructs(obj,nPeriods,varargin)
        %SIMTERMSTRUCTS Simulate Term Structures
        %
        % Syntax:
        %
        %   [ZeroRates] = simTermStructs(nPeriods)
        %   [ZeroRates, ForwardRates] = simTermStructs(nPeriods,'name1','val1')
        %
        % Description:
        %
        %   Simulate future zero curve paths using the specified 2 factor
        %   additive Gaussian interest rate model
        %
        % Required Input Arguments:
        %
        %   nPeriods - Number of simulation periods
        %
        % Option Input Arguments:
        %
        %   deltaTime - scalar time step betwen periods. Default is 1.
        %
        %   nTrials - scalar number of trials. Default is 1
        % 
        %   antithetic - Boolean scalar flag indicating whether antithetic
        %                sampling is used to generate the Gaussian random
        %                variates that drive the zero-drift, unit-variance
        %                rate Brownian vector dW(t). See
        %                hwv/simBySolution for more information.
        %
        %   Z - Direct specification of the dependent random noise process
        %       used to generate the zero-drift, unit-variance rate Brownian
        %       vector dW(t) that drives the simulation. See
        %       hwv/simBySolution for more information.
        %
        %   Tenor - numeric vector of maturities to be computed at each
        %           time step. Default is the tenor of the object's zero
        %           curve.
        % 
        % Output Arguments:
        %
        %   ZeroRates - nPeriods X nTenors X nTrials matrix of simulated
        %               zero rate term structures.
        %
        %   ForwardRates - nPeriods X nTenors X nTrials matrix of simulated
        %               forward rate term structures.
        %
        % Example:
        %
        %   CurveTimes = [1:5 7 10 20]';
        %   ZeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
        %   
        %   a = .07;
        %   b = .5;
        %   sigma = .01;
        %   eta = .006;
        %   rho = -.7;
        %
        %   G2PP = LinearGaussian2F([CurveTimes ZeroRates],a,b,sigma,eta,rho);
        % 
        %   SimPaths = G2PP.simTermStructs(10,'nTrials',100);
            
            narginchk(2,12);
            
            p = inputParser;
            
            p.addParameter('ntrials',1);
            p.addParameter('deltatime',1);
            p.addParameter('tenor',[]);
            p.addParameter('antithetic',false);
            p.addParameter('Z',[]);
            
            try
                p.parse(varargin{:});
            catch ME
                newMsg = message('fininst:LinearGaussian2F:optionalInputError');
                newME = MException(newMsg.Identifier, getString(newMsg));
                newME = addCause(newME, ME);
                throw(newME)
            end
            
            nTrials = p.Results.ntrials;
            deltaTime = p.Results.deltatime;
            Tenor = p.Results.tenor;
            Antithetic = p.Results.antithetic;
            Z = p.Results.Z;
            
            if isempty(Tenor)
                Tenor = yearfrac(obj.ZeroCurve.Settle,obj.ZeroCurve.Dates,obj.ZeroCurve.Basis);
                InitZeroRates = obj.ZeroCurve.Data;
                InitForwardRates = obj.ZeroCurve.getForwardRates(obj.ZeroCurve.Dates);
            else
                TenorDates = daysadd(obj.ZeroCurve.Settle,round(365*Tenor),3);
                InitZeroRates = obj.ZeroCurve.getZeroRates(TenorDates);
                InitForwardRates = obj.ZeroCurve.getForwardRates(TenorDates);
            end
            
            Tenor = Tenor(:);
            
            % Generate factors and short rates
            [Paths,SimTimes] = obj.SDE.simBySolution(nPeriods,'NTRIALS',nTrials,...
                'DeltaTime',deltaTime,'antithetic',Antithetic,'Z',Z);
            
            nTenors = length(Tenor);
            
            % Allocate interest rate paths
            ShortRates = zeros(nPeriods+1,nTrials);
            ZeroCouponBond = zeros(nPeriods+1,nTrials);
            ZeroRates = zeros(nPeriods+1,nTenors,nTrials);
            ForwardRates = ZeroRates;
            
            ZeroRates(1,:,:) = repmat(InitZeroRates',[1 1 nTrials]);
            ForwardRates(1,:,:) = repmat(InitForwardRates',[1 1 nTrials]);
            
            % Formula for V
            V = @(t,T) obj.sigma(t)^2/obj.a(t)^2*(T - t + 2/obj.a(t)*exp(-obj.a(t)*(T-t)) - 1/(2*obj.a(t))*exp(-2*obj.a(t)*(T-t)) - 3/2/obj.a(t)) + ...
                obj.eta(t)^2/obj.b(t)^2*(T - t + 2/obj.b(t)*exp(-obj.b(t)*(T-t)) - 1/(2*obj.b(t))*exp(-2*obj.b(t)*(T-t)) - 3/2/obj.b(t)) + ...
                2*obj.rho*obj.sigma(t)*obj.eta(t)/(obj.a(t)*obj.b(t))*(T - t + (exp(-obj.a(t)*(T-t)) - 1)/obj.a(t) + (exp(-obj.b(t)*(T-t)) - 1)/obj.b(t) - ...
                (exp(-(obj.a(t) + obj.b(t))*(T-t)) - 1)/(obj.a(t) + obj.b(t)));
            
            A = @(t,T) obj.PM(T)./obj.PM(t) .*exp(1/2*(V(t,T) - V(0,T) + V(0,t)));
            
            B = @(z,t,T) bsxfun(@rdivide,(1 - exp(-bsxfun(@times,z,(T-t)))),z);
            
            ZR = @(t,T,x,y) bsxfun(@rdivide,bsxfun(@plus,-log(A(t,T)),...
                bsxfun(@times,B(obj.a(t),t,T),x) + ...
                bsxfun(@times,B(obj.b(t),t,T),y)),T-t);

            ZCB = @(t,T,x,y) A(t,T) .* exp(-bsxfun(@times,B(obj.a(t),t,T),x)) .* exp(-bsxfun(@times,B(obj.b(t),t,T),y));

            phi = @(T) obj.ZeroCurve.getZeroRates(daysadd(obj.ZeroCurve.Settle,round(365*T),3))' + ...
                obj.sigma(T)^2/(2*obj.a(T)^2) * (1-exp(-obj.a(T)*T))^2 + ...
                obj.eta(T)^2/(2*obj.b(T)^2) * (1-exp(-obj.b(T)*T))^2 + ...
                obj.rho*(obj.sigma(T)*obj.eta(T))/(obj.a(T)*obj.b(T)) * (1-exp(-obj.a(T)*T)) * (1-exp(-obj.b(T)*T));
            
            for iPeriod=2:nPeriods+1
                t = SimTimes(iPeriod);
                ShortRates(iPeriod,:) = Paths(iPeriod,1,:) + Paths(iPeriod,2,:) + phi(t);
                ZeroRates(iPeriod,:,:) = ZR(t,t+Tenor',Paths(iPeriod,1,:),Paths(iPeriod,2,:));
                ZeroCouponBond(iPeriod,:) = ZCB(t,t+10,Paths(iPeriod,1,:),Paths(iPeriod,2,:));
                DF = exp(-bsxfun(@times,ZeroRates(iPeriod,:,:),Tenor'));
                ForwardRates(iPeriod,:,:) = bsxfun(@rdivide,-log(cat(2,DF(1,1,:),...
                    DF(1,2:end,:)./DF(1,1:end-1,:))),[Tenor(1) diff(Tenor')]);
            end
            
        end
    end
end