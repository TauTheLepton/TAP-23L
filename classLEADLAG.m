classdef classLEADLAG < handle
    % classPID - klasa prostego regulatora PID
    %
    % mjc - 2016
    % #####################################################################
    % Private internal
    properties(GetAccess = 'public', SetAccess = 'private')
        % Tuning parameters
        K = 1;     % Gain
        LEAD = 0;  % Lead time
        LAG = 0;   % Lag time
        Tp = 1;    % Sample time
        Hlim = 0; % high limit
        Llim = 0; % low limit               
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        K1 = 0;     %coefficients
        K2 = 0;
        K3 = 0;     
        oldIn = 0;        % internal state
    end
    
    % #####################################################################
    % PUBLIC atributes
    properties(GetAccess = 'public', SetAccess = 'public')
        in = 0;     %input
        out = 0;    %LEADLAG output

    end
    
    % PUBLIC methods ######################################################
    methods(Access = 'public')
        % constructor
        function obj = classLEADLAG(K, LEAD, LAG, Tp, Hlim, Llim) 
            % set internal parameters
            obj.reTune(K, LEAD, LAG, Tp, Hlim, Llim);
        end
        
        % #################################################################
        % calc LEADLAG: obj.calc(PV) or obj.calc(PV, SP)
        function output = calc(obj, varargin)
           switch nargin
               case 2 % one argument
                  obj.in = varargin{1};
               otherwise
                   error('error: wrong number of parameters');
           end

         
          
           %LEADLAG math
           obj.out = obj.K1 * obj.in + obj.K2 * obj.oldIn + obj.K3 * obj.out;
           obj.oldIn = obj.in; 
            
           % limits
           if obj.out > obj.Hlim
              obj.out = obj.Hlim;
           end
           if obj.out < obj.Llim
              obj.out = obj.Llim;
           end
           
          
           % aktualna wartoœæ zwracana jako return metody
           output = obj.out;
        end
        
        % PID tuning ######################################################
        function  reTune(obj, K, LEAD, LAG, Tp, Hlim, Llim)
            obj.K = K;
            obj.LEAD = LEAD;
            obj.LAG = LAG;
            obj.Tp = Tp;
            obj.Hlim = Hlim;
            obj.Llim = Llim;
           
            % coefficients recalculation
            obj.K1 = K * (Tp + 2*LEAD) / (Tp + 2*LAG);
            obj.K2 = K * (Tp - 2*LEAD) / (Tp + 2*LAG);
            obj.K3 = (2*LAG - Tp) / (2*LAG + Tp);  
        end
    end
end

