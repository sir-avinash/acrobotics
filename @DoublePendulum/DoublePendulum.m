classdef DoublePendulum

properties
    % model parameters
    m1@double
    l1@double
    l1com@doublemoment2
    m2@double
    l2@double
    l2com@double
    
    % constants
    g = 9.81;
    e1 = [1;0;0];
    e2 = [0;1;0];
    e3 = [0;0;1];
    
    % bounds
    bounds@struct

    controller@function_handle
    animate@function_handle
    controlParams@struct

end

properties (Constant = true)
    nDof = 6;
    nAct = 6;
    nConst = 2;
end

methods
	% class constructor
	function obj = DoublePendulum(varargin)
        
        if nargin > 0
            params = varargin{1};
            if isfield(params, 'm1')
                obj.m1 = params.m1;
            else
                obj.m1 = 1;
            end

            if isfield(params, 'm2')
                obj.m2 = params.m2;
            else
                obj.m2 = 2;
            end
            
            if isfield(params, 'l1')
                obj.l1 = params.l1;
            else
                obj.l1 = 0.5;
            end
            
            if isfield(params, 'l1')
                obj.l1 = params.l1;
            else
                obj.l1 = 0.5;
            end
            
            if isfield(params, 'l1')
                obj.l1 = params.l1;
            else
                obj.l1 = 0.5;
            end 

            if isfield(params, 'l1')
                obj.l1 = params.l1;
            else
                obj.l1 = 0.5;
            end
            
            if isfield(params, 'JQ')	
                obj.JQ = params.JQ;
            else
                obj.JQ = obj.mQ*obj.lQ^2;
            end

            if isfield(params, 'lQ')
                obj.lQ = params.lQ;
            end
        end
                
    end     
    
    % set property values
    function obj = setProperty(obj, propertyName, value)
        if isprop(obj, propertyName)
            set(obj, propertyName, value);
        else
            error([propertyName ' not a property of class ',class(obj)]);
        end
    end

    function sol = simulate(obj, tspan, x0, solver)
        odefun = @(t,x)systemDynamics(obj, t, x);
        sol = solver(odefun, tspan, x0);
    end
    
    function [f, g] = quadVectorFields(obj, x)
%         params = [mQ;JQ;lQ;g];
        params = [obj.mQ;obj.JQ;obj.lQ;obj.g];
        [f,g] = quadrotorVectorFields(x,params);
    end
    
    function ddx = geometricDynamics(obj, t, x)
        u = obj.controller(obj, t,x);
        [fvec, gvec] = obj.quadVectorFields(x);
        ddx = fvec + gvec*u;
    end

    function ddx = eulerDynamics(obj, t, x)
        u = obj.controller(obj, t,x);
        [fvec, gvec] = obj.quadVectorFields(x);
        ddx = fvec + gvec*u;
    end
    
    function u = calcControlInput(obj, t, x)
        u = obj.controller(obj, t,x);
    end
    
    function [Ad, Bd] = discrLinearizeQuadrotor(obj, Ts, xk, uk)
        [Ad, Bd] = discretizeLinearizeQuadrotor(obj, Ts, xk, uk);
    end
    
    % TODO:
%     obj = setLoadMass(mL); 

    % animation
    animateQuadrotor(obj,t,x,varargin)
    
    
end

end