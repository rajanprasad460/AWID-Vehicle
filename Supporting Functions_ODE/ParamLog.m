classdef ParamLog < handle
    properties
        phi
        lambda = []
        Sys_Input = []
        n_wheel = []
        dis_fact
        R
        w
        t = []
        tau0 = [0;0;0;0;0;0]' % Torque of motors
        tau = [0;0;0;0;0;0]
    end
end