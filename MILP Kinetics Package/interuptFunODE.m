function [status] = interuptFunODE(Time,Conc,flag,interupt_time)
%Function interupts ODE after a certain amount of time.

persistent INIT_TIME;
status = 0;
switch(flag)
    case 'init'
        INIT_TIME = tic;
    case 'done'
        clear INIT_TIME;
    otherwise
        elapsed_time = toc(INIT_TIME);
        if elapsed_time > interupt_time
            clear INIT_TIME;
            str = sprintf('%.6f',elapsed_time);
            error('interuptFun:Interupt',...
                 ['Interupted integration. Elapsed time is ' str ' seconds.'])
             disp('error')
        end
end
end

