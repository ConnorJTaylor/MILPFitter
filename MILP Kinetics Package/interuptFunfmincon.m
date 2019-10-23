function [stop] = interuptFunfmincon(Time,Conc,state,interupt_time_fmincon)
%Function interupts fmincon after a certain amount of time.

persistent INITIAL_TIME;
stop = false;
switch(state)
    case 'init'
        INITIAL_TIME = tic;
    case 'done'
        clear INITIAL_TIME
    otherwise
        elapsed_time = toc(INITIAL_TIME);
        if elapsed_time > interupt_time_fmincon
            clear INITIAL_TIME
            stop = true;
            disp('Optimisation stopped due to optimisation time greater than 200 seconds')
        end
end


