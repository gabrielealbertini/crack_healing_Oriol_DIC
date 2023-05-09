function v = compute_speed(space,time,method)
    if strcmp(method,'back_diff')
        vtmp = diff(space)./diff(time);
        v(1) = 0;
        v(2:numel(time))=vtmp;
    elseif strcmp(method,'forward_diff')
        vtmp = diff(space)./diff(time);
        v(1:numel(time)-1)=vtmp;
        v(numel(time))=0;
    else
        Nn=100;
        time_up = linspace(time(1), ...
                           time(end),length(time)*Nn);
        space_up = interp1(time,space,time_up,method); 
        
        v_up=diff(space_up)./diff(time_up);
        
        %v = interp1(time_up(1:end-1),v_up,time(1:end-1));
        v = interp1(time_up(1:end-1),v_up,time);
    end
end