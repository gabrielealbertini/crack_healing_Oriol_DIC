function v = compute_speed_back_diff(space,time)
%Nn=100;
%time_up = linspace(time(1), ...
%time(end),length(time)*Nn);
%space_up = interp1(time,space,time_up,'pchip'); 
    
    vtmp = diff(space)./diff(time);
    v(1) = 0;
    v(2:numel(time))=vtmp;
    %v = interp1(time_up(1:end-1),v_up,time(1:end-1));
end
