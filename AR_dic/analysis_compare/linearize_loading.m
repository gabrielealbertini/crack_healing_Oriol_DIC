function [time_shimadzu_lin,stroke_shimadzu_lin,force_shimadzu_lin] = ...
        linearize_loading(time_shimadzu,stroke_shimadzu,force_shimadzu, ...
                          stroke_lin_start,stroke_lin_end,plot_fig)
    
    [~,lin_start] = min(abs(stroke_shimadzu-stroke_lin_start));
    [~,lin_end  ] = min(abs(stroke_shimadzu-stroke_lin_end));

    stiffness = (force_shimadzu(lin_end)-force_shimadzu(lin_start))/...
        (stroke_shimadzu(lin_end)-stroke_shimadzu(lin_start));

    intercept = force_shimadzu(lin_start)-stiffness* ...
        stroke_shimadzu(lin_start);

    stroke_shimadzu_lin = stroke_shimadzu;
    force_shimadzu_lin = force_shimadzu;
    time_shimadzu_lin = time_shimadzu;
    force_shimadzu_lin(1:lin_start) = stiffness.*stroke_shimadzu_lin(1:lin_start)+intercept;


    [~,flt_start] = min(abs(force_shimadzu_lin(1:lin_end)-0));



    stroke_shimadzu_lin = stroke_shimadzu_lin(flt_start:end);
    force_shimadzu_lin  =  force_shimadzu_lin(flt_start:end);
    time_shimadzu_lin = time_shimadzu_lin(flt_start:end);

    stroke_shimadzu_lin = stroke_shimadzu_lin-stroke_shimadzu_lin(1);

    if plot_fig
        figure
        hold on

        plot(stroke_shimadzu-stroke_shimadzu(flt_start),force_shimadzu)
        plot(stroke_shimadzu_lin,force_shimadzu_lin)
        plot([stroke_shimadzu(lin_end),stroke_shimadzu(lin_start)]-stroke_shimadzu(flt_start),...
             [force_shimadzu(lin_end),force_shimadzu(lin_start)],'.');
        xlabel('stroke (m)')
        ylabel('Force (N)')
        hold off
    end
end
