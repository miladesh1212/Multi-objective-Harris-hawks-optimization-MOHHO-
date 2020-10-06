

function PlotCosts(pop,rep)

    pop_costs=[pop.Cost];
    plot3(pop_costs(1,:),pop_costs(2,:), pop_costs(3,:),'k.', 'markersize', 12);
    hold on;
    
    rep_costs=[rep.Cost];
    plot3(rep_costs(1,:),rep_costs(2,:), rep_costs(3,:),'r*', 'markersize', 12);
    
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    zlabel('3^{nd} Objective');
    grid on;
    
    hold off;

end