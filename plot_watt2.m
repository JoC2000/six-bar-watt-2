function plot_watt2(xA,xD,a,b,g2,g,f,e,ea,eb,eg2,eg,ef,ee,k)
    set(groot,'defaultLineLineWidth',5.0)
    color = [135 206 250]./255;
    
    % Find the positions of points B,C,D,C2,D2
    AB(:,k) = FindPos( xA, a, ea);
    BC(:,k) = FindPos( AB(:,k), b, eb);
    AF(:,k) = FindPos( xA, g2, eg2);
    DE(:,k) = FindPos( xD, g, ee);
    FG(:,k) = FindPos( AF(:,k), f, eg);  
    EG(:,k) = FindPos( DE(:,k), e, ef);
    
    % Draw the lines
    plot([xA(1) AB(1,k)],[xA(2) AB(2,k)]);hold on;  
    plot([AB(1,k) BC(1,k)],[AB(2,k) BC(2,k)]);
    plot([AF(1,k) FG(1,k)],[AF(2,k) FG(2,k)]);
    plot([DE(1,k) EG(1,k)],[DE(2,k) EG(2,k)]);
    
    fill([xA(1) xD(1) AF(1,k) xA(1)],[xA(2) xD(2) AF(2,k) xA(2)],color);
    fill([xD(1) BC(1,k) DE(1,k) xD(1)],[xD(2) BC(2,k) DE(2,k) xD(2)],color);
      
    hold off;
    xlim([-100 350]);
    ylim([-50 250]);
    xlabel('x');
    ylabel('y');
    title('watt Type 2 Linkage')
    grid on;
    axis equal;
end

