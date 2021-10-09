%% Stephenson Type II linkage simulation using the Newton Raphson Method

% Author: José Canales Avellaneda
% Date: 08/10/2021
%% Variable initializations
clearvars; close all;clc;
a = 70;   
b = 100;  
c = 90;   
d = 110;  

e = 120;  
f=160;
g = 150;    

g2 = 150;    

gamma = 30*pi/180;
delta = 20*pi/180;

%%

N = 361;
theta1 = 0;
xA = [0;0]; % ground pin at A
xD = [d;0]; % ground pin at D

[theta3,theta4,theta5,theta6] = deal(zeros(1,N));

% initial guesses for Newton-Raphson algorithm
t3 = 1.1; t4 = 1.5; t5 =0.5; t6 = 4.7;

%Pre-allocation for optimization
AB = zeros(2,N);
BC = zeros(2,N);
AF = zeros(2,N);
DE = NaN(2,N);
FG = zeros(2,N);
EG = zeros(2,N);
theta2 = zeros(1,N);

%%
figure(1);
for k=1:N
    theta2(k) = (k-1)*(2*pi)/(N-1);
    [ea,na] = UnitVector(theta2(k));
    [ed,nd] = UnitVector(theta1);
    
    %% Newton Raphson Method
    for i = 1:100   % máximo 100 iteraciones 
        
        %First 4-bar linkage
        [eb,nb] = UnitVector(t3);
        [ec,nc] = UnitVector(t4);
        
        %Second 4-bar linkage
        [ee,ne] = UnitVector(t4 - gamma);
        [eg,ng] = UnitVector(t5);
        [ef,nf] = UnitVector(t6);
        
        %Ground references
        [eg2,ng2] = UnitVector(-delta);
        g3eg3=g2*eg2-d*ed;
        
        %Phi matrix
        phi(:,1) = [a*ea+b*eb-c*ec-d*ed;
                    g*ee+e*ef-f*eg-g3eg3];
        %Matrix with variables  
        % q = [theta3, theta4, theta5, theta6]
        %Jacobian matrix
        J = [b*nb,-c*nc,zeros(2,1),zeros(2,1);
             zeros(2,1), -g*ne, -f*ng, e*nf ]; % [4x4]
        
        dq = -J\phi;     % - inv(J)*phi
        t3 = t3 + dq(1);
        t4 = t4 + dq(2);
        t5 = t5 + dq(3);
        t6 = t6 + dq(4); 
        
        if norm(phi) < 0.00001
            disp(strcat('Convergió en la iteración:',num2str(i)));
            break
        end        
    end
    %%
    
    theta3(k) = t3;
    theta4(k) = t4;
    theta5(k) = t5;
    theta6(k) = t6;
    
    %% Simulation
    %plot_watt2(xA,xD,a,b,g2,g,f,e,ea,eb,eg2,eg,ef,ee,k);
    %Avoid using the plot function to draw the trajectory
    
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
    plot([xA(1) AB(1,k)],[xA(2) AB(2,k)]);
    hold on;
    plot([AB(1,k) BC(1,k)],[AB(2,k) BC(2,k)]);
    plot([AF(1,k) FG(1,k)],[AF(2,k) FG(2,k)]);
    plot([DE(1,k) EG(1,k)],[DE(2,k) EG(2,k)]);
    
    %3 pin bars
    fill([xA(1) xD(1) AF(1,k) xA(1)],...
        [xA(2) xD(2) AF(2,k) xA(2)],color);
    fill([xD(1) BC(1,k) DE(1,k) xD(1)],...
        [xD(2) BC(2,k) DE(2,k) xD(2)],color);
    
    %Drawing trajectory
    p1=plot(AB(1,:),AB(2,:),'.r');
    p2=plot(DE(1,:),DE(2,:),'g');
    p3=plot(FG(1,:),FG(2,:),'.b');
    
    p1.LineWidth=2;p2.LineWidth=2;p3.LineWidth=2;
    
    hold off;
    xlim([-100 350]);
    ylim([-50 250]);
    xlabel('x');
    ylabel('y');
    title('Watt Type 2 Linkage')
    grid on;
    axis equal;
    pause(0.01)
end

figure;
plot(theta3);hold on;
plot(theta4);plot(theta5);plot(theta6);
xlabel('N° samples')
ylabel('angle (°)')
axis tight
legend('\theta_3','\theta_4','\theta_5','\theta_6')
grid on