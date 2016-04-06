clc
clear all


load('HoustonWeatherData.mat');

delta_t = 1/32;

T_interior = 75;  %%°F


t_start = 42005; %OLE Date from Excel 
t_end = 42370;


alpha = 0.05;   %BTU/(hr-ft²-°F)
k = 0.852;      %BTU/hr-ft-°F
L_wall = 1; %Thickness of wall - ft

delta_x = 1/8; 

numNodes = L_wall/delta_x + 3;

T = zeros(numNodes,1);
Tnew = T;
T(:) = T_interior;

T(1) = Tdb(1,2);

Tss = zeros(4,1);
Tss(:) = T_interior;

h_conv_i = 1.47; %BTU/hr - ft² - °F
h_conv_o = 4;    %BTU/hr - ft² - °F

times = t_start:delta_t:t_end - delta_t;

Fo = alpha*delta_t/(delta_x*delta_x);

Bi_o = h_conv_o * delta_x / k;
Bi_i = h_conv_i * delta_x / k;

figure1 = figure(1);
h1 = plot(Tnew);
hold on;
h2 = plot([1,2,numNodes-1,numNodes],Tss,'Color','red');
set(gca, 'YLimMode', 'manual')
ylim([30 100]);

line([2,2],[0,100], 'Color', 'black');
line([numNodes-1,numNodes-1],[0,100], 'Color', 'black');

text(numNodes/2, 90, 'Steady-State','Color','red');
text(numNodes/2 ,95,'Dynamic','Color','blue');

%For Steady-State Calculations
Rint = 1/h_conv_i
Rwall = L_wall/k;
Rout = 1/h_conv_o;
Rtot = Rint + Rwall + Rout;

xlabel('Node');
ylabel('Temperature (°F)');

annotation('textbox',[0.4 0.2 0 0.1],'String','Wall','FitBoxToText','on');
annotation('textbox',[0.13 0.8 0.2 0.1],'Interpreter','latex','String','$$T_{\infty}$$','FitBoxToText','on');
annotation('textbox',[0.83 0.8 0.2 0.1],'Interpreter','latex','String','$$T_{i}$$','FitBoxToText','on');

%legend('Dynamic','Steady-State','Location','north');



for i = 1:length(times)
	Tnew(1) = interp1(Tdb(:,1),Tdb(:,2),times(i));
	Tnew(2) = 2*Fo*(T(3)+Bi_o*T(1))+(1-2*Fo-2*Bi_o*Fo)*T(2);
	for j = 3:numNodes-2
		Tnew(j) = Fo*(T(j+1)+T(j-1))+(1-2*Fo)*T(j);	
	end
	Tnew(numNodes-1) = 2*Fo*(T(numNodes-2)+Bi_o*T(numNodes))+(1-2*Fo-2*Bi_o*Fo)*T(numNodes-1);

	Tnew(numNodes) = T_interior;
	

	qSS = (Tnew(1) - T_interior) / Rtot;
	Tss(1) = Tnew(1);
	Tss(2) = Tnew(1) - qSS*Rout;
	Tss(3) = Tss(2) - qSS*Rwall;
	Tss(4) = T_interior;

    matlabDate = x2mdate(times(i),0);
              
	set(h1, 'YData', Tnew);
	set(h2, 'YData', Tss);

    title(datestr(matlabDate));
    drawnow;
	T = Tnew;
        
end
