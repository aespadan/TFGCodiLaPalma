%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clear clc
clear cfg
clear figures
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHYSICAL PARAMETERS


dt = 0.9; % In seconds
w=3.5; % Dyke halfwidth
fw=10; % Factor for computing width
dx=0.001;x=-fw*w:dx:fw*w; % Chart and axis width (times dyke width)
x=x';
Nx=length (x(:,1));

To=710;
Tc=1200;
Ts=950;
Tl=1250;
Tm=1100;

K=5.3*10^-7;
rho=2800;
c=1480;
L=4*10^5;
z = L/(rho*c*(Tl-Ts));

bltzm = 1.38*10^-23;
res0 = 231.73946;
E = 1.1473*10^-19;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL STATE

% Initial state controls
ini_Subpl = false;   % Plot initial state with the rest of selected plots in a single window
ini_Ind = false;   % Plot in. state independently

ini_xlim = [-fw*w, fw*w];  % Plot axis limits

ini_ylim = [0, 1375];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERMEDIATE STATES

% Intermediate states controls
int_Subpl = false;   % Plot intermed. states with rest of selected plots
int_Ind = false;   % Plot independently

int_xlim = [-fw*w, fw*w];  % Plot axis limits

int_ylim = [0, 1375];



% REINJECTION

%Reinjection cycles (No reinjection == 0)
reinject = 0; 

% REINJECTION TIMES (years) 
% One array element is needed for each reinjection specified above 
% Each element of the array is the time for the successive reinjection and 
% stops at that point from further evolution until reinjection is set. 
% If no reinjection is set the first element will be the last computed time
% and further evolution will be stopped.
time_reinject = [5, 10];


% Final time (years) (stopping time for last reinj. cycle)
time_limit = 1;


% Number of sequential figures for intermediate states (times of graphing chosen below)
f = 11;


% Corresponding graphing times for intermediat states (years)
% There must be as many susbsets as reinjection cycles chosen and in each of those, one time value for each figure specified above 
% If some graphing times are bigger than the succesive reinjection time, they won't be plotted
t = [0.0833  0.1667  0.25  0.3333  0.4167  0.5  0.5833  0.6667  0.75  0.833  0.9167];
%[0.083  0.125  0.167  0.2083  0.25  0.29167  0.333  0.375  0.4167]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEADY STATE

% Steady state controls
est_Analytic = false;    % Add analytic solution
est_Subpl = false;     % Plot intermed. states with rest of selected plots
est_Ind = false;    % Plot in. state independently

est_xlim = [-fw*w, fw*w];  % Plot axis limits

est_ylim = [0, 1375];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPLICIT FINITE DIFERENCE PARAMETERS

% Mathematical stability and limitations. Boundary and starting conditions

err_lim = 1E-8;

lambda1=K/(1+z)*dt/(dx^2);
lambda2=K*dt/(dx^2);

if lambda1>=0.5 
    error('Stability requirements not met. lambda1 parameter larger than 0.5')
end

if lambda2>=0.5 
    error('Stability requirements not met. lambda2 parameter larger than 0.5')
end


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for k=1:Nx
    T(k)= To;
    Tr(k)=To; %Second reinjection method temperatures
    T_s(k)= Ts; % Solidus temp. line
    T_l(k)= Tl; % Liquidus temp. line
    T_m(k)= Tm;
    Tlast(k)= 0; 
    Trlast(k)=0;
    res(k)=0;
end

% Starting figure
if ini_Subpl|int_Subpl|est_Subpl == 1
    figure(1)
    tiledlayout(reinject+1,f+2) 
end

figure(2)
%tiledlayout(1,2)

figure(3)
%tiledlayout(1,2)

% INITIAL STATE

r_t = 0; %Array index for selecting graph time


for R = 0:reinject 
    c_t =1;
    r_t = r_t+1;

% Temperature field for the initial state (or initial state after reinjection)
    
    for j=1:Nx
        if R==1
            if T(j)>Tm   % Mobile melt boundary selection
                Tr(j)=Tc;
            else
                Tr(j)=Tlast(j);
            end
        else
            if Tr(j)>Tm   % Mobile melt boundary selection
                Tr(j)=Tc;
            else
                Tr(j)=Tr(j);
            end
        end

        if abs(x(j,1))<w  % Initial state
            T(j)=Tc; 
        end
    end
    

    if ini_Subpl == true
        figure(1)  % Plotting initial state
        nexttile(R*(f+2)+1)
        if R<1
            plot(x,T','g', x,T_s,'b',x,T_l,'r');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Initial Temp.') ; xline([w -w],'-.',{'Right Wall','Left Wall'}); xlim(ini_xlim); ylim(ini_ylim);
            hold on
        else
            plot(x,T','g', x,Tr,'c', x,T_s,'b', x,T_l,'r');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Initial Temp. upon reinjection R = ' num2str(R)]) ; xline([w -w],'-.',{'Right Wall','Left Wall'}); xlim(ini_xlim); ylim(ini_ylim);
            hold on
        end
    end

    if ini_Ind == true  % (independently)
        figure()
        if R<1
            plot(x,T','g', x,T_s,'b',x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Initial Temp.') ; xline([w -w],'-.',{'Right Wall','Left Wall'}); legend ('Original dyke injection', 'T solidus', 'T liquidus','Location','eastoutside'); xlim(ini_xlim); ylim(ini_ylim);
        else
            plot(x,T','g', x,Tr,'c', x,T_s,'b', x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Initial Temp. upon reinjection R = ' num2str(R)]) ; xline([w -w],'-.',{'Right Wall','Left Wall'}); legend ('Original dyke reinjection', 'Mobile melt reinjection', 'T solidus', 'T liquidus','Location','eastoutside'); xlim(ini_xlim); ylim(ini_ylim);
        end
    end
    
    


   if R<1
        figure(2)
        %nexttile(1)
        plot(x,T');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Temperature evolution with time') ;  xlim(ini_xlim); ylim(ini_ylim);
        hold on
        for k = 1:Nx
            res(k) = 1/(res0*exp(-E/(bltzm*(T(k)+273.15))));
        end
        figure(3)
        %nexttile(2)
        plot(x,res);xlabel ('Distance (m)'); ylabel ('Resistivity') ; yscale log; title('Resistivity evolution with time') ;  xlim(ini_xlim); axis ([-fw*w fw*w 0.1 100]);
        hold on
   else
        figure(2)
        %nexttile(1)
        plot(x,Tr);xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Temperature evolution with time') ;  xlim(ini_xlim); ylim(ini_ylim);
        hold on
        for k = 1:Nx
            res(k) = 1/(res0*exp(-E/(bltzm*(T(k)+273.15))));
        end
        figure(3)
        %nexttile(2)
        plot(x,res);xlabel ('Distance (m)'); ylabel ('Resistivity'); title('Resistivity evolution with time') ;  axis ([-fw*w fw*w 0.1 100]);
        hold on
   end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate states                                                     

    if R<reinject
        maxiter = floor((time_reinject(R+1)*365*24*3600)/dt); 
    else
        maxiter = (time_limit*365*24*3600)/dt; 
    end

    pos = R*(f+2)+2; %Figure indexing for next graph
    

    for iter = 1:maxiter 
        Tlast=T;
        if R<1  
            Trlast=T;
        else
            Trlast=Tr;
        end

        for i=2:Nx-1
          
            depth_2D = (Tlast(i-1)-2*Tlast(i)+Tlast(i+1))/dx^2;  % PDE implementation
            if (Ts<Tlast(i))&&(Tlast(i)<Tl)
                time_1D = (K/(1+z))*depth_2D;
                T(i) = time_1D*dt + Tlast(i);
            else
               time_1D = K*depth_2D;
                T(i) = time_1D*dt + Tlast(i);
            end
            
            if R>0
                depth_2Dr = (Trlast(i-1)-2*Trlast(i)+Trlast(i+1))/dx^2; % PDE for second reinjection method
                if (Ts<Trlast(i))&&(Trlast(i)<Tl)
                    time_1Dr = (K/(1+z))*depth_2Dr;
                    Tr(i) = time_1Dr*dt + Trlast(i);
                else
                   time_1Dr = K*depth_2Dr;
                    Tr(i) = time_1Dr*dt + Trlast(i); 
                end
            end
        end
        
        if c_t <= f
            if iter > floor((t(r_t,c_t)*365*24*3600)/dt)
                tempstotal = iter*dt/3600; 
                anys=tempstotal/(24*365);
                sanys = num2str(anys);
            

% Intermediate state plotting (for times selected at the start)
                
                if R<1
                    % Adding plot to multiplot
                    if int_Subpl == true
                        figure(1)
                        nexttile(pos)
                        plot(x,T','g', x,T_s,'b',x,T_l,'r');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanys ' yrs']) ; xline([w -w],'-.'); xlim(int_xlim); ylim(int_ylim);
                        hold on
                    end
    
                    % Plotting separately
                    if int_Ind == true
                        figure() 
                        plot(x,T','g',x,T_s,'b',x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanys ' yrs']) ; xline([w -w],'-.',{'Right Wall','Left Wall'});  xlim(int_xlim); ylim(int_ylim);
                        legend ('Temperature profile', 'T solidus', 'T liquidus','Location','eastoutside');
                    end

                else
                    % Adding plot to multiplot
                    if int_Subpl == true
                        figure(1)
                        nexttile(pos)
                        plot(x,T','g', x,Tr','c', x,T_s,'b',x,T_l,'r');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanys ' yrs (R = ' num2str(R) ')']) ; xline([w -w],'-.');  xlim(int_xlim); ylim(int_ylim);
                        hold on
                    end
    
                    % Plotting separately
                    if int_Ind == true
                        figure() 
                        plot(x,T','g', x,Tr','c', x,T_s,'b',x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanys ' yrs (R = ' num2str(R) ')']) ; xline([w -w],'-.',{'Right Wall','Left Wall'});  xlim(int_xlim); ylim(int_ylim);
                        legend ('Original dyke reinjection', 'Mobile melt reinjection', 'T solidus', 'T liquidus','Location','eastoutside');
                    end
                end 
                
                  if R<1
                        figure(2)
                        %nexttile(1)
                        plot(x,T');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Temperature evolution with time') ;  xlim(ini_xlim); ylim(ini_ylim);
                        hold on
                        for k = 1:Nx
                            res(k) = 1/(res0*exp(-E/(bltzm*(T(k)+273.15))));
                        end
                        figure(3)
                        %nexttile(2)
                        plot(x,res);xlabel ('Distance (m)'); ylabel ('Resistivity') ; yscale log; title('Resistivity evolution with time') ; xlim(ini_xlim); axis ([-fw*w fw*w 0.1 100]);
                        hold on
                   else
                        figure(2)
                        %nexttile(1)
                        plot(x,Tr);xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Temperature evolution with time') ;  xlim(ini_xlim); ylim(ini_ylim);
                        hold on
                        for k = 1:Nx
                            res(k) = 1/(res0*exp(-E/(bltzm*(T(k)+273.15))));
                        end
                        figure(3)
                        %nexttile(2)
                        plot(x,res);xlabel ('Distance (m)'); ylabel ('Resistivity'); title('Resistivity evolution with time') ;  xlim(ini_xlim); axis ([-fw*w fw*w 0.1 100]);
                        hold on
                   end            
                
                c_t = c_t+1;
                pos = pos+1;
                
            end
        end
    

        % Error value and convergence
        err = max(abs(T-Tlast)); %Find difference between last two solutions
    
        if err < err_lim 
        break; % Stop if solutions very similar, we have convergence
        end
    
    end
    
    
    if iter==maxiter
        %warning('Convergence not reached')
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Steady state
    
    tempstotalfinal= iter*dt/3600;
    tempstotalfinal; 
    anysfinal=tempstotalfinal/(24*365);
    sanysfinal=num2str(anysfinal);
    leyendfinal=['Steady state t=' sanysfinal 'years'];
    
    % Plotting the steady state
    if R<1
        if est_Analytic == true
            ASsteady =  Analytic_sol(tempstotalfinal*3600,x,Nx,To,Tc,w,K);
            
            if est_Subpl == true
                figure(1)
                nexttile(pos)
                plot(x,T,'g', x,T_s,'b' ,x,T_l,'r', x, ASsteady,'m');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanysfinal ' yrs (steady state/reinjection time)']) ; xline([w -w],'-.');legend ('Intrusion temperature profile', 'T solidus', 'T liquidus','Analytic solution','Location','eastoutside'); xlim(est_xlim); ylim(est_ylim);
                hold on
            end

            if est_Ind == true
                figure()
                nexttile(pos)
                plot(x,T,'g', x,T_s,'b' ,x,T_l,'r', x, ASsteady,'m', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanysfinal ' yrs (steady state/reinjection time)']) ; xline([w -w],'-.',{'Right Wall','Left Wall'});legend ('Intrusion temperature profile', 'T solidus', 'T liquidus','Analytic solution','Location','eastoutside'); xlim(est_xlim); ylim(est_ylim);
            end
        
        else
            if est_Subpl == true
                figure(1)
                nexttile(pos)
                plot(x,T,'g', x,T_s,'b' ,x,T_l,'r');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanysfinal ' yrs (steady state/reinjection time)']) ; xline([w -w],'-.');legend ('Intrusion temperature profile', 'T solidus', 'T liquidus','Location','eastoutside'); xlim(est_xlim); ylim(est_ylim);
                hold on
            end

            if est_Ind == true
                figure()
                nexttile(pos)
                plot(x,T,'g', x,T_s,'b' ,x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanysfinal ' yrs (steady state/reinjection time)']) ; xline([w -w],'-.',{'Right Wall','Left Wall'});legend ('Intrusion temperature profile', 'T solidus', 'T liquidus','Location','eastoutside'); xlim(est_xlim); ylim(est_ylim);
            end
        end

    else
        if est_Subpl == true
            figure(1)
            nexttile(pos)
            plot(x,T,'g', x,Tr','c', x,T_s,'b' ,x,T_l,'r');xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanysfinal ' yrs (steady state/reinjection time)(R = ' num2str(R) ')']) ; xline([w -w],'-.'); legend ('Original dyke reinjection', 'Mobile melt reinjection', 'T solidus', 'T liquidus','Location','eastoutside'); xlim(est_xlim); ylim(est_ylim);
            hold on
        end

        if est_Ind == true
            figure()
            nexttile(pos)
            plot(x,T,'g', x,Tr','c', x,T_s,'b' ,x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title(['Temp. at t = ' sanysfinal ' yrs (steady state/reinjection time)(R = ' num2str(R) ')']) ; xline([w -w],'-.',{'Right Wall','Left Wall'});legend ('Original dyke reinjection', 'Mobile melt reinjection', 'T solidus', 'T liquidus','Location','eastoutside'); xlim(est_xlim); ylim(est_ylim);
        end
    end

   if R<1
        figure(2)
        %nexttile(1)
        plot(x,T', x,T_s,'b',x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Temperature evolution with time') ; xline([w -w],'-.',{'Right Wall','Left Wall'}); xlim(ini_xlim); ylim(ini_ylim);legend ('Initial profile','1 month', '2 months', '3 months','4 months','5 months', '6 months','7 months','8 months', '9 months', '10 months','11 months','12 months','Solidus Temp.','Liquidus Temp.', 'MME Temp.');
        hold on
        for k = 1:Nx
            res(k) = 1/(res0*exp(-E/(bltzm*(T(k)+273.15))));
        end
        figure(3)
        %nexttile(2)
        plot(x,res);xlabel ('Distance (m)'); ylabel ('Resistivity (\Omega\cdot m)') ; yscale log; title('Resistivity evolution with time') ; xline([w -w],'-.',{'Right Wall','Left Wall'}); xlim(ini_xlim); axis ([-fw*w fw*w 0.1 100]);legend ('Initial profile','1 month', '2 months', '3 months','4 months','5 months', '6 months','7 months','8 months', '9 months', '10 months','11 months','12 months');
        hold on
   else
        figure(2)
        %nexttile(1)
        plot(x,Tr, x,T_s,'b', x,T_l,'r', x, T_m,"-.k");xlabel ('Distance (m)'); ylabel ('T (ºC)') ; axis ([-fw*w fw*w To-10 Tc+10]); title('Temperature evolution with time') ; xline([w -w],'-.',{'Right Wall','Left Wall'}); xlim(ini_xlim); ylim(ini_ylim);legend ('Initial profile','1 month', '2 months', '3 months','4 months','5 months', '6 months','7 months','8 months', '9 months', '10 months','11 months','12 months','Solidus Temp.','Liquidus Temp.', 'MME Temp.');
        hold on
        for k = 1:Nx
            res(k) = 1/(res0*exp(-E/(bltzm*(T(k)+273.15))));
        end
        figure(3)
        %nexttile(2)
        plot(x,res);xlabel ('Distance (m)'); ylabel ('Resistivity (\Omega\cdot m)'); title('Resistivity evolution with time') ; xline([w -w],'-.',{'Right Wall','Left Wall'}); xlim(ini_xlim); axis ([-fw*w fw*w 0.1 100]);legend ('Initial profile','1 month', '2 months', '3 months','4 months','5 months', '6 months','7 months','8 months', '9 months', '10 months','11 months','12 months');
        hold on
   end    

    pos = pos+1;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCIONS

% Solució Analítica
function T_a = Analytic_sol(t_b,x,Nx,To,Tc,w,K)
    x=x';
    for i=1:Nx
        T_a(i)=To+((Tc-To)/2)*(erf((x(i)+w)/(2*(K*(t_b))^(1/2)))-erf((x(i)-w)/(2*(K*(t_b))^(1/2)))); 
    end
end
