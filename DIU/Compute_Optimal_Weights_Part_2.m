clc;
d = 2.05;  % element spacing is 2.05"
n = 1:1:8;
c = 300e6; 	% speed of light in m/s
d = d*2.54/100; % convert d to meters
i = 1;
tau = zeros(1,8);
Ideal_tau = zeros(181, 8);  %initialize array of ideal TTDs
Hardware_tau = zeros(181, 8); %initialize array of achievable TTDs (for error analysis)

Error_Low = 0;
Error_High = 0;

Max_Error_Low_Index = zeros(1,181); %initialize array
Max_Error_High_Index = zeros(1,181); %initialize array
Error_Final = zeros(181,8);
Index_Angle = 1; 

for angle_degrees = -90:1:90

    % determine ideal (perfect) time delays
    tau = round(((n-1)*(d/c)*sin(angle_degrees*pi/180))*1e12) ; % ideal element time delays in ps
    Min_tau = min(tau); %find the minimum value of all ideal delays
    tau = tau - Min_tau; %Offset tau's by this value to shift all delays to positive values
                         %either the lowest or highest tau will then be 0,
                         %depending on the AoA.
    if(angle_degrees == -5)
        i = i;
    end
    % added array to store all ideal tau values for all angles
    Ideal_tau(Index_Angle, :) = tau;
    
    Max_Error_Low = 1000000;  %initialize to a high number
    Max_Error_High = 1000000;
    for VPD_count = 1:1:VPD_Total
        if VPD_count == 85
                i = i;
        end
        Error_Low = tau(1:4) - path_delay( Valid_Path_Delays(VPD_count),1:4);
        Error_Low = Error_Low.^2;
        Error_Low_Sum_of_Squares = sum(Error_Low);
        if (Error_Low_Sum_of_Squares < Max_Error_Low)
            Max_Error_Low = Error_Low_Sum_of_Squares;
            Max_Error_Low_Index(Index_Angle) = VPD_count;
%           Hardware_tau(Index_Angle, 1:4) = path_delay( Valid_Path_Delays(VPD_count),1:4);
            Best_Low_VPD = Valid_Path_Delays(VPD_count); %Keep track of the lowest error for later
        end
         
        Error_High = tau(5:8) - path_delay(Valid_Path_Delays(VPD_count), 1:4);
        Error_High = Error_High.^2;
        Error_High_Sum_of_Squares = sum(Error_High);
        if (Error_High_Sum_of_Squares < Max_Error_High)
            Max_Error_High = Error_High_Sum_of_Squares;
            Max_Error_High_Index(Index_Angle) = VPD_count; 
%            Hardware_tau(Index_Angle, 5:8) = path_delay( Valid_Path_Delays(VPD_count),1:4);
            Best_High_VPD = Valid_Path_Delays(VPD_count); %Keep track of the lowest error for later
        end
        i = i;
    end
    Hardware_tau(Index_Angle, 1:4) = path_delay( Best_Low_VPD,1:4);
    Hardware_tau(Index_Angle, 5:8) = path_delay( Best_High_VPD,1:4);
    Error_Final(Index_Angle,:) = Ideal_tau(Index_Angle,:) - Hardware_tau(Index_Angle,:);
    Index_Angle = Index_Angle + 1;
    display(angle_degrees)
end

%Now we want to create the final output array of TTD Switch settings 
% at the very end we will invert all teh bits to be consistent with the
% hardware control convention.
Output_Table = zeros(181, 37); %The extra column is for the angle index
Index_Angle = 1;
for angle_degrees = -90:1:90
    %display(angle_degrees)
    Output_Table(Index_Angle, 1) = angle_degrees;
    % Map the 4 delay path results (18 TTDs) into the final 8 paths (36
    % TTDs for the Low and High portions of the BFN.
  
    for ii = 1:1:12
        Output_Table(Index_Angle, ii+1) = All_Numbers(Valid_Path_Delays(Max_Error_Low_Index(Index_Angle)),ii);
        Output_Table(Index_Angle, ii+12+1) = All_Numbers(Valid_Path_Delays(Max_Error_High_Index(Index_Angle)),ii);
    end
    for ii = 13:1:16
        Output_Table(Index_Angle, ii+12+1) = All_Numbers(Valid_Path_Delays(Max_Error_Low_Index(Index_Angle)),ii);
        Output_Table(Index_Angle, ii+16+1) = All_Numbers(Valid_Path_Delays(Max_Error_High_Index(Index_Angle)),ii);
    end
    for ii = 17:1:18
        Output_Table(Index_Angle, ii+16+1) = All_Numbers(Valid_Path_Delays(Max_Error_Low_Index(Index_Angle)),ii);
        Output_Table(Index_Angle, ii+18+1) = All_Numbers(Valid_Path_Delays(Max_Error_High_Index(Index_Angle)),ii);
    end 
    Index_Angle = Index_Angle + 1;
end
% This is the bit flip loop (TTD active bit state is "0")
Index_Angle = 1;
for angle_degrees = -90:1:90
    for ii = 2:1:37
        if Output_Table(Index_Angle, ii) == 1
            Output_Table(Index_Angle, ii) = 0;
        else
            Output_Table(Index_Angle, ii) = 1;
        end
    end
    Index_Angle = Index_Angle + 1;
end
T = table(Output_Table);
% One way to write results to a text file...
writetable(T,'Output_Table.txt');

% ... and another way...
dlmwrite('TTD_data.txt',Output_Table);

%I think that's it!!!
        



