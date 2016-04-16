% 
% % % % % % % % %% % % % % % % % %% % % % % % % % %% % % % % % % % % % %% % % % % % %
%
% Author: D. Belmar
% Date: 7/7/15 
% This program calculates the True Time Delay table for the Beam Forming
% Network (BFN) on the K-3 program.  The original algorithm provided by the
% Electronics Supplier was not correct, and even after corrections were
% made turned out to be sub-optimal (John Albanese coded this version of 
% the algorithm. This was due to the non-uniform and non-monotonically 
% increasing distribution of the true time delay (TTD) values in the BFN.  
% There are two (maybe more) methods of deriving the BFN TTD values.  The 
% first involves a binary search algorithm, which has the 
% potential to be very fast, but becomes very complex due to the BFN
% architecture and the TTD values.  The second approach, and the one that
% is implemented here, is a brute force (BF) approach, where every possible
% combination of switch states is generated, and at every angle, the error
% is determined.  This is possible since, although there are a total of 36 
% TTD switches in the BFN (which would result in 2^36 possible combinations),
% the BFN is architected so that the upper and lower 4 delay paths are 
% independent (meaning they do not share any TTDs). as such only 2^18 
% combinations need to be generated. There are some efficiencies 
% incorporated in this BF algorithm where many of the possible switch 
% combinations that are invalid are discarded since they result in 
% non-monotonically increasing or decreasing time delays across the array).  
% This algorithm executes fairly quickly as a result, so there 
% is no need to implement the first method (except as an academic exercise).  
% 
% Time Delay Unit (TDU). A zero indicates the particular TDU is enabled
% errors - final time delay error from ideal in ps. an array of N values.
% Index one represents element 1
clc;
All_TTD_Switch_Settings = zeros(10,18);
j=1;

MAX_NUM = 262144;   % This is 2^18

All_Numbers = zeros(MAX_NUM ,18);
Delay_Path = zeros(MAX_NUM ,4);
HW_TTD = [30, 60, 120, 142, 282, 272, 542] ; 	% Hardware time delays in same order as path delay index
Path_Delay_Index = [3 2 1 14 13 18 17; 6 5 4 14 13 18 17; 9 8 7 16 15 18 17; 12 11 10 16 15 18 17];

% Part 1: Allocate and generate the array of all numbers from 0 to 262144 (18-bits)
% 
for j=1:18
    a(j) = 2^(18-j);
end

for N = 0 : MAX_NUM -1  % loop through all possible 18-bit numbers
    
    Number = N;
    for j = 18:-1:1
        %create an 18 column row vector that is a binary representation of N
        k = 18-j+1;
        if Number >= a(k)
            result = 1;
            Number = Number - a(k);
            All_Numbers(N+1,k) = result;
        else
            result = 0;
            All_Numbers(N+1,k) = result;
        end
    end
end
 %Part 2: For the low and high 4-group paths compute the total delay for
 %each number and put into a 262144x4 array
path_delay = zeros(MAX_NUM , 4);

for N = 0 : MAX_NUM -1  % loop through all possible 18-bit numbers
    % First do the lower 4 paths
    for Path_Number = 1 : 4
        path_delay(N+1, Path_Number) = 0;        
        for kk = 1 : 7
           path_delay(N+1, Path_Number) = path_delay(N+1, Path_Number) + All_Numbers(N+1, Path_Delay_Index(Path_Number, kk))*HW_TTD(kk);
        end
    end
end

% need to create a new matrix of valid Numbers from All_Numbers based on
% the path_delay results. Store in a vector that will only include indexes
% to valid numbers.
Num_of_Valid_Path_Delays = 1;  % the "All Zero" case is always valid.
VPD_Total = 1;
Valid_Path_Delays(VPD_Total) = 1;
% Now loop through all possible 18-bit numbers to eliminate all values that
% cannot be valid i.e, they are not monotonically increasing or
% monotonically decreasing, including the all equal case. This reduces the
% possible number of valid delays that need to be tested at each angle by
% roughly an order of magnitude   
for N = 2 : MAX_NUM  
    a1 = path_delay(N, 1);
    b1 = path_delay(N, 2);
    c1 = path_delay(N, 3);
    d1 = path_delay(N, 4);
    if ((a1 <= b1) && (b1 <= c1) && (c1 <= d1)) || ((a1 >= b1) && (b1 >= c1) && (c1 >= d1))
        VPD_Total = VPD_Total + 1;
        Valid_Path_Delays(VPD_Total) = N;
    else
        N=N;  
    end
end
% Go to Part 2
























    
    
    
    
    
%262144