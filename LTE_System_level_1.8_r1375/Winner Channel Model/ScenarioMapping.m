function output = ScenarioMapping(input)
%SCENARIOMAPPING Scenario is given as numerical vector with one element representing 
%   each link. This function maps the element to the corresponding letters as 
%   1=A1, 2=A2, 3=B1, 4=B2, 5=B3, 6=B4, 7=B5a, 8=B5c, 9=B5f, 10=C1, 
%   11=C2, 12=C3, 13=C4, 14=D1 and 15=D2a 

%   Authors: Mikko Alatossava (CWC/UOULU)
%   Created: 12.2.2007

switch (input)
    case {1}
        output = 'A1';
    case {2}
        output = 'A2';
    case {3}
        output = 'B1';
    case {4}
        output = 'B2';
    case {5}
        output = 'B3';
    case {6}
        output = 'B4';
    case {7}
        output = 'B5a';
    case {8}
        output = 'B5c';
    case {9}
        output = 'B5f';
    case {10}
        output = 'C1';
    case {11}
        output = 'C2';
    case {12}
        output = 'C3';
    case {13}
        output = 'C4';
    case {14}
        output = 'D1';
    case {15}
        output = 'D2a';
end