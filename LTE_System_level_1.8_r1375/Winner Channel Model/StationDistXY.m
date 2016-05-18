% [dist]=StationDistXY(S1, S2 [,flag])
% calculates the distance in the XY-plane between stations S1 and S2
% S1, S2 can be vectors of stations
% if flag is omitted or 0 S1/S2 must be of same size and dist contains the
% distance between S1(i) and S2(i)
% if flag is 1 S1/S2 can be of any size and dist is a length(S1)xlength(S2)
% matrix, where dist(i,j) is the distance between S1(i) ans S2(j)
%
% Author: Martin Käske (TUI)
function [dist]=StationDistXY(S1, S2, varargin)
NS1=length(S1); NS2=length(S2);

flag=0;
if((nargin > 2) && varargin{1}==1)
    flag=1;
end;

PosS1=[S1.Pos]; PosS2=[S2.Pos];
if(flag==1)
    dist=sqrt((repmat(PosS1(1,:)',1,NS2)-repmat(PosS2(1,:),NS1,1)).^2+...
              (repmat(PosS1(2,:)',1,NS2)-repmat(PosS2(2,:),NS1,1)).^2);
else
    if(NS1~=NS2)
        error('S1 and S2 must be of same size');
    end;
    dist=sqrt((PosS1(1,:)-PosS2(1,:)).^2+(PosS1(2,:)-PosS2(2,:)).^2);
end;
