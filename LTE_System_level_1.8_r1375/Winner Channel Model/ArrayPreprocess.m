% a_pre=ArrayPreprocess(a)
% performs pre-processing of the array 'a'
% a can also be a vector of arrays
% rotate element field-pattern in ACS
% create one aperture for the whole array
%
% Author: Martin Käske (TUI)
function a_pre=ArrayPreprocess(a)
% there are a lot of similarities to AntennaResponse() so use this function
% preProcessArray is pretty much the same as ArrayGain with
% array.Rot=[0;0;0]
for i=1:length(a)
    a2=a(i);
    a2.Rot=[0;0;0];
    % even if no element is rotated, the following steps are performed to
    % create one aperture for the whole array
    NumElem=length(a2.Element);
    phi=linspace(-180,180-1/360,360); % 1 degree sampling should be dense enough
    g=AntennaResponse(a2,phi*pi/180);
    % leave out phase term
    BP=g{1}(:,1:2,:);
    sz_BP=size(BP);
    % reshape to make dimensions of BP compatible with BP2Aperture1D
    % i.e. add the (singleton) elevation dimension
    a2.Aperture=BP2Aperture1D(reshape(BP,[sz_BP(1:end-1) 1 sz_BP(end)]),phi);

    %set Rot of all elements to zero, remove Aperture
    for k=1:NumElem
        a2.Element(k).Rot=zeros(3,1);
    end;
    % no need to check if neither CommonAperture nor Element.Aperture is
    % present, because in this case it is not allowed to call
    % ArrayPreprocess (a2.Aperture should be there and the FP already in ACS)
    % i.e. it is okay for ArrayPreprocess to fail (actually the call to
    % AntennaResponse would fail)
    if(isfield(a2,'CommonAperture'))
        a2=rmfield(a2, 'CommonAperture');
    else
        a2.Element=rmfield(a2.Element,'Aperture');
    end;
    a2.Rot=a(i).Rot;
    a_pre(i)=a2;
end;
