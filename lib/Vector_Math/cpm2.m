function [CrossProductMatrices] = cpm2(Vector)
% This is a dynamically sized version of the CPM.  It may be useful to find
% a vectorized version of this if possible?

    for ii = 1:3:size(Vector,1)
    CrossProductMatrices(ii:ii+2,1:3) = [      0      , -Vector(ii+2),  Vector(ii+1);
                                          Vector(ii+2),      0       , -Vector(ii);
                                         -Vector(ii+1),  Vector(ii)  ,       0      ];
    end           
end