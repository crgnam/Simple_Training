% converts date vector to a UNIX timestamp
function t = Time2UNIX(v)
    % Jeb's R2014b installation doesn't work for various reasons, and R2014a doesn't have these two functions, so this elegant line is useless
    %t = posixtime(datetime(v));
    
    % less-elegant method (doesn't account for leap seconds)
    t = int32(floor(86400*(datenum(v) - datenum('01-Jan-1970'))));
end
