function [hexColor] = color2Hex(color, sat, opac)
    
    if length(sat) == 1
        sat = ['0' sat];
    end

    if (color == 'y')
        color = ['00' sat sat];
    elseif (color == 'm')
        color = [sat '00' sat];
    elseif (color == 'c')
        color = [sat sat '00'];
    elseif (color == 'r')
        color = ['0000' sat];
    elseif (color == 'g')
        color = ['00' sat '00'];
    elseif (color == 'b')
        color = [sat '0000'];
    elseif (color == 'w')
        color = [sat sat sat];
    else
        color = '000000';
    end
    
    hexColor = [opac color];

end

