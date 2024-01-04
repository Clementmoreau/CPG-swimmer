function psi = psi_custom(t,t_switch,switch_width)

test = t_switch*(t/t_switch - floor(t/t_switch));

% test whether we are in the linear switch part
if (test < switch_width/2) || (test > (t_switch - switch_width/2))
    nb_sw = floor((t-switch_width)/t_switch);
    if mod(nb_sw,2) == 1
        psi = -pi*(t-(nb_sw+1)*t_switch)/switch_width + pi/2;
    else
        psi = pi*(t-(nb_sw+1)*t_switch)/switch_width + pi/2;
    end
else
    psi = pi*(mod(floor(t/t_switch),2));
end

end