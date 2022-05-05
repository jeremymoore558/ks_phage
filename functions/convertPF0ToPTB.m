function [TB,PTB] = convertPF0ToPTB(F0,PF0)
% deal with nans

TB = 1./(1+exp(F0));
dF0dTB = -1./(TB.*(1-TB));
PF0 = PF0/trapz(F0,PF0);

PTB = PF0.*abs(dF0dTB);
PTB = PTB/abs(trapz(TB,PTB));

end

