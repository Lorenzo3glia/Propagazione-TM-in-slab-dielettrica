function Htot = onde_piane_tm(x,z,kx,kzin_p,kzrif_p,kztra1,kztra3,Hinc_p,H1,Hrif_p,H3,d)
%Sto calcolando i campi magnetici
%% Costruisco Gli If

if z<-1*d
    Htot=H1*exp(-1*i*(kx*x+kztra1*z));
elseif (-1*d<z)&&(z<d)
    Htot=Hinc_p*exp(-1*i*((kx*x)+(kzin_p*z)))+Hrif_p*exp(-1*i*(kx*x + kzrif_p*z));
elseif z>=d
    Htot=H3*exp(-1*i*(kx*x+kztra3*z));
end