function DspRst (ite,N,Nc_Ns,N_cal)
if ite==1
    fprintf( '---------------------------------------------\n') ;
    fprintf( '     i      Nsam      Nc      Ns      Ncal   \n') ;
    fprintf( '---------------------------------------------\n') ;
end
fprintf( '%6d    %6d  %6d  %6d    %6d \n',ite,N,Nc_Ns(1),Nc_Ns(2),N_cal);
return
