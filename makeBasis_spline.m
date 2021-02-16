function [B] = makeBasis_spline(num,T)

B = spcol(linspace(0,1,num+3),3,linspace(0,1, T)) ;

end

