%% dinklebach algorithm example from "Nonlinear Fractional Programming" appendix

clear
lambda=0;
while 1
    
  cvx_begin 
    variable x
    variable y
    maximize ( (-3*x^2-2*y^2+4*x+8*y-8)-lambda*(x^2+y^2-6*y+8) )
    subject to
         x>=0
         y>=0
         x+3*y<=5
    cvx_end  
    
    F_q= (-3*x^2-2*y^2+4*x+8*y-8)-lambda*(x^2+y^2-6*y+8);
    lambda=(-3*x^2-2*y^2+4*x+8*y-8)/(x^2+y^2-6*y+8);
    
    if abs(F_q)<0.001
        break
    end
end

