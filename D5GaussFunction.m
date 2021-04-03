function F = D5GaussFunction(x,xdata)
 F = x(1)*exp(   -((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) )    )+...
     x(6)*exp(   -((xdata(:,:,1)-x(7)).^2/(2*x(8)^2) + (xdata(:,:,2)-x(9)).^2/(2*x(10)^2) )    )+...
     x(11)*exp(   -((xdata(:,:,1)-x(12)).^2/(2*x(13)^2) + (xdata(:,:,2)-x(14)).^2/(2*x(15)^2) )    )+...
     x(16)*exp(   -((xdata(:,:,1)-x(17)).^2/(2*x(18)^2) + (xdata(:,:,2)-x(19)).^2/(2*x(20)^2) )    );%+