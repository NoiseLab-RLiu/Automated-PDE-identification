function ux = FiniteDiff(u,dx,d)
u = squeeze(u);
if(size(u,2)~=1)
    u=u.'; % u should be a column vec
end
n = size(u,1);
ux = zeros(n,1); % column vector
if d==1
    for i=2:n-1
        ux(i) = (u(i+1)-u(i-1))/(2*dx);
    end
    ux(1) = (-1.5*u(1)+2*u(2)-0.5*u(3))/(dx);
    ux(n) = (1.5*u(n)-2*u(n-1)+0.5*u(n-2))/dx;
elseif d==2
    for i=2:n-1
        ux(i) = (u(i+1)-2*u(i)+u(i-1))/(dx^2);
    end
    ux(1) = (2*u(1)-5*u(2)+4*u(3)-u(4))/(dx^2);
    ux(n) = (2*u(n)-5*u(n-1)+4*u(n-2)-u(n-3))/(dx^2);
elseif d==3
    for i=3:n-2
        ux(i) = (0.5*u(i+2)-u(i+1)+u(i-1)-0.5*u(i-2))/(dx^3);
    end
    ux(1) = (-2.5*u(1)+9*u(2)-12*u(3)+7*u(4)-1.5*u(5))/(dx^3);
    ux(2) = (-2.5*u(2)+9*u(3)-12*u(4)+7*u(5)-1.5*u(6))/(dx^3);
    ux(n) = (2.5*u(n)-9*u(n-1)+12*u(n-2)-7*u(n-3)+1.5*u(n-4))/(dx^3);
    ux(n-1) = (2.5*u(n-1)-9*u(n-2)+12*u(n-3)-7*u(n-4)+1.5*u(n-5))/(dx^3);
else
    msg = 'Error occurred.';
    error(msg)
end
end
    
    

