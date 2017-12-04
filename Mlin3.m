function Mlin = Mlin3(TimeSeries,L1,L2,F,n,dt,jump)

Mlin = eye(n);

for i=1:jump
    Temp = RK3deriv(TimeSeries(:,i),L1,L2,F,n,dt);
    Mlin = (Temp)*Mlin;
end

end

