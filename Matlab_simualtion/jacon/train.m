function F = train(q,t0_3,l2,l3)
    F(1)=l2*cos(q(1))+l3*cos(q(1)+q(2))-sqrt(t0_3(1)^2 + t0_3(2)^2);
    F(2)=l2*sin(q(1))+l3*sin(q(1)+q(2))-t0_3(3);
end