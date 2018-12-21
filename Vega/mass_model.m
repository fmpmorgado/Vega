function mass = mass_model(Isp, delta_v, Thrust)

e = [0.0001, 0.0001, 0.0001];
%e = [0.088, 0.094, 0.11];
m_pl = 1805;
g = 9.81;
%k=[96243/8533, 26300/2486, 12000/1433];
%k=[11.27, 10.5, 8.3];
k=[0, 0, 0];
m_struct = [0, 0, 0];
m_prop = [0, 0, 0];
cont = 1;
delta_v_stage = [0.25*delta_v, 0.25*delta_v, 0.5*delta_v];

for n= 3:-1:1
    a = 0;
    %Ve = Isp(1,n)*g
    k(1,n)=exp(delta_v_stage(1,n)/(g*Isp(1,n)));
    %k(1,n) = exp(delta_v(1,n)/Ve);
    %en = e(1,n);
    while a == 0 
        
        %k(1,n)=exp(delta_v/(g*Isp(1,n)));
        m_prop(1,n) = m_pl*(k(1,n)-1)*(1-e(1,n))/(1-(k(1,n)*e(1,n)));
        m_struct(1,n) = (e(1,n)*(k(1,n)-1))/((1-(e(1,n)*k(1,n))))*m_pl;
        
        const = m_struct(1,n)/(m_struct(1,n)+m_prop(1,n));
        b = 1.2*(m_struct(1,n)+m_prop(1,n))*9.81;
        
        if const >= 0.99*e(1,n) && const <= 1.01*e(1,n) && Thrust(1,n) >= b
            a = 1;
            %fprintf ('\n %4.2f',n)
            %fprintf ('\n %4.2f',cont)
            %fprintf ('\n %4.2f',Thrust(1,n))
            %fprintf ('\n %4.2f',b)
            
        else
            e(1,n) = const;
            cont = cont+1;
            
        end
    end
    m_pl= m_pl+m_struct(1,n)+m_prop(1,n);

end

mass(:,1) = m_struct;
mass(:,2) = m_prop;
%fprintf('\n %4.2f',e)
%fprintf('\n %4.2f',k)
end