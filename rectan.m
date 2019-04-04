function rect=rectan(t,Amp,T)


rect = zeros(size(t));
t_max = t(end);
t_min = t(1);
k = length(t);

zz = (t_max - t_min)/T;
% 
% for i=1:length(t)
%     if mod(t(i),0.5*T) == 0
%         rect(i) = 0;
%     elseif mod(t(i),T) < 5
%         rect(i) = 1;
%     elseif mod(t(i), T) > 5
%         rect(i) = -1;
%     end
% end

for j=1:t_max/T
    for i=1:T/t_max*k
        if i>1 && i<0.5*T/t_max*k + 1 
            rect((j-1)*k/t_max*T+i) = 1*Amp;
        elseif i==0.5*T/t_max*k + 1 || i==T/t_max*k + 1 || i == 1
            rect((j-1)*k/t_max*T+i) = 0;
        elseif i<T/t_max*k + 1 && i>0.5*T/t_max*k + 1
            rect((j-1)*k/t_max*T+i) = -1*Amp;
        end
    end
end
