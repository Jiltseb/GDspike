function [ OUT ] = Inc_Spk1( IN , Theta)
    N=size(IN,1);
    %Peak as Spike--------------------------------------------
    OUT=zeros(N,1);
    for i=2:N-1
        if IN(i-1,1)<=IN(i,1)
           if IN(i+1,1)<=IN(i,1)
               if IN(i,1)>=Theta
                  OUT(i,1)=1; 
               end 
           end
        end
    end
end