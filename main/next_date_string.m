function date_plus_one = next_date_string(date)
    date_breaked = date -'0';
    if (date_breaked(1:4)==[0,2,2,8])
        year = date_breaked(5)*10 + date_breaked(6);
        if mod(year, 4)>0
            date_breaked(1:4) = [0,3,0,1];
        else
            date_breaked(1:4) = [0,2, 2, 9];
        end
    elseif (date_breaked(1:4)==[0,2,2,9])
        date_breaked(2:4) = [0,3,0,1];
    elseif (date_breaked(3)==3)&&(date_breaked(4)==1)&&(sum(ismember(date_breaked(1:2),[[0,1];[0,3];[0,5];[0,7];[0,8];[1,0]]))>=2)
        month = date_breaked(1)*10 + date_breaked(2);
        month = month + 1;        
        date_breaked(1:2) = str2num_add_zero(month);
        date_breaked(3:4) = [0,1];
    elseif (date_breaked(1:4)==[1,2,3,1])
        date_breaked(1:4) = [0,1,0,1];       
        year = date_breaked(5)*10 + date_breaked(6);
        year = year + 1;   
        date_breaked(5:6) = str2num_add_zero(year);
    else
        day = date_breaked(3)*10 + date_breaked(4);
        day = day + 1;
        date_breaked(3:4) = str2num_add_zero(day);
        
    end
    date_plus_one = strcat(num2str(date_breaked));
    date_plus_one = date_plus_one(find(~isspace(date_plus_one)));
    
end

function str = str2num_add_zero(num)
if num<10
    str = [0, int2str(num)-'0'];
else
    str = int2str(num)-'0';
end
end