function wire = ProcessLoops(wire)

    loop_start  = wire.loop_start(:).';  
    loop_end    = wire.loop_end(:).';
    loop_length = 0;
    for i = 1:numel(loop_start)

        s = loop_start(i);
        e = loop_end(i);

        c_s = s-loop_length;
        c_e = e-loop_length;
        
        first_point = wire.F_point(s:e,:);
        second_point = wire.S_point(s:e,:);
        loop_length = loop_length+length(wire.F_point(s:e,:));

        if norm(first_point(c_s,:)-second_point(c_e,:)) == 0
            wire.F_point(s:e,:) = [first_point(c_e,:); first_point(c_s:c_e-1,:)];
            wire.S_point(s:e,:) = first_point(c_s:c_e,:);
            wire.T_point(s:e,:) = second_point(c_s:c_e,:);
        else
            wire.F_point(s:e,:) = first_point(c_s:c_e-2,:);
            wire.S_point(s:e,:) = first_point(c_s+1:c_e-1,:);
            wire.T_point(s:e,:) = first_point(c_s+2:c_e,:);
        end

    end
end
