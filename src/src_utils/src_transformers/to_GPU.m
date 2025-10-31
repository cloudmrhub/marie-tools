function[T_gpu] = to_GPU(T,gpu_flag)
    % GPU transfer file
    if isa(T, 'function_handle') || isempty(T)
        T_gpu = T;
        return;
    end

    switch gpu_flag
        case 0
            if isstruct(T)
                % Transfer Tucker compressed array to GPU
                fn = fieldnames(T);
                T_gpu  = repmat(struct(), size(T));
                for i = 1:numel(T)
                    for k = 1:numel(fn)
                        T_gpu(i).(fn{k}) = gpuArray(T(i).(fn{k}));
                    end
                end
            elseif issparse(T)
                % Transfer array to GPU
                T_gpu = gpuArray(T);
            else
                % Transfer array to GPU
                T_gpu = gpuArray(T);
            end
        case 1
            if isstruct(T)
                % Transfer Tucker compressed array to GPU
                fn = fieldnames(T);
                T_gpu  = repmat(struct(), size(T));
                for i = 1:numel(T)
                    for k = 1:numel(fn)
                        T_gpu(i).(fn{k}) = gpuArray(T(i).(fn{k}));
                    end
                end
            elseif issparse(T)
                % Transfer array to GPU
                T_gpu = gpuArray(T);
            else
                % Transfer array to GPU
                gpu_mem = gpuDevice().TotalMemory;
                array_mem = whos('T').bytes;
                if array_mem < 1/10*gpu_mem
                    T_gpu = gpuArray(T);
                else
                    T_gpu = gather(T);
                end
            end
        otherwise
            T_gpu = gather(T);
    end
end