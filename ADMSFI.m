function MFGAF = ADMSFI(data, delta)
    [n, m] = size(data);

    varepsilon = zeros(1, m);
    for j = 1:m
        if min(data(:, j)) == 0 && max(data(:, j)) == 1
            varepsilon(j) = std(data(:, j), 1) / delta; 
        end
    end

    En = zeros(m, m);
    r = cell(1, m);  % 使用 cell 数组存储每个 r1，而不是三维矩阵

    for j = 1:m
        r{j} = mfgad_rm(data(:, j), varepsilon(j)); 
    end
    
    % 计算 En，不再存储所有中间的 R
    for i = 1:m
        for j = 1:m
            R_sum = matrixOperation(r{i}, r{j});  % 直接计算每对矩阵的结果
            En(i, j) = -(1 / n) * sum(log2(sum(R_sum, 2) / n));
        end
    end

    % 继续处理其他步骤
    RE = zeros(n, m);
    for i = 1:m
        [~, En_de] = sort(En(i, :), 'descend');
        [~, En_as] = sort(En(i, :), 'ascend');
        
        weight1 = zeros(n, m);
        EN_de = zeros(n, m);
        EN_as = zeros(n, m);

        for k = 1:m
            FSet = mfgad_rm(data(:, k), varepsilon(k));  
            FSet_de = mfgad_rm(data(:, En_de(1:(m - k + 1))), varepsilon(En_de(1:(m - k + 1))));
            FSet_as = mfgad_rm(data(:, En_as(1:(m - k + 1))), varepsilon(En_as(1:(m - k + 1))));

            weight1(:, k) = sum(FSet, 2) / n;
            EN_de(:, k) = sum(FSet_de, 2);
            EN_as(:, k) = sum(FSet_as, 2);
        end

        ENS_de = (EN_de(:, 2:end) - EN_de(:, 1:end-1)) ./ EN_de(:, 2:end);
        ENS_as = (EN_as(:, 2:end) - EN_as(:, 1:end-1)) ./ EN_as(:, 2:end);
        
        MFGAF1 = 1 - (sum(weight1, 2) / m) .* sqrt((sum(ENS_de, 2) + sum(ENS_as, 2)) / (2 * m - 2));

        RE(:, i) = MFGAF1;
    end

    rowSums = sum(RE, 2) / m;
    MFGAF = rowSums;
end