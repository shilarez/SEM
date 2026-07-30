function [ W ] = WMatrix_New(NelX, NelY, NGLL, wgll2, dxe, dye)
for ey = 1:NelY
    for ex = 1:NelX
        eo = (ey-1)*NelX + ex;
        W(:,:,eo) = wgll2;
    end
end

end

