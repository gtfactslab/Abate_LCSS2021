function out = d_w1_w2(w, wh)

    if all(w <= wh)
        out = min([w(1)*w(2), wh(1)*w(2), w(1)*wh(2), wh(1)*wh(2)]);
    elseif all(wh <= w)
        out = max([w(1)*w(2), wh(1)*w(2), w(1)*wh(2), wh(1)*wh(2)]);
    else
        error('Error: d_w1_w2 requires ordeered inputs');
    end
end