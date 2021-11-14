function out = d_w1_cos_w2(w, wh)

    if all(w <= wh) || all(wh <= w)
        out = d_w1_w2( [w(1); d_cos(w(2), wh(2))] , ...
                       [wh(1); d_cos(wh(2), w(2))]);
    else
        error('Error: d_w1_cos_w2 requires ordeered inputs');
    end
end