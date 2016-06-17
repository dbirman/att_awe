
function out = conModel(con,params)

out = params.Rmax .* ((con.^params.n) ./ (con.^params.n + params.c50.^params.n)); 
