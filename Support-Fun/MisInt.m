function z = MisInt(w_mis,w_tag)
Nite = length(w_mis);
z = zeros(Nite,1);
for ite = 1:Nite
    z(ite) = logmean(w_mis{ite}+w_tag{ite},"all");
end
end

