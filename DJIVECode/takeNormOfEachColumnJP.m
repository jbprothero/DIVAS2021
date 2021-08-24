function outmat = takeNormOfEachColumnJP(inmat)

    [d,~] = size(inmat) ;
    sumsqs = sqrt(sum(inmat.^2)) ;
    outmat = inmat ./ repmat(sumsqs, d, 1) ;

end
