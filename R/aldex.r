aldex <- function(reads, conditions, mc.samples=128, test="t",
    effect=TRUE, include.sample.summary=FALSE, verbose=FALSE, denom="all"){

    # wrapper function for the entire set of
    print("aldex.clr: generating Monte-Carlo instances and clr values")
    x <- aldex.clr(reads=reads, conds=conditions, mc.samples=mc.samples,
        denom=denom, verbose=verbose, useMC=FALSE)
    if(test == "t") {
        print("aldex.ttest: doing t-test")
        x.tt <- aldex.ttest(x, conditions, paired.test=FALSE)
    }
    if(test == "glm"){
        print("aldex.glm: doing Kruskal Wallace and glm test")
        x.tt <- aldex.glm(x, conditions)
    }
    if(effect == TRUE){
        print("aldex.effect: calculating effect sizes ")
        x.effect <- aldex.effect(x, conditions,
            include.sample.summary=include.sample.summary, verbose=verbose)
        z <- data.frame(x.effect, x.tt)
    } else {
        z <- data.frame(x.tt)
    }
    return(z)
}
