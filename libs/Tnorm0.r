mean.Tnorm <- function(mn, sd, a = 0, b = 1) {
    #si <- function(ep) (1/sqrt(2*pi)) * exp((-1/2) * ep^2)
    #zi <- function(x)  (1/2) * (1 + erf(x/sqrt(2)))
    sd = sd - mn
    alpha = (a-mn) / sd
    beta  = (b-mn) / sd

    forV <- function(mni, sdi, ai, bi){
        print("yay")
        Z = pnorm(bi, mni, sdi) - pnorm(ai, mni, sdi)
        Y = dnorm(ai, mni, sdi) - dnorm(bi, mni, sdi)
    
        mni + (Y/Z) * sdi
    }
    mask = !is.na(mn+sd)
    out = mn
    out[mask] = mapply(forV, mn[mask], sd[mask], alpha[mask], beta[mask])
    return(out)
}

