.calculateFingerprintMap <- function(data, maxSQSigma = 5000, epsilon = 0.0000001) {

    fingerprintMap = matrix(0, maxSQSigma+1, length(rawData(data)$reads))
    scSp = t(assay(scaleSpace(data), 1))

    # gauss with different sigma^2
    for (j in 1:(maxSQSigma+1)) {

        base = scSp[j,]
        sode = rep(0, length(base))
        # calculate 2nd order difference equations
        for (k in 2:(length(sode)-1)) {
            sode[k] = base[k+1] - 2*base[k] + base[k-1]
        }
        # inflection points (ip)
        ip = rep(0, length(base))
        for (k in 2:(length(sode)-1)) {
            if (sode[k+1] > epsilon & sode[k-1] < - epsilon) {
                ip[k] = 2  # to better differentiate "-" from "+" in grepl: don't use "-1" and "1"
            } else if (sode[k+1] < - epsilon & sode[k-1] > epsilon) {
            ip[k] = -1
            }
        }
        fingerprintMap[j,] = ip
    }

    assay(scaleSpace(data), 2) = t(fingerprintMap) 

    return(data)
}



setMethod("calculateFingerprintMap",
    signature=signature(data="Scale4C"),
    .calculateFingerprintMap)

