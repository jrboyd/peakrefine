my_writeFASTQ.dt = function(read, quality, name, file, append = FALSE, block_size = 1e5){
    require(data.table)
    stopifnot(all(
        length(read) ==
            c(
                length(quality),
                length(name)
            )))

    for(i in seq_len(ceiling(length(read) / block_size))){
        s = (i - 1) * block_size + 1
        e = (i) * block_size
        e = min(e, length(read))
        dt = data.table(r2 = read[seq(s, e)], r4 = quality[seq(s, e)], r1 = name[seq(s, e)])
        dt[, r1 := paste0("@", r1)]
        dt[, r3 := sub("@", "+", r1)]
        dt[, id := seq(.N)]
        dt = dt[, .(r1, r2, r3, r4, id)]
        dt = melt(dt, id.vars = "id")
        dt = dt[order(variable)][order(id)]


        if(i == 1){
            fwrite(dt[, .(value)], file = file, append = append, col.names = FALSE)
        }else{
            fwrite(dt[, .(value)], file = file, append = TRUE, col.names = FALSE)
        }
    }
}

my_writeFASTQ = function(read, quality, name, file, append = FALSE, block_size = 1e5){
    stopifnot(all(
        length(read) ==
            c(
                length(quality),
                length(name)
            )))

    for(i in seq_len(ceiling(length(read) / block_size))){
        s = (i - 1) * block_size + 1
        e = (i) * block_size
        e = min(e, length(read))
        df = data.frame(r2 = read[seq(s, e)], r4 = quality[seq(s, e)], r1 = name[seq(s, e)])
        df$r1 = paste0("@", df$r1)
        df$r3 = sub("@", "+", df$r1)
        df$id = seq(nrow(df))
        df = df[, c("r1", "r2", "r3", "r4", "id")]
        df = reshape::melt(df, id.vars = "id")
        df = df[order(df$variable), ][order(df$id), ]

        if(i == 1){
            write.table(df[, "value", drop = F], file = file, append = append, col.names = FALSE, row.names = FALSE, quote = FALSE)
        }else{
            write.table(df[, "value", drop = F], file = file, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
        }
    }
}

make_sim = function(genome, n_reads, f_enrich, seed = 1, bind_p = .05, no_binding = FALSE){
    set.seed(seed)
    transition <- list(Binding=c(Background=1), Background=c(Binding= bind_p, Background= 1 - bind_p))
    transition <- lapply(transition, "class<-", "StateDistribution")

    init <- c(Binding=0, Background=1)
    class(init) <- "StateDistribution"

    backgroundFeature <<- function(start, length=1000, shape=1, scale=20){
        weight <- rgamma(1, shape=1, scale=20)
        params <- list(start = start, length = length, weight = weight)
        class(params) <- c("Background", "SimulatedFeature")

        params
    }

    bindingFeature <<- function(start, length=500, shape=1, scale=20, enrichment=f_enrich, r=1.5){
        stopifnot(r > 1)

        avgWeight <- shape * scale * enrichment
        lowerBound <- ((r - 1) * avgWeight)
        weight <- actuar::rpareto1(1, r, lowerBound)

        params <- list(start = start, length = length, weight = weight)
        class(params) <- c("Binding", "SimulatedFeature")

        params
    }

    generator <<- list(Binding=bindingFeature, Background=backgroundFeature)

    constRegion <<- function(weight, length) rep(weight, length)
    featureDensity.Background <<- function(feature, ...) constRegion(feature$weight, feature$length)
    featureDensity.Binding <<- function(feature, ...){
        featDens <- numeric(feature$length)
        featDens[floor(feature$length/2)] <- feature$weight
        featDens
    }

    fragLength <<- function(x, minLength, maxLength, meanLength, ...){
        sd <- (maxLength - minLength)/4
        prob <- dnorm(minLength:maxLength, mean = meanLength, sd = sd)
        prob <- prob/sum(prob)
        prob[x - minLength + 1]
    }

    randomQuality <<- function(read, ...){
        paste(sample(unlist(strsplit(rawToChar(as.raw(73)),"")),
                     nchar(read), replace = TRUE), collapse="")
    }
    dfReads <<- function(readPos, readNames, sequence, readLen, ...){
        ## create vector to hold read sequences and qualities
        readSeq <- character(sum(sapply(readPos, sapply, length)))
        readQual <- character(sum(sapply(readPos, sapply, length)))

        idx <- 1
        ## process read positions for each chromosome and strand
        for(k in length(readPos)){ ## chromosome
            for(i in 1:2){ ## strand
                for(j in 1:length(readPos[[k]][[i]])){
                    ## get (true) sequence
                    readSeq[idx] <- as.character(ChIPsim::readSequence(readPos[[k]][[i]][j], sequence[[k]],
                                                                       strand=ifelse(i==1, 1, -1), readLen=50))
                    ## get quality
                    readQual[idx] <- randomQuality(readSeq[idx])
                    ## introduce sequencing errors
                    readSeq[idx] <- ChIPsim::readError(readSeq[idx], ChIPsim::decodeQuality(readQual[idx]))
                    idx <- idx + 1
                }
            }
        }
        data.frame(name=unlist(readNames), sequence=readSeq, quality=readQual,
                   stringsAsFactors = FALSE)
    }

    myFunctions <- ChIPsim::defaultFunctions()
    myFunctions$readSequence <- dfReads

    # featureArgs <<- list(generator=generator, transition=transition, init=init, start = 0,
    #                      length = 1e6, globals=list(shape=1, scale=20), experimentType="TFExperiment",
    #                      lastFeat=c(Binding = FALSE, Background = TRUE), control=list(Binding=list(length=50)))
    readDensArgs <<- list(fragment=fragLength, bind = 50, minLength = 150, maxLength = 250,
                          meanLength = 200)

    features <- ChIPsim::placeFeatures(generator, transition, init, start = 0, length = 1e6, globals=list(shape=1, scale=20),
                                       experimentType="TFExperiment", lastFeat=c(Binding = FALSE, Background = TRUE),
                                       control=list(Binding=list(length=50)))
    # if(no_binding){#remove binding
    #     k = sapply(features, function(x)class(x)[1]) == "Background"
    #     features = features[k]
    #         }

    myFunctions$features = function(...)features

    # dens <- ChIPsim::feat2dens(features)
    # readDens <- ChIPsim::bindDens2readDens(dens, fragLength, bind = 50, minLength = 150, maxLength = 250,
    #                                        meanLength = 200)


    set.seed(seed)
    simulated <- ChIPsim::simChIP(n_reads, genome, file = "", functions = myFunctions,
                                  control = ChIPsim::defaultControl(readDensity=readDensArgs))

    print(table(sapply(simulated$features[[1]], function(x)class(x)[1])))

    simulated
}

