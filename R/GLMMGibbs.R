#this is the R code of version 0.5 of the GLMMGibbs R package,
# by Jonathan Myles (mylesj@icrf.icnet.uk) 

ordinal.end <- function(n){
  str <- c("st","nd","rd",rep("th",7))
  x <- n %% 10
  if (x == 0){
    ans <- "th"
  }
  else{
    ans <- str[x]
  }
  ans    
}

"glmm" <- 
function (formula, family, data, weights, offset, icm = 50, burnin = 1000, 
    keep = 1000, model.show = FALSE, progress.info = 0, store.results = FALSE, 
    thin = 1, bugsfile ,seed1 = 6762, seed2 = 92928, seed3 = 19729) 
{
    call <- match.call
    mm <- match.call(expand = F)
    aaa <- mm
    mm$family <- mm$icm <- mm$burnin <- mm$keep <- mm$model.show <- mm$bugsfile <- mm$store.results <- mm$progress.info <- mm$thin <- NULL
    mm[[1]] <- as.name("model.frame")
    mm <- eval(mm, sys.parent())
    model.init(mm, family = family, model.show = model.show)
    rtobj <- terms(formula)
    rtp <- rtparse(rtobj, data = data, length(model.extract(mm,response)) )
    for (i in 1:length(rtp)) {
        add.block(rtp[[i]], data = data)
    }
    if (model.show) {
        dummy <- .C("model_show", as.integer(0))
    }
    program.length <- ((keep > 0) + (burnin > 0) + (icm > 0)) * 
        2
    program <- numeric(0)
    if (icm > 0) {
        program <- c(0, icm)
    }
    if (burnin > 0) {
        program <- c(program, 1, burnin)
    }
    if (keep > 0) {
        program <- c(program, 1, keep)
    }
    number.of.effects <- numeric(1)
    number.of.hyperparameters <- numeric(1)
    aa <- .C("nparameters", as.integer(number.of.effects), as.integer(number.of.hyperparameters))
    number.of.effects <- aa[[1]]
    number.of.hyperparameters <- aa[[2]]
    effect.names <- character(number.of.effects)
    hyperparameter.names <- character(number.of.hyperparameters)
    bb <- .C("parameternames", as.character(effect.names), as.character(hyperparameter.names))
    sampled.effect.values <- matrix(0, nr = keep/thin, ncol = number.of.effects)
    storage.mode(sampled.effect.values) <- "double"
    sampled.hyperparameter.values <- matrix(0, nr = keep/thin, ncol = number.of.hyperparameters)
    storage.mode(sampled.hyperparameter.values) <- "double"
    x <- .C("onemodel_sample", as.integer(c(0, icm, 1, burnin, 
        1, keep)), as.integer(program.length), as.integer(seed1), 
        as.integer(seed2), as.integer(seed3), as.integer(progress.info), 
        as.integer(thin), sampled.effect.values, sampled.hyperparameter.values)
    
    res <- list(cbind(1:floor(keep/thin), x[[8]], x[[9]]))
    tmp <- 1/sqrt(x[[9]])
    attr(res, "class") <- "glmmfit"
    attr(res, "number.of.effects") <- number.of.effects
    attr(res, "call") <- aaa
    attr(res, "icm") <- icm
    attr(res, "burnin") <- burnin
    attr(res, "keep") <- keep
    attr(res, "thin") <- thin
    attr(res, "parameter.names") <- c(bb[[1]], bb[[2]])
    attr(res, "number.of.effects") <- number.of.effects
    attr(res, "number.of.hyperparameters") <- number.of.hyperparameters
    stdev <- function(x) {
        sqrt(var(x))
    }
    if (number.of.effects > 0) {
        statse <- cbind(apply(x[[8]], 2, mean), apply(x[[8]], 
            2, stdev), t(apply(x[[8]], 2, quantile, prob = c(0.025, 
            0.25, 0.5, 0.75, 0.975))))
    }
    if (number.of.hyperparameters > 0) {
        statsh <- cbind(apply(x[[9]], 2, mean), apply(x[[9]], 
            2, stdev), t(apply(x[[9]], 2, quantile, prob = c(0.025, 
            0.25, 0.5, 0.75, 0.975))))
        statssigma <- cbind(apply(tmp, 2, mean), apply(tmp, 2, 
            stdev), t(apply(tmp, 2, quantile, prob = c(0.025, 
            0.25, 0.5, 0.75, 0.975))))
    }
    if ((number.of.effects > 0) && (number.of.hyperparameters > 
        0)) {
        stats <- rbind(statse, statsh, statssigma)
    }
    if ((number.of.effects == 0) && (number.of.hyperparameters > 
        0)) {
        stats <- rbind(statsh, statssigma)
    }
    if ((number.of.effects > 0) && (number.of.hyperparameters == 
        0)) {
        stats <- rbind(statse)
    }
    attr(res, "stats") <- stats
    .C("free_glm", as.integer(0))
    if(!(missing(bugsfile))){
      bugs.out(res,bugsfile)
    }
    res
}
"rtparse" <-
function (rt, data, l) 
{
    ### rtparse  takes a terms object
    ### and parses it  into ``blocks'' A 
    ### block is the expression formed  by the sum of all terms including a given 
    ### set of random factors.  For example, if A and B are fixed effect factors, 
    ### R and S are random effect factors, and X and Z are covariates, consider
    ### the model
    ###
    ### 1 + A + B + R*X + Z + R*S
    ###
    ### which is, in full
    ###
    ### 1 + A + B + X + Z + R + S + R:X + R:S
    ###
    ### This needs to be parsed as
    ###
    ### (1 + A + B + X + Z) + R:(1 + X)  + S:(1)  +  R:S:(1)
    ###
    ###
    ### The returned object is a list of length equal to the number
    ### of blocks , each of which represents the
    ### result of a call to "rterms" with argument equal to the fised
    # firstly,we obatin a list of the variables, and
    # make a logical vector of the same length which will
    # hold whether or not each variable is a random factor.
    predictor.variables <- attr(rt, "variables")[-(1:2)]
    if (!(is.null(attr(rt, "offset")))) {
      predictor.variables <- predictor.variables[-(attr(rt, 
            "offset") - 1)]
    }
    number.of.predictor.variables <- length(predictor.variables)
    vars.of.interest <- logical(number.of.predictor.variables)
    # Determine whether each variable i in turn is a random factor, and
    # put the result in randomvars[i].  A factor is random iff it is of 
    # class "glmm.random.factor".  The following is quite ugly, but
    # deals  with the case where the data is in a frame.
    randomvars <- logical(number.of.predictor.variables)
    for (i in 1:number.of.predictor.variables) {
        if (missing(data)) {
            xx <- eval(predictor.variables[i])
        }
        else txt <- paste(as.character(substitute(data)), "$", 
            predictor.variables[i], sep = "")
        xx <- eval(parse(text = txt))
        xxr <- !(is.null(attr(xx, "is.random"))) && (attr(xx, 
            "is.random") == T)
        randomvars[i] <- xxr
        if (xxr) {
            xxi <- attr(xx, "of.interest")
        }
        else {
            xxi <- T
        }
        vars.of.interest[i] <- xxi
    }
    # we create a logical vector, done, with length equal to the number of columns
    # in the factors attribute of rt.  Each column represents a term of the
    # expression.  As each term is parsed into its correct block, the appropriate
    # ellement of done is set  to T. The object "res" wil eventually be the parsed
    # expression returned by the function.  For convenience,  put the "factors"
    # attribute into "mat".  Also  create a submatri "submat" which contains all
    # those rows representing random factors. The variable "l" will count the blocks 
    # as they are being parsed.
    done <- rep(F, ncol(attr(rt, "factors")))
    if  (  (!(is.null(attr(rt, "offset"))))   &&  (attr(rt,"offset"))){
    mat <- attr(rt, "factors")[-(c(attr(rt, "response"),attr(rt,"offset"))), , drop = F]      
    }
    else{
    mat <- attr(rt, "factors")[-(attr(rt, "response")), , drop = F]
  }
    submat <- mat[randomvars, , drop = F]
    random.of.interest <- vars.of.interest[randomvars]
    l <- 1
    res <- list()
   
    #
    # The easy case first.  If "submat" is numeric(0), this means that
    # there are no random factors.  In this case, res consists of just
    # one element,  rt itself, with the added attribute "random.factors"
    # which consists of the single string "none".
    #

    

    if (length(submat) == 0) {
        res[[1]] <- rt
        attr(res[[1]], "random.factors") <- "none"
        attr(res[[1]], "of.interest") <- T
        return(res)
    }

    # this loop parses until we are finished,
    # as indicated by all elements of "done"
    # being T.
    while (!(all(done))) {
        # look for the first element of `done' which has value
        # F.  This corresponds to the first column of 
        # mat which itself corresponds to a term which has not already
        # been parsed.  Then calculate a logical vector choose of 
        # length ncol (mat) which has choose[i]=T iff the term 
        # corresponding to column $i$ contains exactly the same random factors
        # as this term.   
        smallest <- match(F, done)
        matching <- mat[, smallest, drop = F][randomvars]
        z <- submat == matching[row(submat)]
        #
        choose <- apply(z, 2, all)
        # Now construct a matrix with rows 
        fixedsub <- mat[!(randomvars), choose, drop = F]
        form <- paste("1:", length(xx), " ~ 1 ")
        if (!(all(randomvars))) {
            for (j in (1:ncol(fixedsub))) {
                f <- fixedsub[, j, drop = F]
                if (any(f == 2)) {
                  first <- names(f[f == 2])
                  second <- names(f[f == 1])
                  term <- paste(first, "% in %", second)
                }
                else {
                  vars <- names(f[f == 1])
                  term <- paste(vars, collapse = ":")
                }
                if (term != "") {
                  form <- paste(form, " + ", term)
                }
            }
        }
        newtext <- paste("terms(", form, ")", sep = "")
        res[[l]] <- eval(parse(text = newtext))
        if (all(matching == 0)) {
            attr(res[[l]], "random.factors") <- "none"
            attr(res[[l]], "of.interest") <- T
        }
        else {
            attr(res[[l]], "random.factors") <- names(matching[matching == 
                1])
            attr(res[[l]], "of.interest") <- all(random.of.interest[matching])
        }
        attr(res[[l]],"only.ones") <- FALSE
        l <- l + 1
        done <- done | choose
    }

    if(all(randomvars)){
      #everything is now done, unless we're dealing with a formula
      #like
      #
      #
      #  1 + R + S +....
      #
      #  where all the terms are random factors.  In this case we need
      #  to put in a block for the intercetp term.  It is most logical for
      # this to be the first block, so we move all the other blocks one
      # object along in the list "res", so that we can "slip in" a block
      # in res[[1]] representing the intercept .

     for(i in seq(length(res),1,by=-1)){
    res[[2]] <- res[[1]]
    attr(res[[2]],"only.ones") <-  FALSE
  }

     #now put an object into res[[1]] representing the intercept term.
     #In the C code for putting the  block into a model, the first thing we
     #do is check that the "only.one" attribute is TRUE, and if it
     #is the C code deals with it specially, so any terms object with
     # attribute "only.ones" set equal to TRUE and with non-NULL attributes
     # "random.factors" and "of.interest" so as not to crash the C will
     #do fine.     
     res[[1]] <- terms(1~1)
     attr(res[[1]],"only.ones") <- TRUE
     attr(res[[1]],"random.factors") <- "none"
     attr(res[[1]],"of.interest") <- TRUE

   }

    #finished.
    
    res
}

"Ra"<-
function(data, shape = 0.001, scale = 0.001, type = "identity", contrast = 
	"treatment", map = "none",  of.interest = FALSE)
{
   	l <- length(levels(data))
	x <- l
	f1 <- data
	attr(f1, "type") <- type
        if (!(type%in%c("identity","adj.map"))){
          stop("type argument not appropiate")
        }
        attr(f1, "contrast") <- contrast
        if (!(contrast%in%c("identity", "sum", "treatment", 
				"helmert", "diff2"))){
          stop("contrast argument not appropiate")
        }       
	attr(f1, "map") <- map
        if ((!(file.exists(map))) && (!(map=="none"))){
          stop("map argument is the name of a non-existent file")
        }
	attr(f1, "scale") <- scale
	attr(f1, "shape") <- shape
	attr(f1, "is.random") <- T
	attr(f1, "levels") <- levels(as.factor(f1))
	attr(f1, "of.interest") <- of.interest
	f1
}
"model.init"<-
function(mm, family, model.show = T)
{
	y.values <- model.extract(mm, response)
        if(NCOL(y.values) == 2){
          response <- y.values[,1]/( y.values[,1]+ y.values[,2])
          weights <-  y.values[,1] + y.values[,2]
        }
        else{
          response <- y.values
          	weights <- model.extract(mm, weights)
        }
	offset <- model.extract(mm, offset)
        if(is.null(offset)) {offset <- 0}
        if (is.character(family))
          family <- get(family)
        if (is.function(family))
          family <- family()
		switch(family$family,
		"binomial" = fam <- 0,
		"poisson" = fam <- 1,
		"gaussian" = fam <- 2,
		stop("unkown familly"))
	scale.factor <- 1
	l <- .C("onemodel_model_init",
		as.integer(length(response)),
		as.double(response),
		as.integer(!(is.null(weights))),
		as.double(weights),
		as.integer(!(length(offset) == 1 && offset[1] == 0)),
		as.double(offset),
		as.integer(fam),
		as.double(scale.factor))
	invisible()
      }
"add.block"<-
function(rtpi, data)
{
        
	random.factors <- attr(rtpi, "random.factors")
        	if(random.factors == "none") {
		n.random.factors <- 0
	}
	else {
		n.random.factors <- length(random.factors)
	}
	random.nlevels <- numeric(0)
	random.levels <- numeric(0)
	random.describe <- numeric(0)
	random.map.files <- character(0)
	random.irregular.data <- numeric(0)
	contrast.describe <- numeric(0)
       
	shape <- 0
	scale <- 0
	if(n.random.factors > 0) {
		for(i in (1:n.random.factors)) {
			if(missing(data)) {
				xx <- eval(parse(text = random.factors[i]))
			}
			else {
				xx <- eval(parse(text = paste(c(as.character(
				  substitute(data)), "$", random.factors[i]), 
				  collapse = "")))
			}
			random.nlevels <- c(random.nlevels, length(levels(xx)))
			random.levels <- c(random.levels, as.numeric(xx))
			p.matrix <- attr(xx, "type")
			shape <- attr(xx, "shape")
			scale <- attr(xx, "scale")
			random.describe <- c(random.describe, match(p.matrix, c(
				"identity", "diff1", "diff2", "adj.map")))
			random.map.files <- c(random.map.files, (attr(xx, "map"
				)))
			contrast.describe <- c(contrast.describe, match(attr(xx,
				"contrast"), c("identity", "sum", "treatment", 
				"helmert", "diff2")))
       
		} 
	}
	m <- model.matrix(rtpi, data = data)
     
	of.interest <- attr(rtpi, "of.interest")
	if(  (ncol(m) == 1)   && (n.random.factors > 0)) {
		arg1 <- 0
		arg2 <- 0
		arg3 <- 0
	}
	else {
          arg1 <- nrow(m)
		arg2 <- ncol(m)
		arg3 <- m
	}

                only.ones <- attr(rtpi,"only.ones")
	.C("r_to_block",
		as.integer(nrow(m)),
		as.integer(arg2),
		as.double(arg3),
		as.integer(n.random.factors),
		as.integer(random.nlevels),
		as.integer(random.levels),
		as.integer(random.describe),
		as.character(random.map.files),
		as.integer(contrast.describe),
		as.character(dimnames(m)[[2]]),
		as.character(random.factors),
		as.double(shape),
		as.double(scale),
		as.integer(of.interest),
                as.integer(only.ones))
	invisible()
}
"sampled.values" <- 
function (g, names, drop = T) 
{
    model.parameter.names <- attr(g, "parameter.names")
    if ((length(names) == 1) && (drop == T)) {
        num <- match(names, model.parameter.names)
        if (is.na(num)) {
            stop(paste(names, "is not the name of a parameter in the model"))
        }
        res <- g[[1]][, (num + 1)]
    }
    else {
        res <- matrix(0, nr = attr(g, "keep"), nc = length(names))
        i <- 1
        for (n in (names)) {
            num <- match(n, model.parameter.names)
            if (is.na(n, model.parameter.names)) {
                stop(paste(n, "is not the name of a parameter in the model"))
            }
            res[, i] <- g[[1]][, (num + 1)]
            i <- i + 1
        }
    }
    res
}

"plot.glmmfit" <-
  function (g, parameters, type = "b", hscale = "b", plots.per.page = 4, 
            postscript) 
{
  inner.rep <- function(x,k){rep(x,rep(k,(length(x))))}
###The plot function for objects of class glmmfit, i.e.
###an object returned gy the glmm function.
  object.name <- deparse(substitute(g))
  
  ##Firstly we check that type and hscale have possible values
  ##(see the documentation)
  possible.type.values <- c("b", "B", "both", "Both", "t", 
                            "T", "trace", "Trace", "h", "H", "histogram", "Histogram")
  possible.hscale.values <- c("b", "B", "both", "Both", "s", 
                              "S", "sigma", "Sigma", "t", "T", "tau", "Tau")
  if (!((type%in%possible.type.values))) {
    stop("type argument not appropiate")
  }
  if (!((hscale%in%possible.hscale.values))) {
    stop("hscale argument not appropiate")
  }

  ## grab from the glmm fit the following parameters
  ne <- attr(g, "number.of.effects")
  nh <- attr(g, "number.of.hyperparameters")
  first.iteration <- attr(g, "burnin") + 1
  last.iteration <- attr(g, "burnin") + attr(g, "keep")
  thin <- attr(g, "thin")
  
  ##
  ##The first real task is to fill construct the vectors
  ##`plot.names', `plot.numbers', and `plot scale', all of length
  ## equal to the number variables to be plotted. These
  ##contain, respectively, the names of the parameters,
  ##the column number within the matrix which forms the
  ##g[[1]] where g is the glmm fit, and a vector of scalars
  ##which take the value 0 if the variable is to be plotted
  ##on the identity scale, and 1 if it is to be plotted on the
  ##1/sqrt scale.
  ##
  ##how we contruct these depends on whether the argument
  ##`parameters' isgiven a value.
  if (!(missing(parameters))) {
    ##
    ##if the parameters argument is given a value, it is a character vector of
    ##the names of those parameters in the glmm fit which are to be plotted.
    ##
    ##

    
    
    n <- length(parameters)
    plot.names <- character(0)
    plot.numbers <- numeric(0)
    plot.scale <- numeric(0)
    ##So for each element of `parameter.names' ....

    for (i in 1:n) {
      ##check that it is the name of a parameter
      ##in the glmm fit, and put in `num' the variable number
      ##(ie the position in the list of variable names)
      name <- parameters[i]
      num <- match(name, attr(g, "parameter.names"))
      if (is.na(num)) {
        stop(paste(parameters[i], "is not the name of a parameter in the model"))
      }
      ##add to `plot.numbers' the relevent column nnumber in g[[1]].
      ##This is num+1 because g[[1]][1,] is a vector of iteration numbers. 
      plot.numbers <- c(plot.numbers, (num + 1))
      ##add to `plot.names' the name of the parameter.
      plot.names <- c(plot.names, name)
      ##add to `plot.scale' the value 0, indicating plotting on the
      ##identity scale.
      plot.scale <- c(plot.scale, 0)
      if (num > ne) {
        if ((hscale%in% c("b", "both", "B", "Both"))) {
          ##if the user has asked for hyperperameters
          ##to be plotted on both the identity scale and
          ##on the 1/sqrt scale, then any hyperparameters have
          ##to be added again, by adding the same values to
          ##`plot.numbers' and `plot.names' but the values
          ## 1 to `plot.scale', indicating plotting on the
          ## 1/sqrt scale.
          plot.numbers <- c(plot.numbers, (num + 1))
          plot.names <- c(plot.names, name)
          plot.scale <- c(plot.scale, 1)
        }
      }
    }
    number.of.plots <- length(plot.names)
  }
  if (missing(parameters)) {
    if ((hscale%in% c("b", "both", "B", "Both"))) {
      number.of.plots <- ne + 2 * nh
      plot.names <- character(number.of.plots)
      plot.numbers <- numeric(number.of.plots)
      plot.scale <- numeric(number.of.plots)
      if(ne > 0){
        plot.names[1:ne] <- attr(g, "parameter.names")[1:ne]
        plot.numbers[1:ne] <- (1:ne) + 1
      }
      if (nh > 0) {
        plot.names[(ne + 1):number.of.plots] <- inner.rep(attr(g, 
                                                               "parameter.names")[(ne + 1):(ne + nh)], 2)
        suffices <- rep(c("(tau scale)", "(sigma scale)"), 
                        nh)
        plot.names[(ne + 1):number.of.plots] <- paste(plot.names[(ne + 
                                                                  1):number.of.plots], suffices)
        plot.numbers[(ne + 1):(ne + (2 * nh))] <- (inner.rep(((ne + 
                                                               1):(ne + nh)), 2)) + 1
        plot.scale[(ne + 1):number.of.plots] <- rep((0:1), 
                                                    nh)
      }
    }
    if ((hscale%in% c("t", "tau", "T", "Tau"))) {
      number.of.plots <- ne + nh
      plot.names <- character(number.of.plots)
      plot.names[1:ne] <- attr(g, "parameter.names")[1:ne]
      if (nh > 0) {
        plot.names[(ne + 1):number.of.plots] <- attr(g, 
                                                     "parameter.names")[(ne + 1):(ne + nh)]
        suffices <- rep(c("(tau scale)"), nh)
        plot.names[(ne + 1):number.of.plots] <- paste(plot.names[(ne + 
                                                                  1):number.of.plots], suffices)
      }
      plot.numbers <- (1:(ne + nh)) + 1
      plot.scale <- c(rep(0, ne), rep((0), nh))
    }
    if ((hscale%in%c("s", "sigma", "S", "Sigma"))) {
      number.of.plots <- ne + nh
      plot.names <- character(number.of.plots)
      plot.names[1:ne] <- attr(g, "parameter.names")[1:ne]
      if (nh > 0) {
        plot.names[(ne + 1):number.of.plots] <- attr(g, 
                                                     "parameter.names")[(ne + 1):(ne + nh)]
        suffices <- rep(c("(sigma scale)"), nh)
        plot.names[(ne + 1):number.of.plots] <- paste(plot.names[(ne + 
                                                                  1):number.of.plots], suffices)
      }
      plot.numbers <- (1:(ne + nh)) + 1
      plot.scale <- c(rep(0, ne), rep((1), nh))
    }
  }
  number.of.pages <- number.of.plots%/%plots.per.page
  if ((number.of.plots%%plots.per.page) > 0) {
    number.of.pages <- number.of.pages + 1
  }
  finish <- FALSE
  i<-1
  while(!(finish)) {
    first <- ((i - 1) * plots.per.page) + 1
    last <- min((i * plots.per.page), number.of.plots)
    if (!(missing(postscript))) {
      postscript(paste(postscript, "_page", i, ".ps" , 
                       sep = ""))
    }
    if ((type%in%c("b", "B", "both", "Both"))) {
      number.of.columns <- 2
    }
    else {
      number.of.columns <- 1
    }
    par(mfrow = c((last - first) + 1, number.of.columns))
    for (j in (first:last)) {
      data <- g[[1]][, plot.numbers[j]]
      if (plot.scale[j] == 1) {
        data <- 1/sqrt(data)
      }
      if ((type%in% c("b", "B", "both", "Both", "t", "T", 
                      "trace", "Trace"))) {
        plot(seq(first.iteration, last.iteration, by = thin), 
             data, type = "l", xlab = "iteration", ylab = plot.names[j])
      }
      if ((type%in%c("b", "B", "both", "Both", "h", "H", 
                     "histogram", "Histogram"))) {
        hist(data, xlab = plot.names[j], main = paste("Histogram of", 
                                           plot.names[j]))
      }
    }
  
    finish <- FALSE
    if(i == number.of.pages){
      finish <- TRUE
    }
    else{
      if(missing(postscript)){
      xx <- menu(c("More Variables", "Exit"))
      if (xx == 2){
        finish <- TRUE
      }
    }
      if(!(missing(postscript))){
        finish <- FALSE
      }
    }
    if(finish){
      if(!(missing(postscript))){dev.off()}
    }
    else i<-i+1
  }
}
"print.glmmfit" <-
function (g) 
{
    #print out the important information about a glmm fit. 
    cat("\n")
    cat("\n")
    ne <- attr(g, "number.of.effects")
    nh <- attr(g, "number.of.hyperparameters")
    column.names <- c("mean", "st.dev", " 2.5%", "1st Qu.", "Median", 
        "3rd Qu.", " 97.5%")
    cat("glmm output\n")
    cat("===========\n")
    cat("call:\n")
    cl <- attr(g, "call")
    print(cl)
    cat("\n")
    cat(attr(g, "icm"), " steps of ICM \n")
    cat(attr(g, "burnin"), " iterations of burnin \n")
    if (attr(g, "thin") == 1){
    cat(attr(g, "keep"), " stored iterations \n")
  }
    if (attr(g,"thin") > 1){
      x <- attr(g, "thin")
    cat("Then a further", attr(g,"keep"),"iterations \n")
    cat("storing every",x,ordinal.end(x),"\n")
     }
    cat("\n")
    cat("\n")
    if (ne > 0) {
        cat("effect parameters\n")
        cat("-----------------\n")
        effect.parameters <- attr(g, "stats")[1:ne, , drop = F]
        dimnames(effect.parameters) <- list(attr(g, "parameter.names")[1:ne], 
            column.names)
        print(effect.parameters, digits = 3)
    }
    if (nh > 0) {
        hyperparameters.tau <- attr(g, "stats")[(ne + 1):(ne + 
            nh), , drop = F]
        dimnames(hyperparameters.tau) <- list(attr(g, "parameter.names")[(ne + 
            1):(ne + nh)], column.names)
        cat("\n")
        cat("hyperparameters (tau scale)\n")
        cat("--------------------------\n")
        print(hyperparameters.tau, digits = 3)
        hyperparameters.sigma <- attr(g, "stats")[(ne + nh + 
            1):(ne + 2 * nh), , drop = F]
        dimnames(hyperparameters.sigma) <- list(attr(g, "parameter.names")[(ne + 
            1):(ne + nh)], column.names)
        cat("\n")
        cat("hyperparameters (sigma scale)\n")
        cat("-----------------------------\n")
        print(hyperparameters.sigma, digits = 3)
    }
}

coda.object <- function(g){
  x <- g[[1]][,-1]
  dimnames(x) <- list(NULL,attr(g,"parameter.names"))
  attr(x,"mcpar") <- c(1,(attr(g,"keep")/(attr(g,"thin"))),1)
  attr(x,"class") <- "mcmc"
  x
}
 bugs.out <-
function(g,file){
  indfile.name <- paste(file,".ind",sep="")
  outfile.name <- paste(file,".out",sep="")
  .C("bugs_output",
     as.integer(length(attr(g,"parameter.names"))),
     as.integer(attr(g,"keep")/attr(g,"thin")),
     as.character(attr(g,"parameter.names")),
     as.double(as.numeric(g[[1]])),
     as.character(indfile.name),
     as.character(outfile.name)
     )
     invisible()
}





