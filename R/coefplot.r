siglev <- function(p)
{
  if (p <= 0.1)
    if (p <= 0.05)
      if (p <= 0.01)
        if (p <= 0.001)
          sig <- "***"
  else
    sig <- "**"
  else
    sig <- "*"
  else
    sig <- "."
  else
    sig <- ""
}

siglev <- function(p)
{
  if (p <= 0.1)
    if (p <= 0.05)
      if (p <= 0.01)
        if (p <= 0.001)
          sig <- "***"
  else
    sig <- "**"
  else
    sig <- "*"
  else
    sig <- "."
  else
    sig <- ""
}

# This functions makes a plot of the COEFFICIENTS of significant factors.
plot.coef <- function(glm, alpha = 0.05, pr="Pr(>|z|)", cex.labels = 1, pos.labels = 0.4, ...)
{
  # Extract p values and coefs.
  p <- coef(summary(glm))[,pr]
  c <- coef(glm)
  
  # Separate intercept from other coefs.
  p.int <- as.numeric(p[1])
  c.int <- as.numeric(c[1])
  p <- p[-1]
  c <- c[-1]
  
  # Find the order in which to put coefficients.
  c <- c[which(p <= alpha)]
  p <- p[which(p <= alpha)]
  
  # Find the order according to coef.
  order <- order(c)
  c <- c[order]
  p <- p[order]
  names <- names(c)
  
  # Put intercept and coefs back together.
  names <- c("(Intercept)", names)
  vals <- c(round(c.int, 3), round(c, 3))
  
  c <- c(as.numeric(c.int), as.numeric(c))
  p <- c(as.numeric(p.int), as.numeric(p))
  vals <- paste(vals, lapply(p, siglev), sep=" ")
  names <- paste(names, vals, sep="\n")
  
  range <- max(c)-min(c)
  ylim <- c(min(c)-0.25*range, max(c)+0.25*range)
  bars <- barplot(c, ylim=ylim, col=c("white", rep("lightgray", length(c)-1)), ...)
  pos <- ifelse(c > 0, c+pos.labels*cex.labels, c-pos.labels*cex.labels)
  text(bars, pos, names, cex=cex.labels)
}

# This functions makes a plot of the ODDS RATIOS of significant factors.
plot.or <- function(glm, alpha = 0.05, pr="Pr(>|z|)", cex.labels = 1, pos.labels = 3, ...)
{
  # Extract p values and coefs.
  p <- coef(summary(glm))[,pr]
  c <- exp(coef(glm))
  
  # Separate intercept from other coefs.
  p.int <- as.numeric(p[1])
  c.int <- as.numeric(c[1])
  p <- p[-1]
  c <- c[-1]
  
  # Extract the significant coefficients.
  c <- c[which(p <= alpha)]
  p <- p[which(p <= alpha)]
  
  # Find the order according to coef.
  order <- order(c)
  c <- c[order]
  p <- p[order]
  names <- names(c)
  
  # Put intercept and coefs back together.
  names <- c("(Intercept)", names)
  vals <- c(round(c.int, 3), round(c, 3))
  
  c <- c(as.numeric(c.int), as.numeric(c))
  p <- c(as.numeric(p.int), as.numeric(p))
  vals <- paste(vals, lapply(p, siglev), sep=" ")
  names <- paste(names, vals, sep="\n")
  
  range <- max(c)-min(c)
  ylim <- c(0, max(c)+0.25*range)
  bars <- barplot(c, ylim=ylim, col=c("white", rep("lightgray", length(c)-1)), ...)
  pos <- ifelse(c > 0, c+pos.labels*cex.labels, c-pos.labels*cex.labels)
  text(bars, pos, names, cex=cex.labels)
}

lr.test <- function(glm, glm0)
{
  ll <- -2*logLik(glm)
  ll0 <- -2*logLik(glm0)
  lr <- ll0 - ll
  df <- glm$rank - glm0$rank
  p <- pchisq(lr, df)
  list(lr=as.numeric(lr), df=as.numeric(df), p=as.numeric(p))
}

phi.glm <- function(glm)
{
  sum(resid(glm, type="pearson")^2 / df.residual(glm))
}
