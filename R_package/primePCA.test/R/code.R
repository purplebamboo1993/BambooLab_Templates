#' @import softImpute
#' @import Matrix
#' @import MASS
#' @importFrom methods as
#' @importFrom stats runif
NULL

#' Frobenius norm sin theta distance between two column spaces
#'
#' @param V1 a matrix with orthonormal columns
#' @param V2 a matrix of the same dimension as V1 with orthonormal columns 
#' @return the Frobenius norm sin theta distance between two V1 and V2
#' @export
sin_theta_distance = function(V1, V2){
    K = min(dim(V1)[2], dim(V2)[2])
    sigmas = svd(t(V1) %*% V2)$d
    return(sqrt(K - sum(sigmas ^ 2)))
}

sample_filter_op = function(thresh, V, Omega, prob=1){
    # Omega is RsparseMatrix

    d = dim(V)[1]
    K = dim(V)[2]
    n = dim(Omega)[1]
    selected = sapply(1:n, function(i){

        if (Omega@p[i] == Omega@p[i + 1]){
            obs_num = 0                             
        }
        else{
            omega_i_col = Omega@j[(Omega@p[i] + 1):Omega@p[i + 1]] + 1
            obs_num = length(omega_i_col)            
        }

        if (obs_num <= K){
            return(FALSE)
        }
        else {
            # print(i)
            # print(omega_i_col)
            # print(V[omega_i_col, ])
            sigmas = svd(V[omega_i_col, ])$d
            if (tail(sigmas, 1) >=  sqrt(obs_num / d) / thresh){
                rv = runif(1)
                if (rv < prob){
                    return(TRUE)                                               
                }
                else{
                    return(FALSE)
                }
            }
            else {
                return(FALSE)
            }
        }
    })
    return(which(selected))

}

get_omega = function(X){

    Omega = X
    num_nonzero = length(Omega@x)    
    Omega@x = rep(1, num_nonzero)    
    return(Omega)
}


#' Center and/or normalize each column of a matrix
#' @param X a numeric matrix with NAs or "Incomplete" matrix object (see softImpute package) 
#' @param center center each column of \code{X} if \code{center == TRUE}. 
#' The default value is \code{TRUE}. 
#' @param normalize normalize each column of \code{X} such that its sample variance is 1 if \code{normalize == TRUE}. 
#' The default value is \code{False}. 
#' @return a centered and/or normalized matrix of the same dimension as \eqn{X}.
#' @export 
col_scale = function(X, center=T, normalize=F){

    X_center = as(X, "Incomplete")
    X_center = biScale(X_center, row.center=FALSE, row.scale=FALSE, col.center=center, col.scale=normalize)
    X_center = as(X_center, "CsparseMatrix")    
    return(X_center)
}


#' Inverse probability weighted method for estimating the top K eigenspaces
#' @param X a numeric matrix with \eqn{NA}s or "Incomplete" matrix object (see softImpute package)
#' @param K the number of principal components of interest
#' @param trace.it report the progress if \code{trace.it == TRUE}
#' @param center center each column of \code{X} if \code{center == TRUE}. 
#' The default value is \code{TRUE}. 
#' @param normalize normalize each column of \code{X} such that its sample variance is 1 if \code{normalize == TRUE}. 
#' The default value is \code{False}. 
#' @return Columnwise centered matrix of the same dimension as \eqn{X}.
#' @export 
#' @examples
#' X = matrix(1:30 + .1 * rnorm(30), 10, 3)
#' X[1, 1] = NA
#' X[2, 3] = NA
#' v_hat = inverse_prob_method(X, 1)

inverse_prob_method = function(X, K, trace.it=F, center=T, normalize=F){
    # X is a regular matrix with NAs or Incomplete matrix. 

    if (all_na_column(X)){
        cat("There exists at least one all-NA column. Please screen this out first.\n")
        return(NULL)
    }
    if ((any(apply(!is.na(X), 2, sum) == 1)) & (normalize)){
        cat("There exists at least one column with only one NA, so that the normalisation cannot be applied. Please screen this out first.\n")        
        return(NULL)
    }
    X_center = col_scale(X, center, normalize)
    # if (!col_centered){
    #     X_center = col_center(X)
    # }
    # else{
    #     X_center = as(X, "Incomplete")
    # }
    Omega = get_omega(X_center)
    Sigma_x = t(X_center) %*% X_center
    Sigma_o = t(Omega) %*% Omega
    weight = 1 / Sigma_o@x
    weight[weight == Inf] = 0
    Sigma_tilde = Sigma_x
    Sigma_tilde@x = Sigma_x@x * weight
    Sigma_tilde = as(Sigma_tilde, "sparseMatrix")
    V_hat = svd.als(Sigma_tilde, rank.max=K, trace.it=trace.it, maxit=1000)$v
    return(V_hat)

}


complete = function(V_cur, X, Omega, K, trace.it=F){
    # X and Omega are RsparseMatrix
    # return all the principal scores as a list U_hat

    n = dim(X)[1]
    results = sapply(1:n, function(i){
        omega_col = Omega@j[(Omega@p[i] + 1):Omega@p[i + 1]] + 1
        x = X@x[(Omega@p[i] + 1):Omega@p[i + 1]]
        V_part = V_cur[omega_col, ]
        u_hat = ginv(V_part) %*% x
        return(list(u=u_hat, residual=x - V_part %*% u_hat))
    }, simplify=F)        
    U_hat = sapply(results, function(i) i$u)
    if (K > 1){
        U_hat = t(U_hat)
    }
    else{
        U_hat = matrix(U_hat, ncol=1)
    }
    residual_all = sapply(results, function(i) i$residual, simplify=F)
    return(list(u=U_hat, res=residual_all))

}


select_rows = function(X, rows){
    # X is a RsparseMatrix

    index = do.call(c, sapply(rows, function(row) (X@p[row] + 1):X@p[row + 1], simplify=F))
    Xs = X
    Xs@j = X@j[index]
    Xs@x = X@x[index]
    obs_per_row = sapply(rows, function(row) X@p[row + 1] - X@p[row])
    Xs@p = as.integer(c(0, cumsum(obs_per_row)))
    Xs@Dim[1] = length(rows)
    return(Xs)

}


select_columns = function(X, cols){
    # X is a CsparseMatrix

    index = do.call(c, sapply(cols, function(col) (X@p[col] + 1):X@p[col + 1], simplify=F))
    Xs = X
    Xs@i = X@i[index]
    Xs@x = X@x[index]
    obs_per_col = sapply(cols, function(col) X@p[col + 1] - X@p[col])
    Xs@p = as.integer(c(0, cumsum(obs_per_col)))
    Xs@Dim[2] = length(cols)
    return(Xs)

}

eigen_refine = function(V_cur, X, Omega, svd_cur=NULL, trace.it=F, thresh=1e-10){
    # X and Omega are RsparseMatrix

    n <- dim(X)[1]
    d <- dim(V_cur)[1]
    K <- dim(V_cur)[2]
    results <- complete(V_cur, X, Omega, K, trace.it=trace.it)
    U_hat <- results$u
    residual_all <- results$res
    S_hat <- X
    S_hat@x <- as.vector(do.call(rbind, residual_all))
    S_hat <- as(S_hat, "CsparseMatrix")
    X_hat <- splr(S_hat, U_hat, V_cur)
    if (is.null(svd_cur)){
        svd_update <- svd.als(X_hat, rank.max=K, trace.it=trace.it, thresh=thresh, maxit=1000)
    }
    else {
        svd_update <- svd.als(X_hat, rank.max=K, trace.it=trace.it, thresh=thresh, warm.start=svd_cur, maxit=1000)
    }
    # V_update <- svd.als(X_hat, rank.max=K, trace.it=trace.it, maxit=1000)$v
    return(svd_update)

}


mse_eval = function(V_new, X, Omega){
    K = dim(V_new)[2]
    residual_all = complete(V_new, X, Omega, K)$res
    res_vector = as.vector(do.call(rbind, residual_all))
    return(mean(res_vector ^ 2))
}


all_na_column = function(X){

    X1 = as(X, "Incomplete")
    X1 = as(X1, "CsparseMatrix")
    if (length(X1@p) > length(unique(X1@p))){
        return(TRUE)
    }
    else{
        return(FALSE)
    }
}


#' primePCA algorithm
#'
#' @param X an \eqn{n}-by-\eqn{d} data matrix with \code{NA} values
#' @param K the number of the principal components of interest
#' @param V_init an initial estimate of the top \eqn{K} eigenspaces of the covariance matrix of \code{X}.
#' By default, primePCA will be initialized by the inverse probability method. 
#' @param thresh_sigma used to select the "good" rows of \eqn{X} to update the principal eigenspaces (\eqn{\sigma_*} in the paper). 
#' @param max_iter maximum number of iterations of refinement 
#' @param thresh_convergence The algorithm is halted if the Frobenius-norm sine-theta distance between the two consecutive iterates
#' @param thresh_als This is fed into \code{thresh} in \code{svd.als} of
#' \code{softImpute}.   
#' is less than \code{thresh_convergence}. 
#' @param trace.it report the progress if \code{trace.it} = \code{TRUE}
#' @param prob probability of reserving the "good" rows. \code{prob == 1} means to reserve all the "good" rows. 
#' @param save_file the location that saves the intermediate results, including \code{V_cur}, \code{step_cur} and \code{loss_all}, 
#' which are introduced in the section of returned values. The algorithm will not save any intermediate result 
#' if \code{save_file == ""}. 
#' @param center center each column of \code{X} if \code{center == TRUE}. 
#' The default value is \code{TRUE}. 
#' @param normalize normalize each column of \code{X} such that its sample variance is 1 if \code{normalize == TRUE}. 
#' The default value is \code{False}. 
#' @return a list is returned, with components \code{V_cur}, \code{step_cur} and \code{loss_all}. 
#' \code{V_cur} is a \eqn{d}-by-\eqn{K} matrix of the top \eqn{K} eigenvectors. \code{step_cur} is the number of iterations. 
#' \code{loss_all} is an array of the trajectory of MSE. 
#' @export 
#' @examples
#' X = matrix(1:30 + .1 * rnorm(30), 10, 3)
#' X[1, 1] = NA
#' X[2, 3] = NA
#' v_tilde = primePCA(X, 1)$V_cur
primePCA = function(X, K, V_init=NULL, thresh_sigma=10, max_iter=1000, thresh_convergence=1e-5, thresh_als=1e-10, trace.it=F, prob=1, save_file="", center=T, normalize=F){
    # X: a regular matrix with NAs or an Incomplete matrix 
    # save_file = "" means that you will not save intermediate step during the optimization procedure
    # V_init: d-by-K orthonormal matrix. 

    if (all_na_column(X)){
        cat("There exists at least one all-NA column. Please screen this out first.\n")
        return(NULL)
    }
    if ((any(apply(!is.na(X), 2, sum) == 1)) & (normalize)){
        cat("There exists at least one column with only one observation, so that the normalisation cannot be applied. Please screen this out first.\n")        
        return(NULL)
    }
    sfn = save_file        
    X_center = col_scale(X, center, normalize)
    loss_all = rep(Inf, max_iter)
    if (is.null(V_init)){
        V_cur = inverse_prob_method(X_center, K, center=center, normalize=normalize)
    }
    else{
        V_cur = V_init        
    }    
    X_center = as(X_center, "RsparseMatrix")
    Omega = get_omega(X_center)
    i = 0
    svd_cur <- NULL
    while (i < max_iter){   
        i = i + 1
        rows = sample_filter_op(thresh_sigma, V_cur, Omega, prob=prob)
        Xs = select_rows(X_center, rows)        
        Os = select_rows(Omega, rows)
        if (is.null(svd_cur)){
            svd_new <- eigen_refine(V_cur, Xs, Os, thresh=thresh_als)
        }
        else {
            svd_new <- eigen_refine(V_cur, Xs, Os, svd_cur=svd_cur, thresh=thresh_als)
        }
        V_new <- svd_new$v
        difference = sin_theta_distance(V_new, V_cur)        
        loss_cur = mse_eval(V_new, Xs, Os)  
        loss_all[i] = loss_cur
        if (trace.it){
            cat("Step ", i, "\n", dim(Xs)[1], " rows selected\n", "MSE: ", loss_cur, "\n")  
            cat("F-norm sine theta difference: ", difference, "\n")                  
        }
        svd_cur <- svd_new
        V_cur <- svd_cur$v
        step_cur = i        
        if (difference < thresh_convergence){
            cat("Convergence threshold is hit.\n")
            break
        }
        if (save_file != ""){
            save(file=save_file, V_new, step_cur, loss_all)            
        }
    }
    if (i >= max_iter){
        cat("Max iteration number is hit.\n")
    }
    return(list(V_cur=V_cur, step_cur=step_cur, loss_all=loss_all[1:step_cur]))
}
