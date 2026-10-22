# Copyright 2026 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

ordinal_data <- function(K=4L) {
  y <- matrix(rep(seq_len(K)-1L, length.out=48L), 12L, 4L)
  y[1,1] <- NA_real_
  z <- matrix(0,12L,4L); z[1:6,3:4] <- 1
  list(y=y, x=matrix(seq_len(24)/24,12L,2L), z=z, t=1:4)
}
ordinal_options <- list(outcome="ordinal",num_categories=4L,
  num_burnin=2L,num_sweeps=3L,n_skip=1L,num_chains=2L,
  num_trees_pr=2L,num_trees_trt=2L,max_depth_pr=3L,max_depth_trt=3L,
  min_points_per_leaf_pr=2L,min_points_per_leaf_trt=2L,random_seed=931L)

test_that("ordinal R arguments never truncate counts or labels", {
  d <- ordinal_data()
  for (K in list(NULL, TRUE, 1, 2.5, Inf, "3", NA_real_)) {
    expect_error(do.call(longbet,c(d,list(outcome="ordinal",num_categories=K))), "num_categories")
  }
  for (s in list(TRUE,0,-1,Inf,NA_real_,"5")) {
    expect_error(do.call(longbet,c(d,list(cutpoint_prior_scale=s))), "cutpoint_prior_scale")
  }
  expect_error(do.call(longbet,c(d,list(num_categories=3))), "nonordinal")
  for (bad in c(-1,1.5,4,Inf,-Inf)) {
    bad_d <- d; bad_d$y[2,2] <- bad
    expect_error(do.call(longbet,c(bad_d,ordinal_options)), "integer labels|infinity")
  }
  bad_d <- d; bad_d$y[,] <- NA_real_
  expect_error(do.call(longbet,c(bad_d,ordinal_options)), "no observed")
  bad_d$y <- matrix(as.character(d$y),12L,4L)
  expect_error(do.call(longbet,c(bad_d,ordinal_options)), "numeric category")
  bad_d$y <- factor(rep(0:3,12L))
  expect_error(do.call(longbet,c(bad_d,ordinal_options)), "numeric category")
})

test_that("ordinal probability and ATT methods agree exactly with Python", {
  skip_without_engine()
  d <- ordinal_data()
  fit <- do.call(longbet,c(d,ordinal_options))
  expect_identical(fit$num_categories,4L)
  expect_equal(fit$meany,0); expect_equal(fit$sdy,1)
  expect_output(print(fit),"latent probit")
  np <- reticulate::import("numpy",convert=TRUE)
  for (summary_only in c(FALSE,TRUE)) {
    p <- predict(fit,d$x,d$z,summary_only=summary_only,alpha=.1)
    expect_equal(dim(p$prob_y_summary$mean),c(12L,4L,4L))
    expect_equal(dim(p$cutpoints_samples),c(6L,2L))
    expect_equal(dim(p$att_prob_full),c(2L,4L,6L))
    expect_output(print(p),"category probabilities")
    if (summary_only) {
      expect_null(p$prob_y)
      expect_error(predict_probabilities(p,summary=FALSE),"discarded")
    } else {
      expect_equal(dim(p$prob_y),c(12L,4L,4L,6L))
      expect_equal(predict_probabilities(p,summary=FALSE),p$prob_y,tolerance=0)
    }
    for (arm in c("factual","control","effect")) {
      got <- predict_probabilities(p,arm)
      expected <- p$py_pred$predict_probabilities(arm=arm)
      for (name in names(got)) expect_equal(got[[name]],as.array(np$asarray(reticulate::py_get_attr(expected,name))),tolerance=0)
    }
    category <- att_probabilities(p,.1)
    expected <- p$py_pred$att_probabilities(alpha=.1)
    for (name in c("att","intervals","att_full")) expect_equal(category[[name]],as.array(np$asarray(expected[[name]])),tolerance=0)
    score <- att_expected_score(p,c(0,0,0,1),.1)
    expected <- p$py_pred$att_expected_score(weights=c(0,0,0,1),alpha=.1)
    expect_equal(score$att_full,as.array(np$asarray(expected[["att_full"]])),tolerance=0)
    expect_equal(att_expected_score(p,rep(1,4))$att_full,matrix(0,2,6),tolerance=1e-14)
    expect_error(effect_draws(p),"category or score")
    # All probability work still calls Python after a dangling cached pointer.
    saved <- unserialize(serialize(p,NULL)); saved$py_pred <- NULL
    expect_equal(predict_probabilities(saved),predict_probabilities(p),tolerance=0)
    expect_equal(att_probabilities(saved),category <- att_probabilities(p),tolerance=0)
    expect_equal(att_expected_score(saved,c(0,0,0,1)),att_expected_score(p,c(0,0,0,1)),tolerance=0)
  }
  saved_fit <- unserialize(serialize(fit,NULL)); saved_fit$handle$py_model <- NULL
  expect_equal(predict(saved_fit,d$x,d$z)$prob_y,predict(fit,d$x,d$z)$prob_y,tolerance=0)
})

test_that("ordinal binary category and singleton axes survive R wrapping", {
  skip_without_engine()
  d <- ordinal_data(2L)
  options <- ordinal_options
  options$num_categories <- 2L; options$num_sweeps <- 1L; options$num_chains <- 1L
  options$random_intercept <- FALSE
  fit <- do.call(longbet,c(d,options))
  # Python warnings are emitted on stderr rather than as R warning conditions.
  p <- predict(fit,d$x[1,,drop=FALSE],matrix(0,1,1),t=1)
  expect_equal(dim(p$prob_y),c(1L,1L,2L,1L))
  expect_equal(dim(p$cutpoints_samples),c(1L,0L))
  expect_equal(dim(predict_probabilities(p)$mean),c(1L,1L,2L))
  expect_equal(dim(att_probabilities(p)$att_full),c(1L,2L,1L))
  expect_equal(dim(att_expected_score(p)$att_full),c(1L,1L))
  expect_true(all(is.nan(att_probabilities(p)$att)))
})

test_that("ordinal multi counts and shared parent/child rehydration agree", {
  skip_without_engine()
  d <- ordinal_data()
  yy <- list(cont=matrix(sin(seq_len(48)),12,4),ord3=d$y %% 3,
             binary=d$y %% 2,ord4=d$y)
  yy$cont[is.na(d$y)] <- NA_real_
  opts <- ordinal_options
  opts$outcome <- c("continuous","ordinal","binary","ordinal")
  opts$num_categories <- list(cont=NULL,ord3=3L,binary=NULL,ord4=4L)
  opts$sigma_prior_a <- 2; opts$sigma_prior_b <- 1
  for (bad in list(list(3,3,NULL,4),list(NULL,3.5,NULL,4),list(ord3=3,ord4=4),TRUE)) {
    bad_opts <- opts; bad_opts$num_categories <- bad
    expect_error(do.call(longbet_multi,c(list(y=yy,x=d$x,z=d$z),bad_opts)),"num_categories|Nonordinal")
  }
  fit <- do.call(longbet_multi,c(list(y=yy,x=d$x,z=d$z),opts))
  expect_equal(fit$order,c(2L,3L,4L,1L))
  expect_identical(fit$num_categories,opts$num_categories)
  p <- predict(fit,d$x,d$z)
  expect_equal(dim(p["ord3"]$prob_y),c(12L,4L,3L,6L))
  expect_equal(dim(p["ord4"]$prob_y),c(12L,4L,4L,6L))
  expect_equal(effect_draws(p,"cont"),p["cont"]$tauhats)
  expect_error(effect_draws(p,"ord3"),"category or score")
  for (name in c("ord3","ord4")) {
    expect_equal(att_probabilities(p[name])$att_full,p[name]$att_prob_full,tolerance=0)
    child <- unserialize(serialize(fit[name],NULL)); child$handle$py_model <- NULL
    expect_equal(predict(child,d$x,d$z)$prob_y,p[name]$prob_y,tolerance=0)
  }
  # A real fresh R process exercises raw_model and saved prediction arrays.
  payload <- tempfile(fileext=".rds"); output <- tempfile(fileext=".rds")
  on.exit(unlink(c(payload,output)),add=TRUE)
  saveRDS(list(fit=fit,pred=p,d=d),payload)
  loader <- if (pkgload::is_dev_package("longbet")) {
    sprintf("pkgload::load_all(%s,quiet=TRUE)",deparse(getNamespaceInfo(asNamespace("longbet"),"path")))
  } else "library(longbet)"
  script <- sprintf(paste0("%s; a<-readRDS(%s); p<-predict(a$fit,a$d$x,a$d$z,summary_only=TRUE); ",
    "saveRDS(list(prob=predict_probabilities(p['ord4']),score=att_expected_score(a$pred['ord4'],c(0,0,0,1))),%s)"),
    loader,deparse(payload),deparse(output))
  logs <- system2(file.path(R.home("bin"),"Rscript"),c("-e",shQuote(script)),stdout=TRUE,stderr=TRUE)
  expect_null(attr(logs,"status"),info=paste(logs,collapse="\n"))
  fresh <- readRDS(output)
  expect_equal(fresh$prob,predict_probabilities(p["ord4"]),tolerance=0)
  expect_equal(fresh$score,att_expected_score(p["ord4"],c(0,0,0,1)),tolerance=0)
})
