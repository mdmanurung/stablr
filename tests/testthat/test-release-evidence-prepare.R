test_that("RE-03 candidate archives require one safe root", {
  env <- .release_module_env()

  inspected <- env$.release_validate_archive_members(c(
    "stablr/DESCRIPTION",
    "stablr/inst/release/release-contract.json"
  ))
  expect_identical(inspected$root, "stablr")
  expect_error(
    env$.release_validate_archive_members(c("stablr/DESCRIPTION", "other/x")),
    class = "stablr_release_identity_error"
  )
  expect_error(
    env$.release_validate_archive_members("stablr/../escape"),
    class = "stablr_release_identity_error"
  )
  expect_error(
    env$.release_validate_archive_members(c("stablr/R", "stablr/R/")),
    class = "stablr_release_identity_error"
  )
  expect_error(
    env$.release_validate_archive_members("/stablr/DESCRIPTION"),
    class = "stablr_release_identity_error"
  )
  expect_error(
    env$.release_validate_archive_members("stablr\\escape"),
    class = "stablr_release_identity_error"
  )
})

test_that("RE-03 tree identities detect mutation and links", {
  env <- .release_module_env()
  root <- tempfile("release-tree-")
  dir.create(root)
  writeLines("first", file.path(root, "value.txt"), useBytes = TRUE)

  before <- env$.release_tree_record(root)
  expect_match(before$sha256, "^[0-9a-f]{64}$")
  expect_identical(env$.release_tree_record(root)$sha256, before$sha256)

  writeLines("second", file.path(root, "value.txt"), useBytes = TRUE)
  expect_false(identical(
    env$.release_tree_record(root)$sha256,
    before$sha256
  ))

  source_tree <- env$.release_tree_record(root)
  copied_root <- tempfile("release-tree-copy-")
  dir.create(copied_root)
  expect_true(file.copy(root, copied_root, recursive = TRUE, copy.mode = TRUE))
  copied <- file.path(copied_root, basename(root))
  Sys.chmod(file.path(copied, "value.txt"), mode = "0600")
  env$.release_apply_tree_modes(copied, source_tree$manifest, "copied tree")
  expect_identical(
    env$.release_tree_record(copied)$sha256,
    source_tree$sha256
  )

  link <- file.path(root, "linked.txt")
  if (isTRUE(file.symlink(file.path(root, "value.txt"), link))) {
    expect_error(
      env$.release_tree_record(root),
      class = "stablr_release_identity_error"
    )
  }
})

test_that("RE-03 candidate adoption is content-addressed and no-replace", {
  env <- .release_module_env()
  store <- tempfile("release-store-")
  dir.create(store)
  tarball <- tempfile("candidate-", fileext = ".tar.gz")
  writeBin(charToRaw("candidate bytes"), tarball)
  sha256 <- env$.release_hash_file(tarball)

  expect_true(env$.release_probe_store_atomicity(store))
  stored <- env$.release_store_candidate(tarball, sha256, store)
  expect_identical(env$.release_hash_file(stored), sha256)
  expect_true(dir.exists(file.path(dirname(stored), "SEALED")))
  expect_identical(env$.release_store_candidate(tarball, sha256, store), stored)

  Sys.chmod(stored, mode = "0644")
  writeBin(charToRaw("changed"), stored)
  expect_error(
    env$.release_store_candidate(tarball, sha256, store),
    class = "stablr_release_identity_error"
  )
})

test_that("RE-03 Runtime Artifacts are typed and self-identifying", {
  env <- .release_module_env()
  hash <- paste(rep("a", 64L), collapse = "")
  base <- list(
    schema_version = "stablr.runtime-artifact/v1",
    rscript = list(path = "/runtime/Rscript", sha256 = hash),
    r_executable = list(path = "/runtime/R", sha256 = hash),
    r_home = "/runtime",
    r_version = R.version.string,
    platform = R.version$platform,
    locale = list(
      all = "C", collate = "C", ctype = "C", numeric = "C", time = "C"
    ),
    blas = "/runtime/blas.so",
    lapack = "/runtime/lapack.so",
    library = list(path = "/store/library", tree_sha256 = hash),
    packages = list(list(
      package = "stablr",
      version = "0.1.1",
      source_path = "/store/candidate.tar.gz",
      installed_path = "/store/library/stablr",
      tree_sha256 = hash
    )),
    native_libraries = list(list(path = "/runtime/blas.so", sha256 = hash))
  )
  artifact <- c(list(runtime_id = env$.release_hash_record(base)), base)

  expect_invisible(env$.release_validate_runtime_artifact(artifact))
  artifact$packages[[1L]]$installed_path <- "relative/stablr"
  expect_error(
    env$.release_validate_runtime_artifact(artifact),
    class = "stablr_release_runtime_error"
  )
})

test_that("RE-03 sanitized child environments pin libraries and startup", {
  env <- .release_module_env()
  contract <- env$.release_read_contract(.release_contract_source_path())
  child <- env$.release_sanitized_environment(contract, "/isolated/library")

  expect_true("R_LIBS=/isolated/library" %in% child)
  expect_true("R_LIBS_USER=/isolated/library" %in% child)
  expect_true("R_LIBS_SITE=" %in% child)
  expect_true("R_DEFAULT_PACKAGES=NULL" %in% child)
  expect_true("LC_ALL=C" %in% child)
  expect_true("OMP_NUM_THREADS=1" %in% child)
})

test_that("RE-03 path and preparation contracts fail closed", {
  env <- .release_module_env()

  expect_error(
    env$.release_normalize_absolute("relative", "candidate"),
    class = "stablr_release_identity_error"
  )
  expect_error(
    env$prepare_release_context("relative", "relative", "relative"),
    class = "stablr_release_identity_error"
  )
})
