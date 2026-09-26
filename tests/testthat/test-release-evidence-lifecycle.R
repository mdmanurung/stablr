test_that("a filesystem tree record detects added, removed and modified files", {
  skip_on_os("windows")
  env <- .release_module_env()
  root <- withr::local_tempdir(.local_envir = parent.frame())

  dir.create(file.path(root, "nested"), showWarnings = FALSE)
  writeLines("a", file.path(root, "one.txt"))
  writeLines("b", file.path(root, "nested", "two.txt"))

  paths_of <- function(record) vapply(record$manifest, `[[`, character(1L), "path")
  before <- env$.release_tree_record(root)
  expect_true(nzchar(before$sha256))
  expect_setequal(paths_of(before), c("one.txt", "nested", "nested/two.txt"))

  # A modified file must change the record.
  writeLines("a-modified", file.path(root, "one.txt"))
  after_modify <- env$.release_tree_record(root)
  expect_false(identical(before$sha256, after_modify$sha256))

  # A removed file must change it back to something else, not silently match.
  after_remove <- withr::with_dir(root, {
    file.remove(file.path(root, "one.txt"))
    env$.release_tree_record(root)
  })
  expect_false(identical(after_modify$sha256, after_remove$sha256))
  expect_setequal(paths_of(after_remove), c("nested", "nested/two.txt"))

  # An added file must be detected too.
  writeLines("c", file.path(root, "three.txt"))
  expect_length(paths_of(withr::with_dir(root, env$.release_tree_record(root))), 3L)
})

test_that("the tree record ignores nothing, so a stray file cannot hide", {
  skip_on_os("windows")
  env <- .release_module_env()
  root <- withr::local_tempdir(.local_envir = parent.frame())
  writeLines("a", file.path(root, "kept.txt"))
  before <- env$.release_tree_record(root)

  # A dotfile and a build byproduct are both real files on disk; if the record
  # filtered them, a tampered packet could smuggle content past verification.
  writeLines("hidden", file.path(root, ".hidden"))
  writeLines("junk", file.path(root, "stablr-Ex.pdf"))
  after <- env$.release_tree_record(root)

  paths_of <- function(record) vapply(record$manifest, `[[`, character(1L), "path")
  expect_false(identical(before$sha256, after$sha256))
  expect_true(all(c(".hidden", "stablr-Ex.pdf") %in% paths_of(after)))
})

test_that("archive member validation rejects traversal and absolute paths", {
  env <- .release_module_env()
  # A member name that cannot escape the single archive root is accepted; the
  # single-root rule below is what constrains it, not a per-member ban.
  expect_silent(env$.release_validate_archive_members("stablr/DESCRIPTION"))

  for (bad in list(
    "../escape",
    "stablr/../../escape",
    "/absolute/path",
    "C:/windows",
    "stablr/a//b",
    "stablr/./b",
    "stablr/a\\b",
    ""
  )) {
    expect_error(
      env$.release_validate_archive_members(bad),
      class = "stablr_release_identity_error",
      label = paste("must reject:", bad)
    )
  }
  expect_error(
    env$.release_validate_archive_members(character(0)),
    class = "stablr_release_identity_error"
  )
  expect_error(
    env$.release_validate_archive_members(c("stablr/DESCRIPTION", NA_character_)),
    class = "stablr_release_identity_error"
  )
  # A duplicated member would let two archive entries collapse to one path.
  expect_error(
    env$.release_validate_archive_members(c("stablr/DESCRIPTION", "stablr/DESCRIPTION")),
    class = "stablr_release_identity_error"
  )
  # More than one archive root means the candidate is ambiguous.
  expect_error(
    env$.release_validate_archive_members(c("stablr/DESCRIPTION", "other/DESCRIPTION")),
    class = "stablr_release_identity_error"
  )
})

test_that("canonical JSON is byte-stable under key reordering", {
  env <- .release_module_env()
  a <- env$.release_canonical_json(list(b = 1L, a = 2L, c = list(z = 1L, y = 2L)))
  b <- env$.release_canonical_json(list(c = list(y = 2L, z = 1L), a = 2L, b = 1L))
  expect_identical(a, b)
  expect_identical(
    env$.release_canonical_json(list(a = 1L)),
    env$.release_canonical_json(list(a = 1L))
  )
})

test_that("writing a canonical record refuses to replace an existing file", {
  skip_on_os("windows")
  env <- .release_module_env()
  root <- withr::local_tempdir(.local_envir = parent.frame())
  path <- file.path(root, "record.json")

  env$.release_write_canonical_json(list(a = 1L), path)
  expect_true(file.exists(path))

  # No-replace is what makes a sealed record immutable, so a second write must
  # fail loudly rather than overwrite.
  expect_error(
    env$.release_write_canonical_json(list(a = 2L), path),
    class = "stablr_release_identity_error"
  )
  expect_identical(
    env$.release_read_json(path, "record")$a,
    1L
  )
})

test_that("a record cannot be written outside an existing parent directory", {
  env <- .release_module_env()
  root <- withr::local_tempdir(.local_envir = parent.frame())
  expect_error(
    env$.release_write_canonical_json(
      list(a = 1L),
      file.path(root, "absent", "record.json")
    ),
    class = "stablr_release_identity_error"
  )
})

test_that("an explicit packet set must match the contract exactly", {
  env <- .release_module_env()

  verified <- function(promotable = TRUE, blockers = list()) {
    list(
      contract = list(
        required_validations = list(
          list(validation_id = "methodology"),
          list(validation_id = "late_fusion")
        ),
        promotable = promotable,
        blockers = blockers
      )
    )
  }
  packet_ref <- function(id) {
    structure(list(validation_id = id), class = "stablr_packet_ref")
  }
  good <- function(id, status = "passing", promotable = TRUE, sealed = TRUE) {
    list(
      packet = list(validation_id = id, status = status,
                    promotable = promotable, sealed = sealed),
      reference = if (sealed) packet_ref(id) else structure(list(), class = "list"),
      gates = list(list(gate_id = "g", decision_status = "accepted"))
    )
  }
  ok_set <- list(good("methodology"), good("late_fusion"))

  verified <- verified()
  ordered <- env$.release_validate_explicit_packet_set(verified, ok_set)
  expect_type(ordered, "list")
  # The set is returned in CONTRACT order, not in the order supplied.
  expect_identical(
    vapply(ordered, function(p) p$packet$validation_id, character(1L)),
    c("methodology", "late_fusion")
  )
  # Supplying them out of contract order must still normalise.
  expect_identical(
    vapply(
      env$.release_validate_explicit_packet_set(verified, rev(ok_set)),
      function(p) p$packet$validation_id, character(1L)
    ),
    c("methodology", "late_fusion")
  )

  # Too few, too many, and a bare single reference must all fail closed.
  expect_error(
    env$.release_validate_explicit_packet_set(verified, list(good("methodology"))),
    class = "stablr_release_assembly_error"
  )
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(good("methodology"), good("late_fusion"), good("methodology"))
    ),
    class = "stablr_release_assembly_error"
  )
  # A duplicated validation id is not a set.
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(good("methodology"), good("methodology"))
    ),
    class = "stablr_release_assembly_error"
  )
  # A failed packet may not be assembled.
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(good("methodology", status = "failed"), good("late_fusion"))
    ),
    class = "stablr_release_assembly_error"
  )
  # A non-promotable packet may not be assembled.
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(good("methodology", promotable = FALSE), good("late_fusion"))
    ),
    class = "stablr_release_assembly_error"
  )
  # An unsealed packet is a crashed or still-running attempt, not evidence.
  expect_error(
    env$.release_validate_explicit_packet_set(
      verified, list(good("methodology", sealed = FALSE), good("late_fusion"))
    ),
    class = "stablr_release_assembly_error"
  )
})

test_that("the expected contract rejects caller-supplied adapters and settings", {
  # The contract is the authority; a caller must not be able to swap the
  # adapter or weaken the fixed release settings.
  env <- .release_module_env()
  contract <- env$.release_expected_contract()

  for (validation in contract$required_validations) {
    expect_error(
      env$.release_validate_settings(
        validation$settings, contract$required_validations[[1]]$settings,
        validation$validation_id
      ),
      label = paste("settings for", validation$validation_id, "must not validate against another")
    )
  }

  first <- contract$required_validations[[1L]]
  mutated <- .release_contract_copy(first$settings)
  mutated$replicates <- 1L
  expect_error(
    env$.release_validate_settings(mutated, first$settings, first$validation_id),
    class = "stablr_release_contract_error"
  )
})

test_that("a scientific gate that reports pass without acceptance is downgraded", {
  # This is the mechanism that keeps the non-promotable SCI-04 gate honest: a
  # v1 adapter emitting pass = TRUE must not be able to promote a packet.
  env <- .release_module_env()
  verified <- list(
    contract = list(required_validations = list(list(validation_id = "methodology")))
  )
  gates <- data.frame(
    gate_id = "null_selected_fraction",
    gate_version = "SCI-04-null-selected-fraction-v2-draft",
    decision_status = "proposed_unaccepted",
    pass = TRUE,
    stringsAsFactors = FALSE
  )
  expect_error(
    env$.release_validate_gates(
      gates, "methodology", "null_selected_fraction",
      "SCI-04-null-selected-fraction-v2-draft", "proposed_unaccepted"
    ),
    class = "stablr_release_contract_error"
  )
})

test_that("the in-code contract projection matches release-contract.json", {
  # .release_expected_contract() is a hand-maintained R literal that FLATTENS
  # the nested JSON: runtime.dependency_roots and runtime.environment are
  # hoisted, and each required validation contributes its id, adapter, settings,
  # artifact names and gate ids/versions/statuses as parallel vectors. Whole
  # object identity is therefore the wrong assertion; the projection is what
  # must agree, and nothing in the package checked it, so a contract edit could
  # silently diverge from the code that enforces it.
  skip_if_not_installed("jsonlite")
  env <- .release_module_env()
  lit <- env$.release_expected_contract()
  j <- jsonlite::fromJSON(.release_contract_source_path(), simplifyVector = FALSE)
  rv <- j$required_validations

  expect_identical(lit$schema_version, j$schema_version)
  expect_identical(lit$contract_id, j$contract_id)
  expect_identical(lit$hash_algorithm, j$hash_algorithm)
  expect_identical(lit$promotable, j$promotable)
  expect_identical(unlist(lit$blockers), unlist(j$blockers))
  expect_identical(unlist(lit$candidate), unlist(j$candidate))

  expect_setequal(unlist(lit$dependency_roots), unlist(j$runtime$dependency_roots))
  expect_identical(unlist(lit$environment), unlist(j$runtime$environment))

  expect_identical(unlist(lit$validation_ids), unlist(lapply(rv, `[[`, "validation_id")))
  expect_identical(unlist(lit$adapters), unlist(lapply(rv, `[[`, "adapter")))

  for (field in c("settings", "artifacts", "gates")) {
    from_json <- lapply(rv, function(v) v[[field]])
    from_lit <- lit[[field]]
    if (field == "gates") {
      expect_identical(
        unlist(lapply(from_lit, function(g) g$gate_id)),
        unlist(lapply(from_json, function(g) g$gate_id))
      )
      expect_identical(
        unlist(lapply(from_lit, function(g) g$gate_version)),
        unlist(lapply(from_json, function(g) g$gate_version))
      )
      expect_identical(
        unlist(lapply(from_lit, function(g) g$decision_status)),
        unlist(lapply(from_json, function(g) g$decision_status))
      )
    } else if (field == "artifacts") {
      expect_identical(
        unlist(lapply(from_lit, function(a) a$name)),
        unlist(lapply(from_json, function(a) a$name))
      )
    } else {
      # The in-code literal is a simplified projection (character vectors where
      # the JSON holds arrays), so normalise both through the same JSON round
      # trip before comparing rather than comparing representations.
      normalise <- function(x) {
        jsonlite::fromJSON(
          jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"),
          simplifyVector = FALSE
        )
      }
      for (v in seq_along(from_json)) {
        expect_identical(
          normalise(from_lit[[v]]), normalise(from_json[[v]]),
          label = paste(field, names(from_lit)[v])
        )
      }
    }
  }
})

test_that("a changed contract file is detected against the in-code contract", {
  skip_if_not_installed("jsonlite")
  env <- .release_module_env()
  on_disk <- jsonlite::fromJSON(
    .release_contract_source_path(), simplifyVector = FALSE
  )
  mutated <- .release_contract_copy(on_disk)
  mutated$promotable <- !mutated$promotable
  expect_false(env$.release_json_identical(
    jsonlite::toJSON(mutated, auto_unbox = TRUE),
    jsonlite::toJSON(on_disk, auto_unbox = TRUE)
  ))
})
