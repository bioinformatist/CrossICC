context("test-crossicc")

data("demo.platforms")
CrossICC.obj <- CrossICC(demo.platforms, skip.mfs = TRUE, max.iter = 100, use.shiny = FALSE, cross = "cluster", overwrite = TRUE, output.dir = tempdir())
test_that("CrossICC works", {
  expect_equal(length(CrossICC.obj), 10)
})
