test_that("signed_series validates and decomposes transitions", {
  s <- small_senate()
  expect_s3_class(s, "signed_series")
  expect_equal(s$n, 30)
  expect_named(s$vertex.attr, c("ID", "Senator", "State", "Party"))
  tn <- transition_networks(s)
  for (t in 1:2) {
    a <- s$y[[t]]; b <- s$y[[t + 1]]
    expect_true(all((tn$F[[t]] != 0) == ((a != 0) | (b != 0))))
    expect_true(all((tn$P[[t]] != 0) == ((a != 0) & (b != 0))))
    expect_true(all((tn$F[[t]] > 0) == ((a > 0) | (b > 0))))
    expect_true(all((tn$P[[t]] > 0) == ((a > 0) & (b > 0))))
  }
  expect_output(print(s), "signed networks")
  expect_output(print(summary(s)), "Signed triangles")
  expect_error(signed_series(list(matrix(0, 3, 3))), "two time points")
  expect_error(signed_series(list(matrix(0, 3, 3), matrix(0, 4, 4))), "same number")
  expect_error(signed_series(list(matrix(0, 3, 3), matrix(0, 3, 3)), vertex.attr = data.frame(a = 1:2)),
               "one row per node")
})
