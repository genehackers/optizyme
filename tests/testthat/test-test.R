test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that('MCGD example has length 6',{
  expect_vector(example_optimization(),size=6)
})
