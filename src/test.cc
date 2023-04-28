#include "ant/sfleta_ant_algorithm.h"
#include "gauss/sfleta_gauss.h"
#include "gtest/gtest.h"
#include "winograd/sfleta_winograd_algorithm.h"

namespace sfleta {
std::unique_ptr<Matrix> res(std::make_unique<Matrix>(1, 1));

TEST(fill_from_file, identical_matrices) {
  sfleta::Matrix mat1("./winograd/examples/matrix_5x5.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_5x5.txt");

  ASSERT_TRUE(mat1.getBuffer() == mat2.getBuffer());
}

TEST(fill_from_file, not_identical_matrices) {
  sfleta::Matrix mat1("./winograd/examples/matrix_5x5.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_5x5_2.txt");

  ASSERT_FALSE(mat1.getBuffer() == mat2.getBuffer());
}

TEST(winograd, test_1) {
  sfleta::Matrix mat1("./winograd/examples/matrix_1x5.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_5x1.txt");

  sfleta::Matrix ref("./winograd/examples/ref_1x5_mul_5x1.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_2) {
  sfleta::Matrix mat1("./winograd/examples/matrix_5x1.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_1x5.txt");

  sfleta::Matrix ref("./winograd/examples/ref_5x1_mul_1x5.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_3) {
  sfleta::Matrix mat1("./winograd/examples/matrix_5x2.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_2x5.txt");

  sfleta::Matrix ref("./winograd/examples/ref_5x2_mul_2x5.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_4) {
  sfleta::Matrix mat1("./winograd/examples/matrix_2x5.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_5x2.txt");

  sfleta::Matrix ref("./winograd/examples/ref_2x5_mul_5x2.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_5) {
  sfleta::Matrix mat1("./winograd/examples/matrix_5x5.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_5x5_2.txt");

  sfleta::Matrix ref("./winograd/examples/ref_5x5_mul_5x5_2.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_6) {
  sfleta::Matrix mat1("./winograd/examples/matrix_5x5_2.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_5x5.txt");

  sfleta::Matrix ref("./winograd/examples/ref_5x5_2_mul_5x5.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_7) {
  sfleta::Matrix mat1("./winograd/examples/matrix_5x10.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_10x5.txt");

  sfleta::Matrix ref("./winograd/examples/ref_5x10_mul_10x5.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_8) {
  sfleta::Matrix mat1("./winograd/examples/matrix_10x5.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_5x10.txt");

  sfleta::Matrix ref("./winograd/examples/ref_10x5_mul_5x10.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_9) {
  sfleta::Matrix mat1("./winograd/examples/matrix_10x10.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_10x10_2.txt");

  sfleta::Matrix ref("./winograd/examples/ref_10x10_mul_10x10_2.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_10) {
  sfleta::Matrix mat1("./winograd/examples/matrix_10x10_2.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_10x10.txt");

  sfleta::Matrix ref("./winograd/examples/ref_10x10_2_mul_10x10.txt");
  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));

  ASSERT_TRUE(ref.getBuffer() == res->getBuffer());
}

TEST(winograd, test_11) {
  sfleta::Matrix mat1("./winograd/examples/matrix_33x49.txt");
  sfleta::Matrix mat2("./winograd/examples/matrix_49x33.txt");

  res.reset(sfleta::Winograd::winograd_algo_multithread(mat1, mat2));
  mat1.mul_matrix(mat2);

  ASSERT_TRUE(mat1.getBuffer() == res->getBuffer());
}

TEST(winograd, test_12) {
  sfleta::Matrix mat1(5, 5);
  sfleta::Matrix mat2(5, 5);

  mat1.RandomFill();
  mat2.RandomFill();

  res.reset(sfleta::Winograd::winograd_algo(mat1, mat2));
  mat1.mul_matrix(mat2);

  for (size_t i = 0; i < mat1.getRows(); ++i) {
    for (size_t j = 0; j < mat1.getCols(); ++j) {
      ASSERT_NEAR(mat1.getBuffer()[i][j], res->getBuffer()[i][j], sfleta::kEps);
    }
  }
}

//// THROW_TESTS

//    TEST(test_throw, wrong_path) {
//        ASSERT_ANY_THROW(sfleta::Matrix mat1("./examplus/matrix_5x5.txt"));
//    }
//
//    TEST(test_throw, wrong_size) {
//        ASSERT_ANY_THROW(sfleta::Matrix mat1(0, 0));
//    }

TEST(Gauss, test_1) {
  std::vector checked = {6.102793,  63.749497, 5.676335,
                         -8.476871, -5.542569, -0.6};
  sfleta::Gauss gauss;
  gauss.LoadMatrixFromFile("./gauss/examples/matrix6x6.txt");
  gauss.GaussianAlgorithm();
  auto result = gauss.GetResult();
  ASSERT_EQ(checked.size(), result.size());
  for (size_t i = 0; i < result.size(); ++i)
    ASSERT_NEAR(result[i], checked[i], sfleta::kEps);
}

TEST(Gauss, test_2) {
  std::vector checked = {-7.444263, 2.543527, -4.028256, 3.217904, 9.321953};
  sfleta::Gauss gauss;
  gauss.LoadMatrixFromFile("./gauss/examples/matrix5x5.txt");
  gauss.GaussianAlgorithm();
  auto result = gauss.GetResult();
  ASSERT_EQ(checked.size(), result.size());
  for (size_t i = 0; i < result.size(); ++i)
    ASSERT_NEAR(result[i], checked[i], sfleta::kEps);
}

TEST(Gauss, test_3_async) {
  sfleta::Gauss gauss;
  gauss.LoadMatrixFromFile("./gauss/examples/matrix5x5_null.txt");
  gauss.GaussianAlgorithm(sfleta::Gauss::launch::ASYNC);
  auto result = gauss.GetResult();
  ASSERT_TRUE(result.empty());
}

TEST(Gauss, test_1_async) {
  std::vector checked = {6.102793,  63.749497, 5.676335,
                         -8.476871, -5.542569, -0.6};
  sfleta::Gauss gauss;
  gauss.LoadMatrixFromFile("./gauss/examples/matrix6x6.txt");
  gauss.GaussianAlgorithm(sfleta::Gauss::launch::ASYNC);
  auto result = gauss.GetResult();
  ASSERT_EQ(checked.size(), result.size());
  for (size_t i = 0; i < result.size(); ++i)
    ASSERT_NEAR(result[i], checked[i], sfleta::kEps);
}

TEST(Gauss, test_2_async) {
  std::vector checked = {-7.444263, 2.543527, -4.028256, 3.217904, 9.321953};
  sfleta::Gauss gauss;
  gauss.LoadMatrixFromFile("./gauss/examples/matrix5x5.txt");
  gauss.GaussianAlgorithm(sfleta::Gauss::launch::ASYNC);
  auto result = gauss.GetResult();
  ASSERT_EQ(checked.size(), result.size());
  for (size_t i = 0; i < result.size(); ++i)
    ASSERT_NEAR(result[i], checked[i], sfleta::kEps);
}

TEST(Gauss, test_2_asunc) {
  std::vector checked = {-7.444263, 2.543527, -4.028256, 3.217904, 9.321953};
  sfleta::Gauss gauss;
  gauss.LoadMatrixFromFile("./gauss/examples/matrix5x5.txt");
  gauss.GaussianAlgorithm(sfleta::Gauss::launch::ASYNC);
  auto result = gauss.GetResult();
  for (size_t i = 0; i < result.size(); ++i)
    ASSERT_NEAR(result[i], checked[i], sfleta::kEps);
}

TEST(SolveTravelingSalesmanProblem, test_1) {
  sfleta::Graph checking_graph;
  checking_graph.LoadGraphFromFile("./ant/examples/sale_1.txt");

  bool mthread = false;

  sfleta::AntAlgorithm alg;
  const sfleta::TsmResult res =
      alg.SolveTravelingSalesmanProblem(checking_graph, mthread);

  double answer = 253;
  ASSERT_EQ(answer, res.distance);
}

TEST(SolveTravelingSalesmanProblem, test_2) {
  sfleta::Graph checking_graph;
  checking_graph.LoadGraphFromFile("./ant/examples/sale_1.txt");

  bool mthread = true;

  sfleta::AntAlgorithm alg;
  const sfleta::TsmResult res =
      alg.SolveTravelingSalesmanProblem(checking_graph, mthread);

  double answer = 253;
  ASSERT_EQ(answer, res.distance);
}

TEST(SolveTravelingSalesmanProblem, test_3) {
  sfleta::Graph checking_graph;
  checking_graph.LoadGraphFromFile("./ant/examples/sale_2.txt");

  bool mthread = false;

  sfleta::AntAlgorithm alg;
  const sfleta::TsmResult res =
      alg.SolveTravelingSalesmanProblem(checking_graph, mthread);

  double answer = 127;
  ASSERT_EQ(answer, res.distance);
}

TEST(SolveTravelingSalesmanProblem, test_4) {
  sfleta::Graph checking_graph;
  checking_graph.LoadGraphFromFile("./ant/examples/sale_2.txt");

  bool mthread = true;

  sfleta::AntAlgorithm alg;
  const sfleta::TsmResult res =
      alg.SolveTravelingSalesmanProblem(checking_graph, mthread);

  double answer = 127;
  ASSERT_EQ(answer, res.distance);
}

TEST(SolveTravelingSalesmanProblem, test_5) {
  sfleta::Graph checking_graph;
  checking_graph.LoadGraphFromFile("./ant/examples/sale_3.txt");

  bool mthread = false;

  sfleta::AntAlgorithm alg;
  const sfleta::TsmResult res =
      alg.SolveTravelingSalesmanProblem(checking_graph, mthread);

  double answer = 69;
  ASSERT_EQ(answer, res.distance);
}

TEST(SolveTravelingSalesmanProblem, test_6) {
  sfleta::Graph checking_graph;
  checking_graph.LoadGraphFromFile("./ant/examples/sale_3.txt");

  bool mthread = true;

  sfleta::AntAlgorithm alg;
  const sfleta::TsmResult res =
      alg.SolveTravelingSalesmanProblem(checking_graph, mthread);

  double answer = 69;
  ASSERT_EQ(answer, res.distance);
}

//// THROW_TESTS

//    TEST(SolveTravelingSalesmanProblem, test_7) {
//        sfleta::Graph checking_graph;
//        checking_graph.LoadGraphFromFile("./ant/examples/sale_4_err.txt");
//
//        bool mthread = false;
//
//        sfleta::AntAlgorithm alg;
//        EXPECT_THROW(alg.SolveTravelingSalesmanProblem(checking_graph,
//        mthread), std::exception);
//    }
//
//    TEST(SolveTravelingSalesmanProblem, test_8) {
//        sfleta::Graph checking_graph;
//        checking_graph.LoadGraphFromFile("./ant/examples/sale_4_err.txt");
//
//        bool mthread = true;
//
//        sfleta::AntAlgorithm alg;
//        EXPECT_THROW(alg.SolveTravelingSalesmanProblem(checking_graph,
//        mthread), std::exception);
//    }
//
//    TEST(SolveTravelingSalesmanProblem, test_9) {
//        sfleta::Graph checking_graph;
//        checking_graph.LoadGraphFromFile("./ant/examples/sale_5_err.txt");
//
//        bool mthread = false;
//
//        sfleta::AntAlgorithm alg;
//        EXPECT_THROW(alg.SolveTravelingSalesmanProblem(checking_graph,
//        mthread), std::exception);
//    }
//
//    TEST(SolveTravelingSalesmanProblem, test_10) {
//        sfleta::Graph checking_graph;
//        checking_graph.LoadGraphFromFile("./ant/examples/sale_5_err.txt");
//
//        bool mthread = true;
//
//        sfleta::AntAlgorithm alg;
//        EXPECT_THROW(alg.SolveTravelingSalesmanProblem(checking_graph,
//        mthread), std::exception);
//    }
//
//    TEST(SolveTravelingSalesmanProblem, test_11) {
//        sfleta::Graph checking_graph;
//        EXPECT_THROW(checking_graph.LoadGraphFromFile("./ant/examples/matrix_err_1.txt"),
//        std::exception);
//    }
//
//    TEST(SolveTravelingSalesmanProblem, test_12) {
//        sfleta::Graph checking_graph;
//        EXPECT_THROW(checking_graph.LoadGraphFromFile("./ant/examples/matrix_err_2.txt"),
//        std::exception);
//    }
//
//    TEST(SolveTravelingSalesmanProblem, test_13) {
//        sfleta::Graph checking_graph;
//        EXPECT_THROW(checking_graph.LoadGraphFromFile("./ant/examples/no_file.txt"),
//        std::exception);
//    }
}  // namespace sfleta

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
