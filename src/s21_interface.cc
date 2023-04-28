#include "sfleta_interface.h"

namespace sfleta {
Interface::Interface() {
  dictionary[1] = std::bind(&Interface::AntStart, this);
  dictionary[2] = std::bind(&Interface::GaussStart, this);
  dictionary[3] = std::bind(&Interface::WinogradStart, this);
}

void Interface::Show(int input) {
  while (true) {
    switch (input) {
      case 1 ... 3:
        dictionary[input]();
        WaitingForInput();
        break;
      case 0:
        return;
    }
  }
}

void Interface::WaitingForInput() {
  std::cout << "Для продолжения нажмите Enter\n";
  std::cin.sync();
  std::cin.ignore(INT_MAX, '\n');
  system("clear");
}

void Interface::WrongInputAttention() {
  std::cout << "Неверный ввод, недопустимое значение или введено не число\n";
  std::cin.clear();
  WaitingForInput();
}

void Interface::WinogradStart() {
  Matrix *mat1;
  Matrix *mat2;

  std::stringstream message("");
  message << "Выберите способ ввода:\n"
          << "\t1 - Из файла\n"
          << "\t2 - Случайно\n"
          << "\t0 - Выход\n";
  std::cout << message.str();

  int input = -1;
  std::cin >> input;
  switch (input) {
    case 1:
      try {
        WinogradFromFile(mat1, mat2);
        WaitingForInput();
      } catch (std::exception &except) {
        return;
      }
      break;
    case 2:
      try {
        WinogradFromRandom(mat1, mat2);
        WaitingForInput();
      } catch (std::exception &except) {
        return;
      }
      break;
    case 0:
      std::cout << "bye-bye";
      exit(0);
      break;
    default:
      system("clear");
      return;
  }

  std::cout << "Введите количество выполнений" << std::endl;
  size_t N;
  std::cin >> N;
  if (std::cin.fail()) return;

  std::cout << "Введите количество потоков" << std::endl;
  size_t thread_amount;
  std::cin >> thread_amount;
  if (std::cin.fail()) return;

  try {
    std::unique_ptr<Matrix> res(std::make_unique<sfleta::Matrix>(1, 1));
    std::cout << "Время последовательного выполнения: ";
    {
      Timer timer;
      for (size_t i = 0; i < N; ++i) {
        res.reset(sfleta::Winograd::winograd_algo(*mat1, *mat2));
      }
    }
    std::cout << "Время параллельного выполнения: ";
    {
      Timer timer;
      for (size_t i = 0; i < N; ++i) {
        res.reset(sfleta::Winograd::winograd_algo_multithread(*mat1, *mat2,
                                                           thread_amount));
      }
    }

    std::cout << "Показать матрицы и результат их перемножения?\n"
              << "\t1 - Да\n"
              << "\t2 - Нет\n";

    input = -1;
    std::cin >> input;
    if (std::cin.fail()) return;

    if (input == 1) {
      std::cout << "Первая матрица: \n";
      mat1->Print();
      std::cout << "Вторая матрица: \n";
      mat2->Print();
      std::cout << "Результат умножения матриц: \n";
      res->Print();
    }

  } catch (std::exception &except) {
    std::cout << except.what() << std::endl;
    return;
  }
}

void Interface::WinogradFromFile(Matrix *&mat1_p, Matrix *&mat2_p) {
  try {
    std::string filepath_1;
    std::cout << "Введите путь до файла с первой матрицей" << std::endl;
    std::cin >> filepath_1;
    mat1_p = new Matrix(filepath_1);

    std::string filepath_2;
    std::cout << "Введите путь до файла с второй матрицей" << std::endl;
    std::cin >> filepath_2;
    mat2_p = new Matrix(filepath_2);

  } catch (std::exception &except) {
    std::cout << except.what() << std::endl;
    throw;
  }
}

void Interface::WinogradFromRandom(Matrix *&mat1_p, Matrix *&mat2_p) {
  try {
    int rows1, cols1, rows2, cols2;

    std::cout << "Введите размеры строк и столбцов первой матрицы" << std::endl;
    std::cin >> rows1 >> cols1;

    std::cout << "Введите размеры строк и столбцов второй матрицы" << std::endl;
    std::cin >> rows2 >> cols2;

    mat1_p = new Matrix(rows1, cols1);
    mat2_p = new Matrix(rows2, cols2);
    mat1_p->RandomFill();
    mat2_p->RandomFill();
  } catch (std::exception &except) {
    std::cout << except.what() << std::endl;
    throw;
  }
}

void Interface::GaussStart() {
  std::cout << "Введите путь до файла с матрицей, описывающей СЛАУ"
            << std::endl;
  std::string file_name;
  std::cin >> file_name;

  std::cout << "Введите количество выполнений" << std::endl;
  int N = 0;
  std::cin >> N;
  if (std::cin.fail() || N < 1) {
    return;
  }

  Gauss gauss;

  try {
    gauss.LoadMatrixFromFile(file_name);
  } catch (std::exception &except) {
    std::cout << except.what() << std::endl;
    return;
  }

  std::cout << "Время последовательного выполнения: ";
  {
    Timer timer;
    for (int i = 0; i < N; ++i) gauss.GaussianAlgorithm();
  }
  gauss.PrintResult();

  std::cout << "Время параллельного выполнения: ";
  {
    Timer timer;
    for (int i = 0; i < N; ++i) gauss.GaussianAlgorithm(Gauss::launch::ASYNC);
  }
  gauss.PrintResult();

  std::stringstream message("");
  message << "Повторим?\n"
          << "\t1 - Погнали\n"
          << "\t0 - Сдаюсь\n";
  std::cout << message.str();

  int input = -1;
  std::cin >> input;
  switch (input) {
    case 1:
      return GaussStart();
    case 0:
      exit(0);
    default:
      system("clear");
      return;
  }
}

void Interface::AntStart() {
  system("clear");

  Graph graph;
  AntAlgorithm ant;

  int choice = -1;
  int N = 1;
  std::string file_name;
  bool mthread = false;

  sfleta::TsmResult result;

  for (;;) {
    PrintAntMenu();
    std::cin >> choice;
    if (std::cin.fail() || (choice != 1 && choice != 2 && choice != 0)) {
      WrongInputAttention();
    } else {
      switch (choice) {
        case 1:
          system("clear");
          std::cout << "Введите путь до файла и его имя\n";
          std::cin >> file_name;
          try {
            graph.LoadGraphFromFile(file_name);
          } catch (const std::exception &e) {
            std::cerr << e.what() << std::endl;
          }
          break;
        case 2:
          system("clear");
          std::cout << "Задайте количество выполнения алогритма не меньше одно "
                       "раза\n";
          std::cin >> N;
          if (std::cin.fail() || N < 1) {
            WrongInputAttention();
          } else {
            if (graph.GetSize()) {
              std::cout << "Вычисление муравьиного алгоритма без применения "
                           "параллельных вычислений\n";
              mthread = false;
              try {
                {
                  sfleta::Timer timer;
                  for (int i = 0; i < N; i++) {
                    result = ant.SolveTravelingSalesmanProblem(graph, mthread);
                  }
                }
                ant.PrintResult(result);
              } catch (const std::exception &e) {
                std::cerr << e.what() << '\n';
              }

              std::cout << "Вычисление муравьиного алгоритма с применением "
                           "параллельных вычислений\n";
              mthread = true;
              try {
                {
                  sfleta::Timer timer;
                  for (int i = 0; i < N; i++) {
                    result = ant.SolveTravelingSalesmanProblem(graph, mthread);
                  }
                }
                ant.PrintResult(result);
              } catch (const std::exception &e) {
                std::cerr << e.what() << '\n';
              }
              N = 1;
            } else {
              std::cout << "Вы не задали матрицу для муравьиного алгоритма\n";
            }
          }
          break;
        case 0:
          exit(0);
      }
    }
  }
}

void Interface::PrintAntMenu() {
  std::cout << "Выберите действие:\n";
  std::cout << "\t1 - Загрузка нового графа из файла\n";
  std::cout << "\t2 - Вычисление муравьиного алгоритма\n";
  std::cout << "\t0 - Выход\n";
}
}  // namespace sfleta
