#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

class Symplex_method {
private:
	double * c; // строка коэффициентов функции
	std::size_t length;

	double ** A; // матрица ограничений
	std::size_t rows;
	std::size_t columns;

	double * b; // стобец ограничений
	std::size_t high;

	double * x; // столбец переменных х
	std::size_t number;

	double * basis; // базисные переменные
	double * free; // свободные переменные

	double ** symplex_table; //симплекс-таблица
	std::size_t rows_;
	std::size_t columns_;

	void creating_symplex_table() {
			rows_ = high + 1; //симплекс-таблица имеет размер на одну больше строку (значения функции)
			columns_ = length + 1; // и на один больше столбец (столбец свободных членов), чем матрица А
			symplex_table = new double * [rows_]; // выделение памяти для таблицы
			for (auto i = 0; i < rows_; i++) {
				symplex_table[i] = new double[columns_];
				for (auto j = 0; j < columns_; j++) {
					symplex_table[i][j] = 0;
				}
			}

			for (auto j = 0; j < high; j++) { // заполнение столбца свободных членов
				symplex_table[j][0] = b[j];
			}

			symplex_table[rows_ - 1][0] = 0;

			for (auto i = 1; i <= length; i++) { // заполнение строки функции
				symplex_table[rows_ - 1][i] = c[i-1];
			}

			for (auto i = 0; i < rows_ - 1; i++) { // заполнение таблицы начальной матрицей
				for (auto j = 1; j < columns_; j++) {
					symplex_table[i][j] = A[i][j - 1];
				}
			}
		}

	void jordan_exceptions(int i, int k) {
		double ** between = new double*[rows_]; //создание вспомогательной симплекс-таблицы,

		for (auto m = 0; m < rows_; m++) { // хранящей значения предыдущей
			between[m] = new double[columns_];
			for (auto n = 0; n < columns_; n++) {
				between[m][n] = symplex_table[m][n];
			}
		}

		symplex_table[i][k] = 1 / between[i][k]; //заполнение симплекс-таблицы новыми значениями

		// Далее используются по два отдельных цикла (до разрещающей строки/столбца
		// и после него). Сделано для быстродействия программы, чтобы в каждой
		// итерации не проверялось условие if(j==K)

		for (auto j = 0; j < k; j++) {
			symplex_table[i][j] = between[i][j] / between[i][k];
		}
		for (auto j = k + 1; j < columns_; j++) {
			symplex_table[i][j] = between[i][j] / between[i][k];
		}

		for (auto j = 0; j < i; j++) {
			symplex_table[j][k] = -between[j][k] / between[i][k];
		}
		for (auto j = i + 1; j < rows_; j++) {
			symplex_table[j][k] = -between[j][k] / between[i][k];
		}


		for (auto m = 0; m < rows_; m++) {
			if (m == i) {
				continue;
			}
			for (auto n = 0; n < columns_; n++) {
				if (n == k) {
					continue;
				}
				symplex_table[m][n] = between[m][n] - (between[m][k] * between[i][n] / between[i][k]);
			}
		}

		for (auto j = 0; j < rows_; j++) { //очищение памяти
			delete[] between[j];
		}
		delete[] between;

		for (auto j = 0; j < high; j++) { // занесение значений переменных в столбец x
			x[int(basis[j]) - 1] = symplex_table[j][0];
		}
		for (auto j = 0; j < length; j++) {
			x[int(free[j]) - 1] = 0;
		}

		//std::swap(x[high + i], x[k-1]); // с учетом смены строки-столбца
		std::swap(x[int(basis[i]) - 1], x[int(free[k - 1]) - 1]);
		std::swap(basis[i], free[k-1]);

	}

	void finding_reference_solution() {
		// Поиск опорного решения согласно известному алгоритму
		int i = -1;
		int k = -1;

		for (auto j = 0; j < high; j++) {
			if (symplex_table[j][0] < 0) {
				i = j;

				for (auto n = 1; n < columns_; n++) { // поиск разрешающего столбца
					if (symplex_table[i][n] < 0) {
						k = n;
						break;
					}
				}
				if (k != -1) {
					break;
				}
			}
		}

		if (i == -1) {
			for (auto j = 0; j < high; j++) {
				x[j + length] = symplex_table[j][0];
			}
			return; // если столбец свободных членов неотрицательный, это и есть опорное решение
		}

		if (i != -1 && k == -1) {
			std::cout << "No solution"; // если ни в одной строке не имеется отрицательных членов, решения не существует
			exit(0);
		}

		double min = symplex_table[i][0] / symplex_table[i][k]; // поиск разрешающей строки
		for (auto j = 0; j < rows_ - 1; j++) {
			double Min_ = symplex_table[j][0] / symplex_table[j][k];
			if (Min_ < min && Min_ > 0) {
				min = Min_;
				i = j;
			}
		}
		write_table();
		jordan_exceptions(i, k); // изменение симплекс-таблицы
		finding_reference_solution();
	}

	bool is_minimum() {
		// проверка, возможно ли более минимальное значение функции
		bool succsess = true;

		for (auto i = 1; i < columns_; i++) {
			if (symplex_table[rows_ - 1][i] > 0) {
				succsess = false;
				break;
			}
		}

		return succsess;
	}

	void find_resolving_column(int & k) {
		// поиск разрешающего столбца
		for (auto i = k + 1; i < columns_; i++) {
			if (symplex_table[rows_ - 1][i] > 0) {
				k = i;
				return;
			}
		}
		k = columns_ - 1;
	}

	auto find_resolving_row(int & k) {
		// поиск разрешающей строки
		unsigned int i = 0;
		double min;
		for (auto j = 0; j < rows_ - 1; j++) {
			if (symplex_table[j][k] > 0 && symplex_table[rows_ - 1][k] > 0) {
				if (symplex_table[j][0] >= 0) {
					min = symplex_table[j][0] / symplex_table[j][k];
					i = j;
					break;
				}
			}
			else if (j == rows_ - 2) {
				if (k != columns_ - 1) {
					find_resolving_column(k);
					j = -1;
				}
				else {
					std::cout << "No maximum"; // если разрешающий столбец отрицателен, то функция неограничена
					exit(0);
				}
			}
		}

		for (auto j = i; j < rows_ - 1; j++) {
			if (symplex_table[j][k] != 0) {
				double Min_ = symplex_table[j][0] / symplex_table[j][k];
				if ( Min_ >= 0 && Min_ < min) {
					min = Min_;
					i = j;
				}
			}
		}

		return i;
	}

	void write_table() {
		// функция выводит на печать форматированную симплекс-таблицу, значение столбца
		// переменных х с учетом мнимых переменных (опорное решение)
		// и значение функции на данной итерации
		std::cout << std::endl << std::setw(61) << std::setfill('-') << ' ' << std::endl << std::setfill(' ');

		for (auto i = 0; i < columns_ + 1; i++) {
			if (i == 0) {
				std::cout << std::setw(12) << '|';
			}
			else {
				if (i == 1) {
					std::cout << std::setw(12) << "Si0|";
				}
				else {
					std::cout << std::setw(10) << 'X' << free[i - 2] << '|';
				}
			}
		}

		std::cout << std::endl << std::setw(61) << std::setfill('-') << ' ' << std::endl << std::setfill(' ');

		for (auto i = 0; i < rows_ - 1; i++) {
			std::cout << std::setw(10) << 'X' << basis[i] << '|';
			for (auto j = 0; j < columns_; j++) {
				std::cout << std::setw(10) << std::setprecision(3) << symplex_table[i][j] << " |";
			}
			std::cout << std::endl << std::setw(61) << std::setfill('-') << ' ' << std::endl << std::setfill(' ');
		}

		for (auto i = 0; i < columns_ + 1; i++) {
			if (i == 0) {
				std::cout << std::setw(12) << "F|";
			}
			else {
				std::cout << std::setw(10) << std::setprecision(3) << symplex_table[rows_ - 1][i - 1] << " |";
			}
		}

		std::cout << std::endl << std::setw(61) << std::setfill('-') << ' ' << std::endl << std::setfill(' ');

		std::cout << std::endl << "Referance solution is:\nX = (" ;

		for (auto j = 0; j < number; j++) {
			if (j != number - 1) {
				std::cout << std::setw(5) << std::setprecision(3) << x[j] << ',';
			}
			else {
				std::cout << std::setw(5) << std::setprecision(3) << x[j] << ')';
			}
		}

		std::cout << std::endl << std::endl << "F maximum is  " << std::setprecision(3) << -symplex_table[rows_ - 1][0] << std::endl << std::endl;
	}

	void write_solution() {
		// функция выводит на печать оптимальное решение
		// и максимальное значение функции
		std::cout << std::endl << "Optimal solution is:\n";

		for (auto j = 0; j < length; j++) {
			for (auto i = 0; i < high; i++) {
				if (basis[i] == j + 1) {
					std::cout << 'X' << j + 1 << " = " << std::setprecision(3) << x[j] << std::endl;
					break;
				}
				if (i == high - 1) {
					std::cout << 'X' << j + 1 << " = " << std::setprecision(3) << 0 << std::endl;
				}		
			}
		}

		std::cout << "F maximum is  " << std::setprecision(3) << -symplex_table[rows_ - 1][0] << std::endl << std::endl;
	}

	double create_new_problem(std::vector<double> & vector, int k, int m = 1) {
		std::ofstream fout(R"(C:\Users\Asus\Desktop\another.txt)");
		fout << length << std::endl;
		for (auto i = 0; i < length; ++i) {
			fout << c[i];
			if (i != high - 1)
				fout << ' ';
		}
		fout << std::endl << ';' << std::endl;

		fout << high + 1 << ", " << length << std::endl;
		for (auto i = 0; i < high; ++i) {
			for (auto j = 0; j < length; ++j) {
				fout << A[i][j];
				if (j != length - 1)
					fout << ' ';
			}
			fout << std::endl;
		}
		fout << 1 * m << ' ';
		for (auto j = 0; j < length - 1; ++j) {
			fout << 0;
			if (j != length - 2)
				fout << ' ';
		}

		fout << std::endl << ';' << std::endl;
		fout << high + 1 << std::endl;
		for (auto i = 0; i < high; ++i) {
			fout << b[i] << ' ';
		}
		fout << m*k;
		fout.close();

		Symplex_method another_problem(R"(C:\Users\Asus\Desktop\another.txt)");

		return another_problem.MBB(vector);

		return 0;
	};

	double MBB(std::vector<double> & vector) {
		algorithm();
		for (auto i = 0; i < length; ++i) {
			vector.push_back(x[i]);
		}

		double record = -symplex_table[rows_ - 1][0], record1, record2;
		int floor, ceil;
		std::vector<double> vector1, vector2;
		for (auto i = 0; i < length; ++i) {
			if (int(x[i]) != x[i]) {
				floor = int(x[i]);
				ceil = int(x[i] + 1);

				record1 = create_new_problem(vector1, floor);

				record2 = create_new_problem(vector2, ceil, -1);

				record1 = record1 < record2 ? record2 : record1;
				if (record1 == record2)
					vector1 = vector2;

				record = record1 < record ? record1 : record;
				if (record == record1)
					vector = vector1;
			}
		}
		return record;
	}

public:
	void dual_problem() {
		std::swap(b, c);
		std::swap(length, high);
		for (auto i = 0; i < high; ++i) {
			b[i] *= -1;
		}
		for (auto i = 0; i < length; ++i) {
			c[i] *= -1;
		}

		double ** B = new double * [columns];

		for (auto i = 0; i < columns; i++) {
			B[i] = new double[rows];
		}

		for (auto i = 0; i < columns; i++) {
			for (auto j = 0; j < rows; ++j) {
				B[i][j] = A[j][i] * -1;
			}
		}

		std::swap(A, B);

		for (auto i = 0; i < rows; i++) {
			delete[] B[i];
		}
		delete[] B;

		for (auto i = 0; i < rows_; i++) {
			delete[] symplex_table[i];
		}
		delete[] symplex_table;

		std::swap(rows, columns);

		for (auto i = 0; i < number; i++) {
			x[i] = 0;
		}
		delete[] free;
		delete[] basis;

		free = new double [length];
		basis = new double [high];

		for (auto i = 0; i < length; i++) {
			free[i] = i + 1;
		}

		for (auto i = 0; i < high; i++) {
			basis[i] = length + i + 1;
		}
		
		creating_symplex_table();
		write_table();

		algorithm();
	}

	void MBB() {
		std::vector<double> vector;
		double record = MBB(vector);

		std::cout << "Optimal integer solution is:" << std::endl;

		for (auto i = 0; i < length; i++) {
			std::cout << "X" << i + 1 << " = " << vector[i] << std::endl;
		}

		std::cout << "F(integer) maximum is " << record;
	}

	Symplex_method(std::string name_of_file = R"(C:\Users\Asus\Desktop\lab01.txt)") {

		std::ifstream fin(name_of_file); // чтение данных из файла
		if ( fin.is_open() ) {
			char symbol;

			if (fin >> length) { // считывание строки коэффициентов функции с
				c = new double[length];
				for (auto i = 0; i < length; i++) {
					fin >> c[i];
				}
			}

			if (fin >> symbol && symbol == ';' &&
				fin >> rows && fin >> symbol &&
				symbol == ',' && fin >> columns) { // считывание матрицы ограничений А
				A = new double * [rows];
				for (auto i = 0; i < rows; i++) {
					A[i] = new double[columns];
					for (auto j = 0; j < columns; j++) {
						fin >> A[i][j];
					}
				}
			}
			
			if (fin >> symbol && symbol == ';' && fin >> high) { // считывание столбца ограничений b
				b = new double[high];
				for (auto i = 0; i < high; i++) {
					fin >> b[i];
				}
			}

			number = length + high;
			x = new double [number]; // создание столбца переменных с учетом мнимых переменных
			for (auto i = 0; i < number; i++) {
				x[i] = 0;
			}
			
			free = new double[length]; // создание массива свободных переменных
			for (auto j = 0; j < length; j++) {
				free[j] = j + 1;
			}

			basis = new double[high]; // создание базисных свободных переменных
			for (auto j = 0; j < high; j++) {
				basis[j] = j + length + 1;
			}

			creating_symplex_table(); // создание симплекс таблицы
			fin.close();
		}
	}

	~Symplex_method() {
		//очистка памяти от занимаемых ресурсов
		delete [] c;
		delete [] b;
		for (auto i = 0; i < rows; i++) {
			delete[] A[i];
		}
		delete[] A;

		delete[] x;
		delete[] free;
		delete[] basis;

		for (auto i = 0; i < rows_; i++) {
			delete[] symplex_table[i];
		}
		delete[] symplex_table;
	}

	void algorithm() {
		finding_reference_solution(); // нахождение опорного решения

		write_table(); // вывод таблицы на печать

		while (!is_minimum()) {

			int k = 0;
			find_resolving_column(k); // поиск разрешающего столбца
			auto i = find_resolving_row(k); // поиск разрешающей строки
			
			std::cout << "Resolving colomn is " << k << "\tResolving row is " << i << std::endl;
			jordan_exceptions(i, k); // изменение симплекс таблицы
			
			write_table(); // вывод таблицы на печать
		}

		write_solution(); // вывод на печать оптимального решения задачи
	}
};

int main()
{
	Symplex_method lab; // создание объекста класса
	//lab.dual_problem(); // решение двойственной задачи лабараторной работы
	//lab.algorithm();	// решение задачи лабараторной работы
	lab.MBB();

	return 0;
}
