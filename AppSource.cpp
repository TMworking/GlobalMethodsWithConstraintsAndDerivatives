#include <iostream>
#include "BaseMethod.h"
#include "TransformMethod.h"
#include "ParabolMethod.h"
#include "BrokenLineMethod.h"
#include "Problems.h"
#include <time.h> 

using namespace std;

enum class TypeOfMethod {
	BrokenLine = 0,
	Parabol = 1,
	Transform = 2
};

shared_ptr<TProblem> ConfigureProblem(int numTask, int constrCount, double mu) {
	//auto problem = make_shared<TProblem>(TypeProblem::TestClass, numTask); // для авто набора задач
	auto problem = make_shared<TProblem>(TypeProblem::AnyTask, numTask); // для ручного набора задач
	cout << endl << "NUM TASK" << numTask << endl;
	problem->setParamsTestTask(TestFuncType::Pinter, TestFuncType::Gibson, constrCount, mu); //1 параметр- Тип целевой функции, 2 параметр- тип функции 1-го ограничения; (остальные - полиномы)
	problem->initializeTask(); // для ручного набора задач
	//problem->generateTestTask(); //Генерация задачи // для авто набора задач
	//problem->seachAlphaValue(); // для авто набора задач
	problem->scanTestTask();
	cout << endl << "minFunc = " << problem->getMin() << " maxFunc = " << problem->getMax() << endl;
	cout << "Число глобальных минимумов = " << problem->getGlobalSolutionCount() << endl;
	for (int i = 0; i < problem->getGlobalSolutionCount(); ++i) //Вывод глобальных минимумов на консоль
		cout << "xMin = " << (*problem)[i].x << " fMin = " << (*problem)[i].y << endl;
	cout << endl;
	return problem;
}

int StartMethod(TypeOfMethod typeOfMethod, shared_ptr<TProblem> problem, bool logNeeded, double eps, int iterationsCount, StopCriterion criterion,
	int N, Gamma gamma1, double v, double d, Gamma gamma2, bool needChangeGamma, int changeGammaStep)
{
	unique_ptr<BaseMethod> method = nullptr;
	string methodName = "";
	switch (typeOfMethod) {
	case TypeOfMethod::BrokenLine:
		method = make_unique<BrokenLineMethod>(problem);
		methodName = "Broken Line Method";
		break;
	case TypeOfMethod::Parabol:
		method = make_unique<ParabolMethod>(problem);
		methodName = "Parabol Method";
		break;
	case TypeOfMethod::Transform:
		method = make_unique<TransformMethod>(problem);
		methodName = "Transform Method";
		break;
	default:
		cout << "Unknown Method" << endl;
		return -1;
	}
	cout << endl << methodName << endl;

	method->setMaxIterationsCount(iterationsCount);
	method->setStopCriterion(criterion);
	method->setParams(N, gamma1, v, d);
	if (needChangeGamma)
		method->setNeedChangeGammaValue(N, gamma2, changeGammaStep);
	method->setLogNeeded(logNeeded);
	method->setEps(eps);
	method->calculate(); // TODO: одновременно два метода выводилось
	int iterarions = method->getIterationsCount();
	cout << "Количество итераций: " << iterarions << "/" << iterationsCount << endl;
	return iterarions;
}

//Создание файла Операционных характеристик
// (ФАЙЛ .TXT, данные в 3 столбца: 1- №задачи, 2- затраченные итерации, 3- maxIterationsCount)
void CreateLogFile(string fileName, vector<int> results_oh, int maxIterationsCount) {
	
	ofstream fout_oh(fileName); //Создаём объект класса ofstream для записи и связываем его с файлом.txt
	bool lineEnd = true;
	for (int elem : results_oh)
	{
		if (lineEnd)
			fout_oh << elem << "\t";
		else
			fout_oh << elem << "\t" << maxIterationsCount << "\n";
		lineEnd = !lineEnd;
	}
	fout_oh.close(); //Закрыли выходной поток в файл
}

int main()
{
	setlocale(LC_ALL, ".1251"); //Русификация

	//=========== Построение тестовой задачи из класса ===========
	int numTask = 1; //Номер генерируемой задачи 
	int constrCount = 4; //Число функций-ограничений
	double mu = 0.15; //Желаемая относительная мера допустимого множества
	//StopCriterion criterion = StopCriterion::GeneratorClose;
	StopCriterion criterion = StopCriterion::CutLength;
	TestFunc::generateFib(); //***
	//=============================================================
	//Параметры для вычисления оценок
	Gamma gamma1(2.3, 1.5, 2.3, 1.5);
	int N = 75; //k* - порог проведенных измерений, когда для завышения глобальной оценки * меняется на +s
	double v = 0.25; //Отступ
	double d = 0.05; //Beta-доля от длины исходного интервала
	double eps = 1.0e-5;
	int maxIterationsCount = 500;
	//Параметры для глобализованной стратегии

	//Беляков
	//double eps = 1.0e-5;
	//int maxIterationsCount = 500;
	//int N = 75; // k*
	//Gamma gamma1(6.0, 1.5, 4.2, 1.5);
	//double d = 2.0e-1; // B-доля интервала
	//double v = 1.0e-1; // отступ

	Gamma gamma2(8.5, 1.2, 5.5, 1.5);
	bool needChangeGamma = true; //Управление балансировкой
	int changeGammaStep = 5; //Параметр балансировки


	//=============================================================
	bool logNeeded = false;

	vector<int> transform_results_oh; //номер задачи и затраченные на решение итерации, последовательное вложение пары данных
	vector<int> parabol_results_oh;
	vector<int> brokenline_results_oh;

	clock_t start = clock(); //Старт таймера
	int successTasksBrokenLine = 0;
	int successTasksParabol = 0;
	int successTasksTransform = 0;

	int StartCount = 123; //С какой задачи начинается тестирование
	int TasksCount = 123; //Размер тестовой выборки
	//for (numTask = StartCount; numTask <= TasksCount; numTask++) {
		//if (StartCount == TasksCount) {
			logNeeded = true;
		//}

		//auto problem = ConfigureProblem(numTask, constrCount, mu); 
		/*int brokenline_iterations_oh = StartMethod(TypeOfMethod::BrokenLine, problem, logNeeded, eps, maxIterationsCount, criterion,
			N, gamma1, v, d, gamma2, needChangeGamma, changeGammaStep);
		if (logNeeded == true) {
			brokenline_results_oh.push_back(numTask);
			brokenline_results_oh.push_back(brokenline_iterations_oh);
		}
			if (brokenline_iterations_oh < maxIterationsCount)
				successTasksBrokenLine++;*/
		for (int i = 2; i < 3; i++) {
			auto problem = ConfigureProblem(i, constrCount, mu);
				int parabol_iterations_oh = StartMethod(TypeOfMethod::Parabol, problem, logNeeded, eps, maxIterationsCount, criterion,
					N, gamma1, v, d, gamma2, needChangeGamma, changeGammaStep);
				if (logNeeded == true) {
					parabol_results_oh.push_back(numTask);
						parabol_results_oh.push_back(parabol_iterations_oh);
				}
			if (parabol_iterations_oh < maxIterationsCount)
				successTasksParabol++;
		}
		

	/*	int transform_iterations_oh = StartMethod(TypeOfMethod::Transform, problem, logNeeded, eps, maxIterationsCount, criterion,
			N, gamma1, v, d, gamma2, needChangeGamma, changeGammaStep);
		if (logNeeded == true) {
			transform_results_oh.push_back(numTask);
			transform_results_oh.push_back(transform_iterations_oh);
		}
			if (transform_iterations_oh < maxIterationsCount)
				successTasksTransform++;*/
		//ОДИНОЧНАЯ ПРОСТАЯ ЗАДАЧА
		/*std::shared_ptr<TProblem> _problem =
			std::shared_ptr<TProblem>(new TProblem(TypeProblem::AnyTask, 1));
		int parabol_iterations_oh = StartMethod(TypeOfMethod::Parabol,
		    _problem, logNeeded, eps, maxIterationsCount, criterion,
			N, gamma1, v, d, gamma2,	needChangeGamma, changeGammaStep);*/
	//}

	clock_t end = clock(); //Завершение работы таймера
	printf("\nЗатраченное время: %f секунд\n", (double)(end - start) / CLOCKS_PER_SEC);
	int numOfTasks = TasksCount - StartCount + 1;
	//cout << "Решено задач BrokenLine Method: " << successTasksBrokenLine << "/" << numOfTasks << ", то есть " << successTasksBrokenLine /((double)numOfTasks / 100) << "%" << endl;
	cout << "Решено задач Parabol Method: " << successTasksParabol << "/" << numOfTasks << ", то есть " << successTasksParabol / ((double)numOfTasks / 100) << "%" << endl;
	//cout << "Решено задач Transform Method: " << successTasksTransform << "/" << numOfTasks << ", то есть " << successTasksTransform / ((double)numOfTasks / 100) << "%" << endl;

	if (logNeeded == true) {
		//CreateLogFile("brokenline_operation_haracteristic.txt", brokenline_results_oh, maxIterationsCount);
		CreateLogFile("parabol_operation_haracteristic.txt", parabol_results_oh, maxIterationsCount);
		//CreateLogFile("transform_operation_haracteristic.txt", transform_results_oh, maxIterationsCount);
	}
	
	cin.get(); //Ожидание завершения работы по клику Enter

	return 0;
}
