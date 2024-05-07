// +++ ВЕРСИЯ ГЕНЕРАТОРА ТЕСТОВЫХ ЗАДАЧ УСЛОВНОЙ ГЛОБАЛЬНОЙ ОПТИМИЗАЦИИ: v2 от 12.07.2019
#pragma once
#include <vector>
#include <iostream>
#include <corecrt_math_defines.h>
#include <algorithm>
using namespace std;

enum class TypeProblem {AnyTask, TestClass};
enum class TestFuncType {Pinter,Sheckel,Gibson,Polinom};

using TGlobalSolution = struct Point { double x; double y; };

class TProblem;

class TestFunc
{
 protected:
	double Fa, Fb;//Принятый диапазон определения данной тестовой функции
	double FaNew, FbNew;//Диапазон, к которому функция должна быть приведена
	int pointsCount{ 1000 };//Число точек в сетке при грубом определении Min и Max
	double epsX{ 1.0e-08 }; //Точность определения глобальных минимумов по х
	double epsY{ 1.0e-05 }; //Точность сравнения значений в глоб.минимумах
	vector<TGlobalSolution> globalSolutions; //Набор глобальных минимумов тестовой задачи
	double minFunc, maxFunc;//Границы колебания функции
	//Временное хранилище грубых оценок локальных минимумов тестовой функции:
	static vector<TGlobalSolution> localTmpSolutions;
	static vector<unsigned int> Fib;//Числа Фибоначчи
 public:
	TestFunc(double aNew, double bNew, double a, double b) 
		: Fa(a), Fb(b), FaNew(aNew), FbNew(bNew)
	{
		generateFib();
		/*if (!Fib.size())
		{ Fib.push_back(1); Fib.push_back(1);
		for (int i = 2; i < 46; ++i)
			Fib.push_back(Fib[i - 2] + Fib[i - 1]);
		}*/
	}
 virtual void paramGeneration() = 0; //Генерация параметров функции
         //Вычисление функции в "родных" координатах из [Fa, Fb]:
 virtual double calcTestFunc(double x) = 0;
    //Оператор вычисления функции в координатах из [FaNew, FbNew]:
	double operator()(double x)
	{	double xx = Fa + (Fb - Fa)*(x - FaNew) / (FbNew - FaNew);
		return calcTestFunc(xx);
	}
	    //Вычисл-е произ-ой функ-и в "родных" координатах из [Fa, Fb]:
virtual double calcDerivativeTestFunc(double x) = 0;
    //Оператор выч-ия произв-й функции в коорд-ах из [FaNew, FbNew]:
	double operator()(double x, char t)
	{ if (t != 'd') return 0.0;
	   double K = (Fb - Fa) / (FbNew - FaNew);
	   double xx = Fa + (Fb - Fa)*(x - FaNew) / (FbNew - FaNew);
	   return calcDerivativeTestFunc(xx)*K;
	}
	double getBaseLeft() { return Fa; }
	double getBaseRight() { return Fb; }
	 //Определение максимального и минимального значений, 
	// формирование списк точно найденных глобальных минимумов:
    void scanTestFunction();
	double getMin() { return minFunc; }
	double getMax() { return maxFunc; }
	int globalSolutionsCount() { return (int) globalSolutions.size(); }
	TGlobalSolution operator[](int i)
	{	i = i < 0 ? 0 : i;
		i = (int)((size_t)i >= globalSolutions.size() ? globalSolutions.size() - 1 : i);
		return globalSolutions[i];
	}
//Поиск минимума на [x0, x1] из [0, 1] методом Фибоначчи с абсолютной точностью epsXabs
// у функции (*obj)(x), при obj != NULL или (*obj2)(x) при obj2 != NULL;
// Результат помещается по ссылке в структуру p:
static	void getLocalMin(double x0, double x1, double epsXabs,
		             TestFunc* obj, TProblem* obj2, TGlobalSolution & p);
static	double getRandomIn(double a, double b)
	{	double per = (double)rand() / (double)RAND_MAX;
		double r = a + per * (b - a);
		return r;
	}
static void generateFib()
	{
		if (!Fib.size())
		{
			Fib.push_back(1); Fib.push_back(1);
			for (int i = 2; i < 46; ++i)
				Fib.push_back(Fib[i - 2] + Fib[i - 1]);
		}
	}
 };

class PinterTask : public TestFunc {
	double xGlobMin; //Безусловный глобальный минимум
public:
	static	double FA, FB;//Принятые границы диапазона для функций Пинтера
public:
	PinterTask(double aNew, double bNew, double a = PinterTask::FA, double b = PinterTask::FB)
		: TestFunc(aNew, bNew, a, b)
	{  // [a, b] - "Родные" границы изменения аргумента
		paramGeneration();
		scanTestFunction();
	}
virtual	void paramGeneration()
	{ xGlobMin = getRandomIn(Fa, Fb);
	 //@ cout << "xGlobMin = " << xGlobMin << " - real coordinate" << endl;
	}
	double calcTestFunc(double x)
	{	double PinterFunc;
		PinterFunc = 0.0025*(x - xGlobMin)*(x - xGlobMin) 
			+ sin(x - xGlobMin)*sin(x - xGlobMin) 
			+ sin((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin))
			 * sin((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin));
	    return PinterFunc;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaPinterFunc;
		derivaPinterFunc = 0.005*(x - xGlobMin)
			+ 2.0*sin(x - xGlobMin)*cos(x - xGlobMin)
			+ 2.0*sin((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin))
			* cos((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin))
			*(1.0+2.0*(x - xGlobMin));
		return derivaPinterFunc;
	}
};
class GibsonTask : public TestFunc {
	double parA[14], parB[14], par{ 0.0 };
public:
	static double FA, FB;//Принятые границы диапазона для функций Гибсона
public:
	GibsonTask(double aNew, double bNew, double a = GibsonTask::FA, double b = GibsonTask::FB)
		: TestFunc(aNew, bNew, a, b)
	{  // [a, b] - "Родные" границы изменения аргумента
		paramGeneration();
		scanTestFunction();
	}
	virtual	void paramGeneration()
	{ for (int i = 0; i < 14; ++i)
		{ parA[i] = getRandomIn(0.0, 1.0);
		  parB[i] = getRandomIn(0.0, 1.0);
		}
	}
	double calcTestFunc(double x)
	{	double GibsonFunc = parB[0];
		for (int i = 1; i < 14; i++)
			GibsonFunc += parA[i] * sin(4 * M_PI*x*i) + parB[i] * cos(4 * M_PI*x*i);
		return GibsonFunc;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaGibsonFunc{ 0.0 };
		for (int i = 1; i < 14; i++)
			derivaGibsonFunc += (4.0 * M_PI*i)*(parA[i] * cos(4 * M_PI*x*i)
			    - parB[i] * sin(4 * M_PI*x*i));
		return derivaGibsonFunc;
	}
};

class SheckelTask : public TestFunc {
	double parC[10], parK[10], parA[10];
public:
	static	double FA, FB;//Принятые границы диапазона для функций Шекеля
public:
	SheckelTask(double aNew, double bNew, double a = SheckelTask::FA, double b = SheckelTask::FB)
		: TestFunc(aNew, bNew, a, b)
	{  // [a, b] - "Родные" границы изменения аргумента
		paramGeneration();
		scanTestFunction();
	}
	virtual	void paramGeneration()
	{
		for (int i = 0; i < 10; ++i)
		{
			parC[i] = getRandomIn(0.001, 2.0);
			parK[i] = getRandomIn(1.0, 9.0);
			parA[i] = getRandomIn(0.0, 10.0);
		}
	}
	double calcTestFunc(double x)
	{ double SheckelFunc{ 0.0 };

		for (int i = 0; i < 10; i++)
			SheckelFunc += - 1.0 / (parC[i] + parK[i]*(x - parA[i])*(x - parA[i]));
		return SheckelFunc;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaSheckelFunc{ 0.0 };
	  for (int i = 1; i < 10; i++)
	  { double tmp = parC[i] + parK[i] * (x - parA[i])*(x - parA[i]);
	    derivaSheckelFunc += 2.0*parK[i] * (x - parA[i]) / (tmp*tmp);
	  }
		return derivaSheckelFunc;
	}
};
class PolinomTask : public TestFunc {
	double parR[10], par{ 1.0 }, sign{ +1.0 };
	int rootCount;
public:
	static	double FA, FB;//Принятые границы диапазона для функций "Полином"
public:
	PolinomTask(double aNew, double bNew, int _rootCount, 
		        double a = PolinomTask::FA, double b = PolinomTask::FB)
		: TestFunc(aNew, bNew, a, b),rootCount(_rootCount)
	{  // [a, b] - "Родные" границы изменения аргумента
		paramGeneration();
		scanTestFunction();
	}
	virtual	void paramGeneration()
	{ for (int i = 1; i < rootCount - 1; ++i)
		{ parR[i] = getRandomIn(Fa + 0.05*(Fb - Fa), Fb - 0.05*(Fb - Fa));
	//@	  cout << "root = " << parR[i] << "; ";
		}
    //@	cout << endl;
	  parR[0] = Fa;  parR[rootCount - 1] = Fb;
	  sign = getRandomIn(-1.0, 1.0);
	  sign = sign < 0.0 ? -1.0 : 1.0; //Знаковый множитель
	  int tmpPar = (rand() % 3) - 1; //Случайное -1 или 0 или +1
	  switch (tmpPar)
	  { case -1: par = 0.1; break;
	    case 0:  par = 1.0; break;
		case +1: par = 10.0; break;
		default: par = 1.0;
	  }
	
	}
	double calcTestFunc(double x)
	{ double PolinomFunc{ 1.0 };
	  for (int i = 0; i < rootCount; i++)
			PolinomFunc = PolinomFunc *(x - parR[i]);
	  return sign * PolinomFunc * par;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaPolinomFunc{ 0.0 };
	  for (int j = 0; j < rootCount; ++j)
	  { double tmp;
	    tmp = 1.0;
		for (int i = 0; i < rootCount; ++i)
			if (i == j) continue;
			else tmp = tmp * (x - parR[i]);
		derivaPolinomFunc += tmp;
	  }
	 return sign * derivaPolinomFunc * par;
	}
};

class TProblem
{
	double Fa{ 0.0 }; double Fb{ 1.0 }; //Отрезок поиска
	TypeProblem typeProblem;//Тип задачи (из стандартного набора или тестового класса)
	int numberTask;//Номер задачи в стандартном наборе или номер задачи в тестовом классе
	int constraintsCount;//Количество ограничений
	vector<int> taskFunctions;//Номера функций, входящих в стандартную задачу
	TestFuncType typeMainConstraint, typeMinFunc;//Тип основ.ограничения и миним.функции
	vector<TestFunc*> testFunctions; //Набор функций тестовой задачи
	vector<double> constraintsValues;//Верхние ограничители функций тестовой задачи
	double alpha{ 0.98 }; //Относительный уровень отсечения всех ограничений
	double mu{ 0.25 }; //Желаемая доля допустимых точек
	int tmpPointCount{ 500 };//Размер сетки для оценивания меры допустимого множества

	int pointsCount{ 5000 };//Число точек в сетке при грубом определении Min и Max
	double epsX{ 1.0e-08 }; //Относит. точность определения глобальных минимумов по х
	double epsY{ 1.0e-05 }; //Относит. точность сравнения значений в глоб.минимумах
	vector<TGlobalSolution> globalSolutions; //Набор решений тестовой задачи
	double minFunc, maxFunc;//Границы колебания целевой функции на допустимом множестве
		//Временное хранилище грубых оценок условных локальных минимумов тестовой задачи
		// в нормированных координатах из [0, 1]:
	static vector<TGlobalSolution> localTmpSolutions;
	  //Текущий условный локальный минимум (который должен уточняться)
	  // в нормированных координатах из [0, 1]:
	TGlobalSolution anyLocalSolution;
public:
	TProblem(TypeProblem type = TypeProblem::AnyTask, int numTask = 1)
	{
		typeProblem = type;
		numberTask = numTask;
		if (type == TypeProblem::AnyTask)
		{  initializeTask();//Заполняет поля описания задачи
		}
		else
		{ // constraintsCount = constrCount;
		 //  generateTestTask();//Генерация тестовой задачи
		}
	}
	~TProblem()
	{	for (auto obj : testFunctions) delete obj;
	    testFunctions.clear();
		taskFunctions.clear();
		globalSolutions.clear();
		constraintsValues.clear();

	}
	TProblem(const TProblem &) = delete;
	TProblem(const TProblem &&) = delete;
	TProblem & operator=(const TProblem&) = delete;
	TProblem & operator=(const TProblem&&) = delete;

	//Установка параметров генерации
	void setParamsTestTask(TestFuncType _typeMinFunc,
		TestFuncType _typeMainConstraint, int constrCount = 2, double _mu = 0.25)
	{	typeMinFunc = _typeMinFunc; 
		typeMainConstraint = _typeMainConstraint;
		constraintsCount = constrCount;
		alpha = 0.98;
		mu = _mu;
	}
	//Оператор вычисления функции с аргументом x из нормированного отрезка [0,1]
	// numFunc = 0 - целевая функция
	double operator()(double x, int numFunc) 
	{
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > constraintsCount ? constraintsCount : numFunc;
		double xx = Fa + (Fb - Fa)*x;//Пересчет в реальные координаты
		if (typeProblem == TypeProblem::AnyTask)
		      return anyFunction(xx, taskFunctions[numFunc]);
		else return testFunction( xx, numFunc) - constraintsValues[numFunc];
	}
	//Оператор вычисления производной функции с аргументом x из нормированного отрезка [0,1]
    // numFunc = 0 - целевая функция
	double operator()(double x, int numFunc, char t)
	{	if (t != 'd') return 0.0;
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > constraintsCount ? constraintsCount : numFunc;
		double xx = Fa + (Fb - Fa)*x;//Пересчет в реальные координаты
		if (typeProblem == TypeProblem::AnyTask)
			return derivaAnyFunction(xx, taskFunctions[numFunc])*(Fb - Fa);
		else return derivaTestFunction(xx, numFunc)*(Fb - Fa);
	}
	//Вычисление функции обобщенного ограничения с аргументом x из нормир.отрезка [0, 1]
	double G(double x)
	{ double gmax = (*this)(x, 1), g;
	  for (int i = 2; i <= constraintsCount; ++i)
	   {  g = (*this)(x, i);
			if (g > gmax) gmax = g;
	   }
	  return gmax;
	}
	//Оператор для вычисления в нормированных координатах вспомогательной функции F(x)
	// при уточнении условного локального минимума (вызывается из scanTestFunction()):
	double operator()(double x)
	{ double g = G(x), F = max<double>(0.0, g);
		if (g <= 0) F += (*this)(x, 0);
		else F += max<double>((*this)(x, 0), anyLocalSolution.y);
		return F;
	}
	//Граници отрезка поиска:
	double left() const { return Fa; }
	double right() const { return Fb; }

	// Выч-е функции из встроенного набора задач 
	// в "родных" координатах задачи [Fa, Fb]:
	double anyFunction(double x, int numStFunc);
	// Вычисление функции сгенерированной тестовой задачи в координатах 
	// задачи [Fa, Fb], равных "родным" координатам минимизируемой функции:
	double testFunction(double x, int numFunc)
	{
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > (int)testFunctions.size()-1 ? (int)testFunctions.size() - 1 : numFunc;
		return (*(testFunctions[numFunc]))(x);
	}
	// Выч-е производной функции из встроенного набора задач
	// в "родных" координатах задачи [Fa, Fb]:
	double derivaAnyFunction(double x, int numStFunc);
	// Вычисление производной функции сгенерированной тестовой задачи в координатах
	// задачи [Fa, Fb], равных "родным" координатам минимизируемой функции:
	double derivaTestFunction(double x, int numFunc)
	{
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > (int)testFunctions.size() - 1 ? 
			                  (int)testFunctions.size() - 1 : numFunc;
		return (*(testFunctions[numFunc]))(x,'d');
	}
	void initializeTask(); //Инициализация встроенной одиночной задачи

	//Фабрика функций тестовой задачи:
	//Первый вызов должен быть при isConstraints == false
	// при isConstraints == false в поля Fa, Fb задачи должны помещаться значения
	// "родных" границ минимизируемой функции,
	// при последующих вызовах для генерации функций ограничений (при isConstraints == true)
	// конструкторы этих функций должны вызываться со значениями границ Fa, Fb,
	// взятых из полей генерируемой задачи
	TestFunc* generateTestFunction(TestFuncType funcType, bool isConstraints, int rootCount);

	void generateTestTask();//Генерация тестовой задачи
	int getConstraintsCount() const //Возвращает количество ограничений
	{  return constraintsCount;
	}
	//Возвращает значение верхнего ограничителя функции с номером numFunc для уровня alpha:
	double getConstraintValue (int numFunc, double alpha) const
	{
	 return testFunctions[numFunc]->getMin()
	  + (testFunctions[numFunc]->getMax() - testFunctions[numFunc]->getMin())*alpha;
	}
	//Обновление значения Alpha в экземпляре класса:
	void setNewAlphaValue(double newAlpha)
	{	alpha = newAlpha;
		for(int i = 1; i <= constraintsCount;  ++i)
			constraintsValues[i] = getConstraintValue(i, alpha);
	}
	//Приближенное определение доли допустимых точек на сетке из M точек:
	double getProportionGoodDomain(int M);
	//Подбор уровня отсечения Alpha для достижения меры mu допустимого множества
	// (вернет достигнутое приближение mu)
	double seachAlphaValue();
	//Узнаем найденный относительный уровень отсечения для функций ограничений:
	double getAlpha() const { return alpha; }
	//Узнаем желаемую относительную меру допустимого множества:
	double getMu() const { return mu; }
	//Определение максимального и минимального условных значений, 
	// формирование списка точно найденных условных глобальных минимумов тестовой задачи:
	void scanTestTask();
	double getMin() const { return minFunc; }
	double getMax() const { return maxFunc; }
	int getGlobalSolutionCount() //Количество глоб.минимумов тестовой задачи 
	{
		return (int)globalSolutions.size();
	}
	//Получение данных об i-м глобальном минимуме задачи:
	TGlobalSolution operator[](int i)
	{
		i = i < 0 ? 0 : i;
		i = (int)((size_t)i >= globalSolutions.size() ? globalSolutions.size() - 1 : i);
		return globalSolutions[i];
	}
	TypeProblem getTypeProblem() const { return typeProblem; }
	int getTaskNumber() const { return numberTask; }
	TestFuncType getTypeMinFunc() const { return typeMinFunc;}//Тип миним.функции
	TestFuncType getTypeMainConstraint() const { return typeMainConstraint;}//Тип основ.ограничения
	
};

