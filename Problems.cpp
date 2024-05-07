// +++ ВЕРСИЯ ГЕНЕРАТОРА ТЕСТОВЫХ ЗАДАЧ УСЛОВНОЙ ГЛОБАЛЬНОЙ ОПТИМИЗАЦИИ: v2 от 12.07.2019
#include "Problems.h"

//Временное хранилище грубых оценок локальных минимумов тестовой функции:
vector<TGlobalSolution> TestFunc::localTmpSolutions;
vector<unsigned int> TestFunc::Fib;
//Временное хранилище грубых оценок условных локальных минимумов тестовой задачи:
vector<TGlobalSolution> TProblem::localTmpSolutions;

//"Родные" диапазоны определения функций различных типов:
double PinterTask::FA = -5.0;
double PinterTask::FB = +5.0;
double GibsonTask::FA = 0.0;
double GibsonTask::FB = 1.0;
double SheckelTask::FA = 0.0;
double SheckelTask::FB = 10.0;
double PolinomTask::FA = 0.0;
double PolinomTask::FB = 1.0;

//Определение максимального и минимального значений, 
// формирование списка точно найденных глобальных минимумов тестовой функции:
void TestFunc::scanTestFunction()
{	localTmpSolutions.clear();
    globalSolutions.clear();
    double x, y, xNext, yNext, h = (FbNew - FaNew) / (pointsCount - 1);
	double yPred;
	 yPred = (*this)(FaNew);
	 minFunc = maxFunc = yPred;
	//@	cout << "-----------Set points:" << endl;
	//@	cout << "x = " << FaNew << " y = " << yPred << endl;
	x = FaNew + h;
	y = (*this)(x);
	if (y < minFunc) minFunc = y;
	if (y > maxFunc) maxFunc = y;
	//@	cout << "x = " << x << " y = " << y << endl;
	//Ищем все локальные "ямки"
	TGlobalSolution p;
	if (yPred < y)
	{	p.x = FaNew; p.y = yPred;
		localTmpSolutions.push_back(p);
	}
	for (int i = 2; i <= pointsCount - 1; i++)
	{ xNext = FaNew + h * i;
	  yNext = (*this)(xNext);
	  //@ cout << "x = " << xNext << " y = " << yNext << endl;
	  if (yNext < minFunc) minFunc = y;
	  if (yNext > maxFunc) maxFunc = y;
	  if (yPred > y && y < yNext)
	  {	  p.x = x; p.y = y;
		  localTmpSolutions.push_back(p);
	  }
	  if (y > yNext && i == pointsCount - 1)
	  {
		  p.x = xNext; p.y = yNext;
		  localTmpSolutions.push_back(p);
	  }
	  yPred = y; x = xNext; y = yNext;
	}
	//@	cout << "-----------Set local points:" << endl;
	//@	for (int j = 0; j < (int)localTmpSolutions.size(); ++j)
	//@		cout << "x = " << localTmpSolutions[j].x
	//@		<< " y = " << localTmpSolutions[j].y << endl;

	double porogFunc = minFunc + (maxFunc - minFunc) / 10;
	//Уточнение локальных минимумов, выделение глобальных:
	for (int i = 0; i < (int)localTmpSolutions.size(); ++i)
	{
		if (localTmpSolutions[i].y > porogFunc) continue;
		double x0, x1;
		x0 = localTmpSolutions[i].x - h;
		x1 = x0 + 2 * h;
		x0 = FaNew > x0 ? FaNew : x0;
		x1 = FbNew < x1 ? FbNew : x1;

		getLocalMin(x0, x1,(FbNew - FaNew)*epsX, this, nullptr, p);
		if (minFunc > p.y) minFunc = p.y;
		if (p.y < minFunc + epsY * (maxFunc - minFunc))
							globalSolutions.push_back( p );
	}
	int count = (int)globalSolutions.size(), minusCount = 0;
	for( int i =0; i < (int)globalSolutions.size() - minusCount; ++i)
		if (globalSolutions[i].y > minFunc + epsY * (maxFunc - minFunc))
		{   globalSolutions[i] = globalSolutions[count - 1 - minusCount];
	         minusCount++;
		}
	globalSolutions.resize(count - minusCount);
}
//Поиск минимума на [x0, x1] из [0, 1] методом Фибоначчи с абсолютной точностью epsXabs
// у функции (*obj)(x), при obj != NULL или (*obj2)(x) при obj2 != NULL;
// Результат помещается по ссылке в структуру p:
void TestFunc::getLocalMin(double x0, double x1, double epsXabs,
	TestFunc* obj, TProblem* obj2, TGlobalSolution & p)
{  	double delta = epsXabs;
    double eps = delta / 10000;
	int N = 2;
	while ((x1 - x0) / Fib[N] + eps > delta) ++N;
	double Lamda, a = x0, b = x1;
	Lamda = (double)Fib[N - 1] / (double)Fib[N];
	double y = a + Lamda * (b - a), x = b - Lamda * (b - a);
	double Fx, Fy;
		if (obj) { Fx = (*obj)(x), Fy = (*obj)(y); }
		else if (obj2) { Fx = (*obj2)(x), Fy = (*obj2)(y); }
	int stepsCount = 2;
	//@ cout << "a = " << a << "\t b = " << b << endl;
	while (stepsCount < N - 1)
	{ ++stepsCount;
	  Lamda = (double)Fib[N - stepsCount+1] / (double)Fib[N - stepsCount + 2];
		if (Fx < Fy)
		{	b = y;
			y = x; Fy = Fx;
			x = b - Lamda * (b - a);
			if (obj) Fx = (*obj)(x);
			else if (obj2) Fx = (*obj2)(x);
		}
		else
		{	a = x;
			x = y; Fx = Fy;
			y = a + Lamda * (b - a);
			if (obj) Fy = (*obj)(y);
			else if (obj2) Fy = (*obj2)(y);
		}
		//@ cout << "a = " << a << "\t b = " << b << endl;
	}
	// Пороводим N-е измерение:
	if (Fx < Fy)
	{	b = y;
		y = x; Fy = Fx;
		x = y - eps * (b - a);
		if (obj) Fx = (*obj)(x);
		else if (obj2) Fx = (*obj2)(x);
	}
	else
	{	a = x;
		x = y; Fx = Fy;
		y = x + eps * (b - a);
		if (obj) Fy = (*obj)(y);
		else if (obj2) Fy = (*obj2)(y);
	}
	//@ cout << "a = " << a << "\t b = " << b << endl;
	//Последняя обработка, формируем решение:
	if (Fx < Fy) { p.x = x;  p.y = Fx;}
	else { p.x = y;  p.y = Fy;}
	//@ cout << "p.x = " << p.x << " p.y = " << p.y <<endl;
}


// Выч-е функции из встроенного набора:
double TProblem::anyFunction(double x, int numStFunc)
{
	double y = 0.0;
	switch (numStFunc)
	{
	case 0:
		y = (x + 2.0) * x * (x - 2.0); break;
	case 1:
		y = (x - 0.5) * (x + 1.0); break;
	case 2:
		y = (x * x + 1.0) * x * x * (x - 1.0); break;
	case 3:
		y = x * (x - 1.0) * (x + 1.0); break;
	case 4:
		y = x * (x + 1.0) * (x + 2.0) * (x - 2.0); break;
	case 5:
		y = (x - 2) * (x + 1.0); break;
	case 6:
		y = (4 * x - 0.5) * (3 * x - 0.8) * (2 * x - 0.9) * (x - 1.0) + 0.7; break;
	case 7:
		y = (3 * x - 0.5) * (x - 0.8) * (6 * x - 1) + 0.7; break;
	case 8:
		y = x; break;
	default:
		y = x;
	}
	return y;
}
// Выч-е производной функции из встроенного набора:
double TProblem::derivaAnyFunction(double x, int numStFunc)
{
	double y = 0.0;
	switch (numStFunc)
	{
	case 0:
		y = x * (x - 2.0) + (x + 2.0) * (x - 2.0)
			+ (x + 2.0) * x; break;
	case 1:
		y = (x + 1.0) + (x - 0.5); break;
	case 2:
		y = 2 * x * x * x * (x - 1.0) + 2 * x * (x * x + 1.0) * (x - 1.0)
			+ (x * x + 1.0) * x * x; break;
	case 3:
		y = x * (x - 1.0) + (x + 1.0) * (x - 1.0)
			+ (x + 1.0) * x; break;
	case 4:
		y = (x + 1.0) * (x + 2.0) * (x - 2.0) + x * (x + 1.0) * (x + 2.0)
			+ x * (x + 2.0) * (x - 2.0) + x * (x + 1.0) * (x - 2.0); break;
	case 5:
		y = (x - 2) + (x + 1.0); break;
	case 6:
		y = 4 * (3 * x - 0.8) * (2 * x - 0.9) * (x - 1.0) + 3 * (4 * x - 0.5) * (2 * x - 0.9) * (x - 1.0) +
			2 * (4 * x - 0.5) * (3 * x - 0.8) * (x - 1.0) + (4 * x - 0.5) * (3 * x - 0.8) * (2 * x - 0.9);  break;
	case 7: 
		y = 3 * (x - 0.8) * (6 * x - 1) + (3 * x - 0.5) * (6 * x - 1) + 6 * (3 * x - 0.5) * (x - 0.8);  break;
	case 8:
		y = 1; break;
	default:
		y = 1;
	}
	return y;
}
//Инициализация встроенной задачи:
void TProblem::initializeTask() 
{
	taskFunctions.clear();
	switch (numberTask)
	{ case 0:
		taskFunctions = { 0,1 };
		Fa = -3; Fb = 2.0; break;
	case 1:
		taskFunctions = { 5,3,4 };
		Fa = -3; Fb = 5.0; break;
	case 2:
		taskFunctions = { 8, 6, 7 };
		Fa = 0; Fb = 1.0; break;
	default:
		taskFunctions = { 0,1 };
		Fa = -3.5; Fb = 2.5; break;

	}
	constraintsCount = (int)(taskFunctions.size() - 1);
}

//Фабрика функций тестовой задачи:
//Первый вызов должен быть при isConstraints == false
// при isConstraints == false в поля Fa, Fb задачи должны помещаться значения
// "родных" границ минимизируемой функции,
// при последующих вызовах для генерации функций ограничений (при isConstraints == true)
// конструкторы этих функций должны вызываться со значениями границ Fa, Fb,
// взятых из полей генерируемой задачи
TestFunc* TProblem::generateTestFunction(TestFuncType funcType, bool isConstraints, int rootCount)
{  switch (funcType)
	{
	case TestFuncType::Pinter:
		if (!isConstraints) { Fa = PinterTask::FA; Fb = PinterTask::FB; }
		return (TestFunc *) new PinterTask(Fa, Fb);
		
	case TestFuncType::Sheckel:
		if (!isConstraints) { Fa = SheckelTask::FA; Fb = SheckelTask::FB; }
		return (TestFunc *) new SheckelTask(Fa, Fb);
	
	case TestFuncType::Gibson:
		if (!isConstraints) { Fa = GibsonTask::FA; Fb = GibsonTask::FB; }
		return (TestFunc *) new GibsonTask(Fa, Fb);
		
	case TestFuncType::Polinom:
		if (!isConstraints) { Fa = PolinomTask::FA; Fb = PolinomTask::FB; }
		return (TestFunc *) new PolinomTask(Fa, Fb, rootCount);
	default:
		if (!isConstraints) { Fa = PolinomTask::FA; Fb = PolinomTask::FB; }
		return (TestFunc *) new PolinomTask(Fa, Fb, rootCount);
	}

}
void TProblem::generateTestTask()//Генерация тестовой задачи
{
	unsigned int startRnd =  10 * (int(typeMainConstraint)+1) 
		                   + 100 * (int(typeMinFunc)+1)
		+ 10000 * numberTask + 1000 * constraintsCount;
	srand(startRnd);

	testFunctions.clear();
	constraintsValues.clear();
	//Минимизируемая функция:
	testFunctions.push_back(generateTestFunction(typeMinFunc, false, 0));
	constraintsValues.push_back(0);
	////Первое ограничение:
	//testFunctions.push_back(generateTestFunction(TestFuncType::Polinom,true, 5));
	//        // Установка значения верхнего ограничителя
	//constraintsValues.push_back(getConstraintValue(1, alpha));
	      //Главное первое ограничение:
	testFunctions.push_back(generateTestFunction(typeMainConstraint, true, 5));
	         // Установка значения верхнего ограничителя
	constraintsValues.push_back(getConstraintValue(1, alpha));
	//Дополнительные ограничения (полиномиальные):
	for (int i = 2; i <= constraintsCount; ++i)
	{ 	//Дополнительное полиномиальное ограничение:
		int constrCount = 6 - i % 2;
		testFunctions.push_back(generateTestFunction(TestFuncType::Polinom, true, constrCount));
		       // Установка значения верхнего ограничителя
		constraintsValues.push_back(getConstraintValue(i, alpha));
	}
}
//Приближенное определение доли допустимых точек на сетке из M точек:
double TProblem::getProportionGoodDomain(int M)
{ int goodPointCount = 0;
	for (int i = 0; i <= M; ++i)
	{
		double x = i * 1.0 / (double)M, y;
		bool dop = true;
		for (int j = 0; j <= getConstraintsCount(); ++j)
		{	y = (*this)(x, j);
			if (j > 0 && y > 0) dop = false;
		}
		if (dop) ++goodPointCount;
	}
	return (double)goodPointCount / (double)M;
}
//Подбор уровня отсечения Alpha для достижения меры mu допустимого множества
// (вернет достигнутое приближение mu)
double TProblem::seachAlphaValue()
{
	double alphR = 1.0, alphL = 1.0, alphStep = 0.01, alphC;
	double muR =1.0, muL, muC;
	do{ 
		alphL -= alphStep;
		setNewAlphaValue(alphL);
		muL = getProportionGoodDomain( tmpPointCount );
		alphStep *= 1.2;
	} while (muL > mu);
	do {
		alphC = 0.5*(alphL + alphR);
		setNewAlphaValue(alphC);
		muC = getProportionGoodDomain(tmpPointCount);
		if (muC < mu) { alphL = alphC; muL = muC;}
		else if (muC == mu) { alpha = alphC; return mu; }
		else { alphR = alphC;   muR = muC; }
	} while (fabs(muC - mu) > 2.0 / tmpPointCount);
	if ((muR - mu) < (mu - muL)) { alpha = alphR; return muR; }
	else { alpha = alphL; return muL; }
}

//Определение максимального и минимального условных значений, 
// формирование списка точно найденных условных глобальных минимумов тестовой задачи.
// Поиск проводится в нормированных координатах, список глобальных минимумов заполняется
// в нормированных координатах:
void TProblem::scanTestTask()
{
	localTmpSolutions.clear();
	globalSolutions.clear();
	double FaNew = 0.0, FbNew = 1.0; //Границы нормированного отрезка
	double x, y, g, xNext, yNext, gNext, h = (FbNew - FaNew) / ((double)(pointsCount - 1));
	double yPred, gPred;
	double minX;
	bool goodPointExist = false; //Встретилась ли допустимая точка?
	yPred = (*this)(FaNew, 0);
	gPred = G(FaNew);
	if ( gPred <=0 )
	{ minFunc = maxFunc = yPred;
	  minX = FaNew;
	  goodPointExist = true;
	}

	//@	cout << "-----------Set points:" << endl;
	//@	cout << "x = " << FaNew << " y = " << yPred << endl;
	x = FaNew + h;
	y = (*this)(x, 0);
	g = G(x);
	if ( g <= 0 )
		if ( ! goodPointExist)
		{	minFunc = maxFunc = yPred;
			minX = x;
			goodPointExist = true;
		}
		else
		{	if (minFunc > y) minX = x;
			minFunc = min<double>(minFunc, y);
			maxFunc = max<double>(maxFunc, y);
		}
	//@	cout << "x = " << x << " y = " << y << endl;
	//Ищем все локальные "ямки"
	TGlobalSolution p;
	if (gPred<=0.0 && ((yPred < y && g <=0) || g >0))
	{
		p.x = FaNew; p.y = yPred;
		localTmpSolutions.push_back(p);
	}
	for (int i = 2; i <= pointsCount - 1; i++)
	{
		xNext = FaNew + h * i;
		yNext = (*this)(xNext, 0);
		gNext = G(xNext);
		//@ cout << "x = " << xNext << " y = " << yNext << endl;
		if (gNext <= 0)
			if (!goodPointExist)
			{	minFunc = maxFunc = yNext;
				minX = xNext;
				goodPointExist = true;
			}
			else
			{	if (minFunc > yNext) minX = xNext;
				minFunc = min<double>(minFunc, yNext);
				maxFunc = max<double>(maxFunc, yNext);
			}
		if (g <= 0.0)
		  if (  ( gPred <= 0.0 && yPred > y || gPred > 0.0 ) 
			  && (y < yNext && gNext <= 0.0 || gNext > 0.0 ) )
			{
				p.x = x; p.y = y;
				localTmpSolutions.push_back(p);
			}

		 if (gNext <= 0.0)
			if ( (y > yNext && g <= 0.0 || g > 0.0 )&& i == pointsCount - 1)
			{
				p.x = xNext; p.y = yNext;
				localTmpSolutions.push_back(p);
			}
		 yPred = y; gPred = g; x = xNext; y = yNext; g = gNext;
	}
	//@	cout << "-----------Set local points:" << endl;
	//@	for (int j = 0; j < (int)localTmpSolutions.size(); ++j)
	//@		cout << "x = " << localTmpSolutions[j].x
	//@		<< " y = " << localTmpSolutions[j].y << endl;

	double porogFunc = minFunc + (maxFunc - minFunc) / 10;
	//Уточнение локальных минимумов, выделение глобальных:
	for (int i = 0; i < (int)localTmpSolutions.size(); ++i)
	{
		if (localTmpSolutions[i].y > porogFunc) continue;
		double x0, x1;
		x0 = localTmpSolutions[i].x - h;
		x1 = x0 + 2 * h;
		x0 = max<double>( FaNew, x0);
		x1 = min<double>( FbNew, x1);
		anyLocalSolution = localTmpSolutions[i];
      //Уточняем условный локальный минимум методом Фибоначчи:
		TestFunc::getLocalMin(x0, x1, (FbNew - FaNew)*epsX, nullptr, this, p);

		if (minFunc > p.y) minFunc = p.y;
		if (p.y < minFunc + epsY * (maxFunc - minFunc))
			globalSolutions.push_back(p);
	}

	int count = (int)globalSolutions.size(), minusCount = 0;
	for (int i = 0; i < (int)globalSolutions.size() - minusCount; ++i)
		if (globalSolutions[i].y > minFunc + epsY * (maxFunc - minFunc))
		{
			globalSolutions[i] = globalSolutions[count - 1 - minusCount];
			minusCount++;
		}
	globalSolutions.resize(count - minusCount);
}